!> Module adding evolution of refined temperature to AFiD
!! Temperature is evolved on a refined grid, and therefore
!! this module depends on the multiple-resolution and
!! interpolation modules
module afid_tempr
    use param
    use mgrd_arrays
    use afid_salinity, only: rays
    use afid_phasefield, only: pf_eps, read_phase_field_params, pf_Tm   
    use decomp_2d, only: xstart, xend, xstartr, xendr, update_halo
    use AuxiliaryRoutines
    use HermiteInterpolations, only: interpolate_xyz_to_coarse, interpolate_xyz_to_coarse_fast
    use ibm_param, only: solidr
    implicit none

    real, allocatable, dimension(:,:,:) :: rutempr    !! RK storage array for temperature (previous substep)
    real, allocatable, dimension(:,:,:) :: htempr     !! RK storage array for temperature
    real, allocatable, dimension(:,:,:) :: tempc     !! Interpolated temperature field on coarse grid

    real, allocatable, dimension(:,:,:) :: temprbp    !! temperature boundary value (lower plate)
    real, allocatable, dimension(:,:,:) :: temprtp    !! temeprature boundary value (upper plate)

    real, allocatable, dimension(:) :: ap3ttkr      !! Upper diagonal derivative coefficient for temperature
    real, allocatable, dimension(:) :: ac3ttkr      !! Diagonal derivative coefficient for temperature
    real, allocatable, dimension(:) :: am3ttkr      !! Lower diagonal derivative coefficient for temperature

contains

!> Subroutine to allocate memory for temperature-related variables
subroutine InitTemprVariables
    
    ! Boundary planes
    call AllocateReal3DArray(temprbp,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    call AllocateReal3DArray(temprtp,1,1,xstartr(2)-lvlhalo,xendr(2)+lvlhalo,xstartr(3)-lvlhalo,xendr(3)+lvlhalo)
    
    ! Runge-Kutta storage arrays (without ghost cells)
    call AllocateReal3DArray(rutempr,1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))
    call AllocateReal3DArray(htempr, 1,nxr,xstartr(2),xendr(2),xstartr(3),xendr(3))

    ! Coarse array
    call AllocateReal3DArray(tempc,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

    ! Second derivative coefficients
    call AllocateReal1DArray(ap3ttkr,1,nxr)
    call AllocateReal1DArray(ac3ttkr,1,nxr)
    call AllocateReal1DArray(am3ttkr,1,nxr)

end subroutine InitTemprVariables

!> Deallocate the variables used for evolving temperature
subroutine DeallocateTemprVariables

    ! Boundary planes
    call DestroyReal3DArray(temprbp)
    call DestroyReal3DArray(temprtp)

    call DestroyReal3DArray(rutempr)
    call DestroyReal3DArray(htempr)

    ! Coarse array
    call DestroyReal3DArray(tempc)

    ! Second derivative coefficients
    call DestroyReal1DArray(ap3ttkr)
    call DestroyReal1DArray(ac3ttkr)
    call DestroyReal1DArray(am3ttkr)

end subroutine DeallocateTemprVariables

!> Set the values for the boundary planes of salinity
subroutine SetTemprBCs
    integer :: i, j

    if (rayt>=0) then ! unstable T gradient
        if (inslwN==0) then !Single heated wall case
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    temprtp(1,j,i)=0.0
                    temprbp(1,j,i)=1.0
                end do
            end do
        else
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    temprtp(1,j,i)=-0.5d0
                    temprbp(1,j,i)=0.5d0
                end do
            end do
        end if
    else              ! stable T gradient
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                temprtp(1,j,i)=0.5d0
                temprbp(1,j,i)=-0.5d0
            end do
        end do
    end if
    
    if (phasefield) then
        if (rayt>=0) then
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    temprtp(1,j,i) = 0.d0
                    temprbp(1,j,i) = 1.d0
                end do
            end do
        else
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    temprtp(1,j,i) = 1.d0
                    temprbp(1,j,i) = 0.d0
                end do
            end do
        end if
    end if

    if (moist) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                temprbp(1,j,i) = 0.0
                temprtp(1,j,i) = beta_q - 1.0
            end do
        end do
    end if
    
    call update_halo(temprtp,lvlhalo)
    call update_halo(temprbp,lvlhalo)

end subroutine SetTemprBCs

!> Set initial conditions for temperature field
!! N.B. This can get overwritten by CreateInitialPhase if also using phase-field
subroutine CreateInitialTempr
    integer :: j,k,i,kmid
    real :: xxx,yyy,zzz,eps,varptb,amp
    real :: t0,Lambda,r, x0, A, B, alpha
    real, dimension(11) :: yh, zh
    
    if ((RayT < 0) .and. (RayS < 0)) then
        !CJH: Stratified shear layer + noise in centre
        eps = 1e-2
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    tempr(k,j,i) = tanh(xmr(k) - 0.5*alx3)
                    call random_number(varptb)
                    tempr(k,j,i) = tempr(k,j,i) + &
                            cosh(xmr(k) - 0.5*alx3)**(-2)*eps*(2.0*varptb - 1.0)
                end do
            end do
        end do
    else
        ! Assign linear temperature profile in the nodes k=1 to k=nxm
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    xxx = xmr(k)
                    tempr(k,j,i) = temprbp(1,j,i) + (temprtp(1,j,i) - temprbp(1,j,i))*xmr(k)/alx3
                end do
            end do
        end do

        ! Add noise in the temperature profile
        eps = 1e-3
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    call random_number(varptb)
                    if (abs(xmr(k)-0.5) + eps > 0.5) then
                        amp = 0.5 - abs(xmr(k)-0.5) ! CJH Prevent values of |T| exceeding 0.5
                        tempr(k,j,i) = tempr(k,j,i) + amp*(2.d0*varptb - 1.d0)
                    else
                    tempr(k,j,i) = tempr(k,j,i) + eps*(2.d0*varptb - 1.d0)
                    end if
                end do
            end do
        end do
    end if

    if (gAxis==3 .and. active_Tr==0) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2) ! Convergence test
                do k=1,nxmr
                    xxx = xmr(k) ! Linear profile + sin perturbation
                    tempr(k,j,i) = temprbp(1,j,i) + (temprtp(1,j,i) - temprbp(1,j,i))*xmr(k)/alx3
                    tempr(k,j,i) = tempr(k,j,i) + sin(2.0*pi*xxx/alx3) - sin(6.0*pi*xxx/alx3)
                end do
            end do
        end do
    end if

    if (gAxis==2 .and. inslwN==0) then  ! Ke et al comparison case
        t0 = 1.4195567
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    amp = 0.0
                    do kmid=0,7
                        amp = amp + sin(2.0**kmid * 2.0*pi*ymr(j)/ylen)
                    end do
                    amp = 1.0 + 1e-3*amp + 1e-3*sin(46.0*pi*zmr(i)/zlen)
                    tempr(k,j,i) = amp*erfc(xmr(k)/2*sqrt(pect/t0))
                end do
            end do
        end do
    end if

    if (IBM .and. dPdy/=0) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    tempr(k,j,i) = 0.0
                end do
            end do
        end do
    end if

    if (moist) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    tempr(k,j,i) = 0.0
                end do
            end do
        end do
    end if

    if (phasefield) then
        ! Most of this is now in `afid_phasefield` in the routine `CreateInitialPhase`

        if (pf_IC==3) then
            do i=xstartr(3),xendr(3)
                do j=xstartr(2),xendr(2)
                    do k=1,nxmr
                        xxx = xmr(k)
                        ! Piecewise linear base profile for Purseed et al
                        if (xxx < h0) then
                            tempr(k,j,i) = 1.0 - (1.0 - pf_Tm)*xxx/h0
                        else
                            tempr(k,j,i) = pf_Tm*(1.0 - xxx)/(1.0 - h0)
                        end if
                    end do
                end do
            end do
        end if

        if (salinity) then
            if (pf_IC==1) then
                call read_phase_field_params(A, B, alpha)
                t0 = 1e-3
                x0 = 0.8
!                h0 = x0 + 2*alpha*sqrt(t0)
                do i=xstartr(3),xendr(3)
                    do j=xstartr(2),xendr(2)
                        do k=1,nxmr
                            if (xmr(k) <= h0) then
                                tempr(k,j,i) = 1 - A*erfc((x0 - xmr(k))/sqrt(t0)/2.0)
                            else
                                tempr(k,j,i) = 1 - A*erfc(-alpha)
                            end if
                        end do
                    end do
                end do
            else if (pf_IC==2) then
                call read_phase_field_params(A, B, alpha)
                t0 = 1e-3
!                h0 = 0.1 - 2*alpha*sqrt(t0)
                eps = 5e-3
                do i=xstartr(3),xendr(3)
                    do j=xstartr(2),xendr(2)
                        do k=1,nxmr
                            call random_number(varptb)
                            if (abs(ymr(j) - ylen/2.0) <= h0) then
                                tempr(k,j,i) = 1.0 - A*erfc(-alpha)
                            else if (ymr(j) < ylen/2.0) then
                                tempr(k,j,i) = 1.0 - A*erfc((ylen/2.0 - h0 - ymr(j))/sqrt(t0)/2.0) &
                                                + eps*(2.d0*varptb - 1.d0)
                            else
                                tempr(k,j,i) = 1.0 - A*erfc((ymr(j) - ylen/2.0 - h0)/sqrt(t0)/2.0) &
                                + eps*(2.d0*varptb - 1.d0)
                            end if
                        end do
                    end do
                end do
            else if (pf_IC==3) then
                call read_phase_field_params(A, B, alpha)
                ! Scallop initial condition
                yh = [0.0, ylen/3, 2*ylen/3, ylen, &
                        ylen/6, ylen/2, 5*ylen/6, &
                        0.0, ylen/3, 2*ylen/3, ylen]
                zh(1:4) = 0.0
                zh(5:7) = zlen/2
                zh(8:11) = zlen
                x0 = 0.8
                amp = 0.9
                eps = 5e-3
                do i=xstartr(3),xendr(3)
                    do j=xstartr(2),xendr(2)
!                        h0 = 0.0
                        do k=1,11
!                            h0 = max(h0, x0 - amp*((ym(j) - yh(k))**2 + (zm(i) - zh(k))**2))
                        end do
                        do k=1,nxmr
                            call random_number(varptb)
                            if (xmr(k) <= h0) then
                                tempr(k,j,i) = 1.0 + eps*(2.d0*varptb - 1.d0)
                            else
                                tempr(k,j,i) = 1.0 - A*erfc(-alpha)
                            end if
                        end do
                    end do
                end do
            else
                kmid = nxmr/2
                do i=xstartr(3),xendr(3)
                    do j=xstartr(2),xendr(2)
                        do k=1,kmid
                            tempr(k,j,i) = 1.0
                        end do
                        do k=kmid+1,nxmr
                            tempr(k,j,i) = 0.0
                        end do
                    end do
                end do
            end if
        end if

    end if

    if (melt) then
        A = 1.08995
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    ! call random_number(varptb)
                    ! temp(k,j,i) = eps*(2.d0*varptb - 1.d0) * exp(-xm(k)/0.1)
                    tempr(k,j,i) = 1.0 - A*erfc(xmr(k)*sqrt(pect)/2.0)
                end do
            end do
        end do
    end if

end subroutine CreateInitialTempr


!> Compute the explicit terms for the salinity evolution
!! and store the result in hsal
subroutine ExplicitSalinity
    integer :: ic, jc, kc
    integer :: im, jm, km
    integer :: ip, jp, kp

    real :: udyr, udzr, udyrq, udzrq
    real :: aldt
    real, dimension(1:nxmr) :: sdx
    real :: hsx, hsy, hsz
    real :: dyys, dzzs

    ! Advection coefficients
    udyr = 0.5d0*dyr
    udzr = 0.5d0*dzr
    ! Diffusion coefficients
    udyrq = dyqr/pecs
    udzrq = dzqr/pecs

    ! x-advection coefficients
    do kc=1,nxmr
        sdx(kc) = 0.5*dxr/g3rmr(kc)
    end do
    ! Time advancing pre-factor
    aldt = 1.0/al/dt

    do ic=xstartr(3),xendr(3)
        im = ic - 1
        ip = ic + 1
        do jc=xstartr(2),xendr(2)
            jm = jc - 1
            jp = jc + 1
            do kc=1,nxmr
                km = kc - 1
                kp = kc + 1

                ! x-advection d/dx (vx * S)
                if (kc==1) then
                    hsx = ( &
                          vxr(kp,jc,ic)*(sal(kp,jc,ic) + sal(kc,jc,ic)) &
                        - vxr(kc,jc,ic)*2.d0*salbp(1,jc,ic) &
                    )*udx3mr(kc)*0.5d0
                elseif (kc==nxmr) then
                    hsx = ( &
                          vxr(kp,jc,ic)*2.d0*saltp(1,jc,ic) &
                        - vxr(kc,jc,ic)*(sal(kc,jc,ic) + sal(km,jc,ic)) &
                    )*udx3mr(kc)*0.5d0
                else
                    hsx = ( &
                          vxr(kp,jc,ic)*(sal(kp,jc,ic) + sal(kc,jc,ic)) &
                        - vxr(kc,jc,ic)*(sal(kc,jc,ic) + sal(km,jc,ic)) &
                    )*udx3mr(kc)*0.5d0
                end if

                ! y-advection d/dy(vy * S)
                hsy = ( &
                      vyr(kc,jp,ic)*(sal(kc,jp,ic) + sal(kc,jc,ic)) &
                    - vyr(kc,jc,ic)*(sal(kc,jc,ic) + sal(kc,jm,ic)) &
                )*udyr

                ! z-advection d/dz(vz * S)
                hsz = ( &
                      vzr(kc,jc,ip)*(sal(kc,jc,ip) + sal(kc,jc,ic)) &
                    - vzr(kc,jc,ic)*(sal(kc,jc,ic) + sal(kc,jc,im)) &
                )*udzr

                !! If using immersed boundary, enforce zero lateral gradient at interface
                if (IBM) then
                    !!! ADD THIS TO IBM MODULE, CALL AS SUBROUTINE
                    ! yy second derivative of salinity
                    if (solidr(kc,jp,ic)) then
                        dyys = (sal(kc,jm,ic) - sal(kc,jc,ic))*udyrq
                    elseif (solidr(kc,jm,ic)) then
                        dyys = (sal(kc,jp,ic) - sal(kc,jc,ic))*udyrq
                    else
                        dyys = (sal(kc,jp,ic) - 2.0*sal(kc,jc,ic) + sal(kc,jm,ic))*udyrq
                    end if

                    ! zz second derivative of salinity
                    if (solidr(kc,jc,ip)) then
                        dzzs = (sal(kc,jc,im) - sal(kc,jc,ic))*udzrq
                    elseif (solidr(kc,jc,im)) then
                        dzzs = (sal(kc,jc,ip) - sal(kc,jc,ic))*udzrq
                    else
                        dzzs = (sal(kc,jc,ip) - 2.0*sal(kc,jc,ic) + sal(kc,jc,im))*udzrq
                    end if
                else
                    ! yy second derivative of salinity
                    dyys = (sal(kc,jp,ic) - 2.0*sal(kc,jc,ic) + sal(kc,jm,ic))*udyrq
                    ! zz second derivative of salinity
                    dzzs = (sal(kc,jc,ip) - 2.0*sal(kc,jc,ic) + sal(kc,jc,im))*udzrq
                end if

                ! Sum explicit terms
                hsal(kc,jc,ic) = -(hsx + hsy + hsz) + dyys + dzzs
            end do
        end do
    end do

end subroutine ExplicitSalinity

!> Compute the implicit terms for the salinity evolution
subroutine ImplicitSalinity
    integer :: ic, jc, kc
    real :: dxxs, alpec

    alpec = al/pecs

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr

                ! Second xx derivative
                ! Apply lower BC
                if (kc==1) then
                    dxxs= sal(kc+1,jc,ic)*ap3sskr(kc) &
                        + sal(kc  ,jc,ic)*ac3sskr(kc) &
                        - (ap3sskr(kc) + ac3sskr(kc))*salbp(1,jc,ic)*SfixS
                ! Apply upper BC
                elseif (kc==nxmr) then
                    dxxs= sal(kc  ,jc,ic)*ac3sskr(kc) &
                        + sal(kc-1,jc,ic)*am3sskr(kc) &
                        - (am3sskr(kc) + ac3sskr(kc))*saltp(1,jc,ic)*SfixN
                else
                    dxxs= sal(kc+1,jc,ic)*ap3sskr(kc) &
                        + sal(kc  ,jc,ic)*ac3sskr(kc) &
                        + sal(kc-1,jc,ic)*am3sskr(kc)
                end if

                rhsr(kc,jc,ic) = (ga*hsal(kc,jc,ic) + ro*rusal(kc,jc,ic) + alpec*dxxs)*dt

                rusal(kc,jc,ic) = hsal(kc,jc,ic)
            end do
        end do
    end do

    if (IBM) then
        call SolveImpEqnUpdate_Sal_ibm
    else
        call SolveImpEqnUpdate_Sal
    end if

end subroutine ImplicitSalinity

!> Solve the implicit system for the salinity
!! and update the global variable sal
subroutine SolveImpEqnUpdate_Sal
    real :: betadx, ackl_b
    integer :: ic, jc, kc, nrhs, ipkv(nxr), info
    real :: amkT(nxmr-1), ackT(nxmr), apkT(nxmr-1), appk(nxmr-2)

    betadx = 0.5d0*al*dt/pecs

    ! Construct tridiagonal matrix for LHS
    do kc=1,nxmr
        ackl_b = 1.0d0/(1. - ac3sskr(kc)*betadx)
        if (kc > 1) amkT(kc-1) = -am3sskr(kc)*betadx*ackl_b
        ackT(kc) = 1.d0
        if (kc < nxmr) apkT(kc) = -ap3sskr(kc)*betadx*ackl_b
    end do

    ! Factor the tridiagonal matrix
    call dgttrf(nxmr,amkT,ackT,apkT,appk,ipkv,info)

    ! Rescale RHS to match rescaling of LHS
    nrhs=(xendr(3)-xstartr(3)+1)*(xendr(2)-xstartr(2)+1)
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                ackl_b = 1.0/(1.0 - ac3sskr(kc)*betadx)
                rhsr(kc,jc,ic) = rhsr(kc,jc,ic)*ackl_b
            end do
        end do
    end do

    ! Solve tridiagonal system
    call dgttrs('N',nxmr,nrhs,amkT,ackT,apkT,appk,ipkv,rhsr,nxmr,info)

    ! Update global variable
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
           do kc=1,nxmr
             sal(kc,jc,ic) = sal(kc,jc,ic) + rhsr(kc,jc,ic)
            end do
         end do
     end do

end subroutine SolveImpEqnUpdate_Sal

!> Interpolate the salinity field onto the coarse grid to
!! provide buoyancy forcing to the momentum equation
subroutine InterpSalMultigrid
    integer :: icr, jcr, kcr

    ! Set coarse salinity array to zero
    salc(:,:,:) = 0.d0

    ! Extend refined array in wall-normal direction to give sufficient points
    ! for cubic interpolation
    do icr=xstartr(3)-lvlhalo,xendr(3)+lvlhalo
        do jcr=xstartr(2)-lvlhalo,xendr(2)+lvlhalo
            do kcr=1,nxmr
                tpdvr(kcr,jcr,icr) = sal(kcr,jcr,icr)
            end do
            if (SfixS==1) then
                tpdvr(0,jcr,icr) = 2.0*salbp(1,jcr,icr) - sal(1,jcr,icr)
            else
                tpdvr(0,jcr,icr) = sal(1,jcr,icr)
            end if
            if (SfixN==1) then
                tpdvr(nxr,jcr,icr) = 2.0*saltp(1,jcr,icr) - sal(nxmr,jcr,icr)
            else
                tpdvr(nxr,jcr,icr) = sal(nxmr,jcr,icr)
            end if
        end do
    end do

    ! Interpolate the refined field to the coarse grid, storing in salc
    if ((xmr(1) < xm(1)) .and. (xmr(nxmr) > xm(nxm))) then
        call interpolate_xyz_to_coarse_fast(tpdvr, salc(1:nxm,:,:), "sal")
    else
        call interpolate_xyz_to_coarse(tpdvr, salc(1:nxm,:,:))
    end if

end subroutine InterpSalMultigrid

!> Add buoyancy contribution from the salinity to one of the
!! momentum forcing arrays
subroutine AddSalBuoyancy(rkv)
    real, dimension(:,xstart(2):,xstart(3):), intent(inout) :: rkv
    integer :: ic, jc, kc

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                rkv(kc,jc,ic) = rkv(kc,jc,ic) + bycs*salc(kc,jc,ic)
            end do
        end do
    end do
end subroutine

!> Calculate and save vertical profiles related to the salinity
!! field, storing the data in means.h5
subroutine CalcSalStats
    real, dimension(nxmr) :: Sbar   !! Horizontally-averaged salinity
    real, dimension(nxmr) :: Srms   !! Horizontally-averaged rms salinity
    real, dimension(nxmr) :: chiS   !! Horizontally-averaged solutal dissipation rate (nu grad(S)^2)
    
    real, dimension(nxmr) :: vxS    !! Advective flux of salinity (x)
    real, dimension(nxmr) :: vyS    !! Advective flux of salinity (y)
    real, dimension(nxmr) :: vzS    !! Advective flux of salinity (z)

    real :: inyzmr      !! 1.0/nymr/nzmr

    character(30) :: dsetname   !! Dataset name for HDF5 file
    character(30) :: filename   !! HDF5 file name for statistic storage
    character( 5) :: nstat      !! Character string of statistic index

    integer :: i, j, k

    inyzmr = 1.0/nymr/nzmr

    filename = trim("outputdir/means.h5")

    Sbar(:) = 0.0;  Srms(:) = 0.0;  chiS(:) = 0.0
    vxS(:) = 0.0;   vyS(:) = 0.0;   vzS(:) = 0.0

    if (IBM) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    ! Only record data from fluid phase
                    if (.not. solidr(k,j,i)) then
                        Sbar(k) = Sbar(k) + sal(k,j,i)
                        Srms(k) = Srms(k) + sal(k,j,i)**2
                    end if
                end do
            end do
        end do
    else
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    Sbar(k) = Sbar(k) + sal(k,j,i)
                    Srms(k) = Srms(k) + sal(k,j,i)**2
                end do
            end do
        end do
    end if

    ! Since velocities are zero in solid, no need to use if statement for IBM here
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                vxS(k) = vxS(k) + 0.5*(vxr(k,j,i)+vxr(k+1,j,i))*sal(k,j,i)
                vyS(k) = vyS(k) + 0.5*(vyr(k,j,i)+vyr(k,j+1,i))*sal(k,j,i)
                vzS(k) = vzS(k) + 0.5*(vzr(k,j,i)+vzr(k,j,i+1))*sal(k,j,i)
            end do
        end do
    end do

    call CalcDissipationSal(chiS)

    call MpiSumReal1D(Sbar, nxmr)
    call MpiSumReal1D(Srms, nxmr)
    call MpiSumReal1D(vxS,  nxmr)
    call MpiSumReal1D(vyS,  nxmr)
    call MpiSumReal1D(vzS,  nxmr)
    call MpiSumReal1D(chiS, nxmr)

    ! Turn sums into averages
    ! (and root Srms & scale chi)
    do k=1,nxmr
        Sbar(k) = Sbar(k)*inyzmr
        Srms(k) = sqrt(Srms(k)*inyzmr)
        vxS(k)  = vxS(k)*inyzmr
        vyS(k)  = vyS(k)*inyzmr
        vzS(k)  = vzS(k)*inyzmr
        chiS(k) = chiS(k)/pecs*inyzmr
    end do

    ! Store index as character string
    write(nstat,"(i5.5)")nint(time/tout)

    if (ismaster) then
        dsetname = trim("Sbar/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, Sbar, nxmr)
        dsetname = trim("Srms/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, Srms, nxmr)
        dsetname = trim("vxS/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vxS, nxmr)
        dsetname = trim("vyS/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vyS, nxmr)
        dsetname = trim("vzS/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vzS, nxmr)
        dsetname = trim("chiS/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, chiS, nxmr)
    end if

    call MpiBarrier

end subroutine CalcSalStats

!> Calculate dissipation rate for salinity (on local process, no MPI action here)
subroutine CalcDissipationSal(chiS)
    real, dimension(:), intent(out) :: chiS

    integer :: i, ip, im
    integer :: j, jp, jm
    integer :: k

    real, dimension(1:nxmr) :: tdxr

    do k=1,nxmr
        tdxr(k) = 0.5*dxr/g3rmr(k)
    end do

    ! If we are using the IBM, do not add contributions from the solid phase,
    ! and enforce zero gradient at the boundary points
    if (IBM) then
        do i=xstartr(3),xendr(3)
            ip = i + 1
            im = i - 1
            do j=xstartr(2),xendr(2)
                jp = j + 1
                jm = j - 1
                do k=1,nxmr
                    if (.not. solidr(k,j,i)) then
                        if (solidr(k,j,ip)) then
                            chiS(k) = chiS(k) + ((sal(k,j,i ) - sal(k,j,im))*0.5*dzr)**2
                        elseif (solidr(k,j,im)) then
                            chiS(k) = chiS(k) + ((sal(k,j,ip) - sal(k,j,i ))*0.5*dzr)**2
                        else
                            chiS(k) = chiS(k) + ((sal(k,j,ip) - sal(k,j,im))*0.5*dzr)**2
                        end if
                        if (solidr(k,jp,i)) then
                            chiS(k) = chiS(k) + ((sal(k,j ,i) - sal(k,jm,i))*0.5*dyr)**2
                        elseif (solidr(k,jm,i)) then
                            chiS(k) = chiS(k) + ((sal(k,jp,i) - sal(k,j ,i))*0.5*dyr)**2
                        else
                            chiS(k) = chiS(k) + ((sal(k,jp,i) - sal(k,jm,i))*0.5*dyr)**2
                        end if
                    end if
                end do
                if (.not. solidr(1,j,i)) then
                    chiS(1) = chiS(1) + (( &
                                sal(2,j,i) - sal(1,j,i) + 2.0*SfixS*(sal(1,j,i)-salbp(1,j,i))&
                            )*tdxr(1))**2
                end if
                do k=2,nxmr-1
                    if (.not. solidr(k,j,i)) then
                        chiS(k) = chiS(k) + ((sal(k+1,j,i) - sal(k-1,j,i))*tdxr(k))**2
                    end if
                end do
                if (.not. solidr(nxmr,j,i)) then
                    chiS(nxmr) = chiS(nxmr) + (( &
                                    sal(nxmr,j,i) - sal(nxmr-1,j,i) + 2.0*SfixN*(saltp(1,j,i)-sal(nxmr,j,i)) &
                                )*tdxr(nxmr))**2
                end if
            end do
        end do
    else
        do i=xstartr(3),xendr(3)
            ip = i + 1
            im = i - 1
            do j=xstartr(2),xendr(2)
                jp = j + 1
                jm = j - 1
                do k=1,nxmr
                    chiS(k) = chiS(k) + ((sal(k,j,ip)-sal(k,j,im))*0.5*dzr)**2
                    chiS(k) = chiS(k) + ((sal(k,jp,i)-sal(k,jm,i))*0.5*dyr)**2
                end do
                chiS(1) = chiS(1) + (( &
                            sal(2,j,i) - sal(1,j,i) + 2.0*SfixS*(sal(1,j,i) - salbp(1,j,i)) &
                        )*tdxr(1))**2
                do k=2,nxmr-1
                    chiS(k) = chiS(k) + ((sal(k+1,j,i) - sal(k-1,j,i))*tdxr(k))**2
                end do
                chiS(nxmr) = chiS(nxmr) + (( &
                            sal(nxmr,j,i) - sal(nxmr-1,j,i) + 2.0*SfixN*(saltp(1,j,i)-sal(nxmr,j,i)) &
                        )*tdxr(nxmr))**2
            end do
        end do
    end if
end subroutine CalcDissipationSal

!> Create the groups in the means.h5 file to store the
!! salinity-related statistics
subroutine CreateSalinityH5Groups(filename)
    use HDF5
    
    character(30), intent(in) :: filename
    integer(HID_T) :: file_id, group_id
    integer :: hdf_error

    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

    call h5gcreate_f(file_id, "Sbar", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "Srms", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vxS", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vyS", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vzS", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "chiS", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)

    call h5fclose_f(file_id, hdf_error)

end subroutine CreateSalinityH5Groups

end module afid_tempr
