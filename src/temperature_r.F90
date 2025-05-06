!> Module adding evolution of refined temperature to AFiD
!! Temperature is evolved on a refined grid, and therefore
!! this module depends on the multiple-resolution and
!! interpolation modules
module afid_tempr
    use param
    use mgrd_arrays
    use afid_salinity, only: rays
    use afid_phasefield, only: pf_eps, read_phase_field_params, pf_Tm, pf_S   
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
    use afid_moisture, only: beta_q       
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


!> Compute the explicit terms for the tempr evolution
!! and store the result in hsal
subroutine ExplicitTempr
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
    udyrq = dyqr/pect
    udzrq = dzqr/pect

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
                          vxr(kp,jc,ic)*(tempr(kp,jc,ic) + tempr(kc,jc,ic)) &
                        - vxr(kc,jc,ic)*2.d0*temprbp(1,jc,ic) &
                    )*udx3mr(kc)*0.5d0
                elseif (kc==nxmr) then
                    hsx = ( &
                          vxr(kp,jc,ic)*2.d0*temprtp(1,jc,ic) &
                        - vxr(kc,jc,ic)*(tempr(kc,jc,ic) + tempr(km,jc,ic)) &
                    )*udx3mr(kc)*0.5d0
                else
                    hsx = ( &
                          vxr(kp,jc,ic)*(tempr(kp,jc,ic) + tempr(kc,jc,ic)) &
                        - vxr(kc,jc,ic)*(tempr(kc,jc,ic) + tempr(km,jc,ic)) &
                    )*udx3mr(kc)*0.5d0
                end if

                ! y-advection d/dy(vy * S)
                hsy = ( &
                      vyr(kc,jp,ic)*(tempr(kc,jp,ic) + tempr(kc,jc,ic)) &
                    - vyr(kc,jc,ic)*(tempr(kc,jc,ic) + tempr(kc,jm,ic)) &
                )*udyr

                ! z-advection d/dz(vz * S)
                hsz = ( &
                      vzr(kc,jc,ip)*(tempr(kc,jc,ip) + tempr(kc,jc,ic)) &
                    - vzr(kc,jc,ic)*(tempr(kc,jc,ic) + tempr(kc,jc,im)) &
                )*udzr

                !! If using immersed boundary, enforce zero lateral gradient at interface
                if (IBM) then
                    !!! ADD THIS TO IBM MODULE, CALL AS SUBROUTINE
                    ! yy second derivative of salinity
                    if (solidr(kc,jp,ic)) then
                        dyys = (tempr(kc,jm,ic) - tempr(kc,jc,ic))*udyrq
                    elseif (solidr(kc,jm,ic)) then
                        dyys = (tempr(kc,jp,ic) - tempr(kc,jc,ic))*udyrq
                    else
                        dyys = (tempr(kc,jp,ic) - 2.0*tempr(kc,jc,ic) + tempr(kc,jm,ic))*udyrq
                    end if

                    ! zz second derivative of salinity
                    if (solidr(kc,jc,ip)) then
                        dzzs = (tempr(kc,jc,im) - tempr(kc,jc,ic))*udzrq
                    elseif (solidr(kc,jc,im)) then
                        dzzs = (tempr(kc,jc,ip) - tempr(kc,jc,ic))*udzrq
                    else
                        dzzs = (tempr(kc,jc,ip) - 2.0*tempr(kc,jc,ic) + tempr(kc,jc,im))*udzrq
                    end if
                else
                    ! yy second derivative of salinity
                    dyys = (tempr(kc,jp,ic) - 2.0*tempr(kc,jc,ic) + tempr(kc,jm,ic))*udyrq
                    ! zz second derivative of salinity
                    dzzs = (tempr(kc,jc,ip) - 2.0*tempr(kc,jc,ic) + tempr(kc,jc,im))*udzrq
                end if

                ! Sum explicit terms
                htempr(kc,jc,ic) = -(hsx + hsy + hsz) + dyys + dzzs
            end do
        end do
    end do

end subroutine ExplicitTempr

!> Compute the implicit terms for the tempr evolution
subroutine ImplicitTempr
    integer :: ic, jc, kc
    real :: dxxs, alpec

    alpec = al/pect

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr

                ! Second xx derivative
                ! Apply lower BC
                if (kc==1) then
                    dxxs= tempr(kc+1,jc,ic)*ap3ttkr(kc) &
                        + tempr(kc  ,jc,ic)*ac3ttkr(kc) &
                        - (ap3ttkr(kc) + ac3ttkr(kc))*temprbp(1,jc,ic)*TfixS
                ! Apply upper BC
                elseif (kc==nxmr) then
                    dxxs= tempr(kc  ,jc,ic)*ac3ttkr(kc) &
                        + tempr(kc-1,jc,ic)*am3ttkr(kc) &
                        - (am3ttkr(kc) + ac3ttkr(kc))*temprtp(1,jc,ic)*TfixN
                else
                    dxxs= tempr(kc+1,jc,ic)*ap3ttkr(kc) &
                        + tempr(kc  ,jc,ic)*ac3ttkr(kc) &
                        + tempr(kc-1,jc,ic)*am3ttkr(kc)
                end if

                rhsr(kc,jc,ic) = (ga*htempr(kc,jc,ic) + ro*rutempr(kc,jc,ic) + alpec*dxxs)*dt

                rutempr(kc,jc,ic) = htempr(kc,jc,ic)
            end do
        end do
    end do

    if (IBM) then
        !call SolveImpEqnUpdate_Sal_ibm
    else
        call SolveImpEqnUpdate_Tempr
    end if

end subroutine ImplicitTempr

!> Solve the implicit system for the salinity
!! and update the global variable sal
subroutine SolveImpEqnUpdate_Tempr
    real :: betadx, ackl_b
    integer :: ic, jc, kc, nrhs, ipkv(nxr), info
    real :: amkT(nxmr-1), ackT(nxmr), apkT(nxmr-1), appk(nxmr-2)

    betadx = 0.5d0*al*dt/pect

    ! Construct tridiagonal matrix for LHS
    do kc=1,nxmr
        ackl_b = 1.0d0/(1. - ac3ttkr(kc)*betadx)
        if (kc > 1) amkT(kc-1) = -am3ttkr(kc)*betadx*ackl_b
        ackT(kc) = 1.d0
        if (kc < nxmr) apkT(kc) = -ap3ttkr(kc)*betadx*ackl_b
    end do

    ! Factor the tridiagonal matrix
    call dgttrf(nxmr,amkT,ackT,apkT,appk,ipkv,info)

    ! Rescale RHS to match rescaling of LHS
    nrhs=(xendr(3)-xstartr(3)+1)*(xendr(2)-xstartr(2)+1)
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                ackl_b = 1.0/(1.0 - ac3ttkr(kc)*betadx)
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
             tempr(kc,jc,ic) = tempr(kc,jc,ic) + rhsr(kc,jc,ic)
            end do
         end do
     end do

end subroutine SolveImpEqnUpdate_Tempr

!> Interpolate the salinity field onto the coarse grid to
!! provide buoyancy forcing to the momentum equation
subroutine InterpTemprMultigrid
    integer :: icr, jcr, kcr

    ! Set coarse salinity array to zero
    tempc(:,:,:) = 0.d0

    ! Extend refined array in wall-normal direction to give sufficient points
    ! for cubic interpolation
    do icr=xstartr(3)-lvlhalo,xendr(3)+lvlhalo
        do jcr=xstartr(2)-lvlhalo,xendr(2)+lvlhalo
            do kcr=1,nxmr
                tpdvr(kcr,jcr,icr) = tempr(kcr,jcr,icr)
            end do
            if (TfixS==1) then
                tpdvr(0,jcr,icr) = 2.0*temprbp(1,jcr,icr) - tempr(1,jcr,icr)
            else
                tpdvr(0,jcr,icr) = tempr(1,jcr,icr)
            end if
            if (TfixN==1) then
                tpdvr(nxr,jcr,icr) = 2.0*temprtp(1,jcr,icr) - tempr(nxmr,jcr,icr)
            else
                tpdvr(nxr,jcr,icr) = tempr(nxmr,jcr,icr)
            end if
        end do
    end do

    ! Interpolate the refined field to the coarse grid, storing in salc
    if ((xmr(1) < xm(1)) .and. (xmr(nxmr) > xm(nxm))) then
        call interpolate_xyz_to_coarse_fast(tpdvr, tempc(1:nxm,:,:), "temp")
    else
        call interpolate_xyz_to_coarse(tpdvr, tempc(1:nxm,:,:))
    end if

end subroutine InterpTemprMultigrid

!> Add buoyancy contribution from tempr to one of the
!! momentum forcing arrays
subroutine AddTemprBuoyancy(rkv)
    real, dimension(:,xstart(2):,xstart(3):), intent(inout) :: rkv
    integer :: ic, jc, kc

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                rkv(kc,jc,ic) = rkv(kc,jc,ic) + byct*tempc(kc,jc,ic)
            end do
        end do
    end do
end subroutine

!> Add "latent tempr" term to the RK forcing array for salinity (hsal),
!! having calculated d/dt(phi) from the implicit solve and stored it in rhsr
subroutine AddLatentHeatr
    real :: aldt
    integer :: ic, jc, kc

    aldt = 1.0/al/dt

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                htempr(kc,jc,ic) = htempr(kc,jc,ic) &
                                + pf_S*rhsr(kc,jc,ic)*aldt
            end do
        end do
    end do
end subroutine AddLatentHeatr


!> Calculate and save vertical profiles related to the temperature
!! field, storing the data in means.h5
subroutine CalcTemprStats
    real, dimension(nxmr) :: Trbar   !! Horizontally-averaged temperature
    real, dimension(nxmr) :: Trrms   !! Horizontally-averaged rms temperature
    real, dimension(nxmr) :: chiTr   !! Horizontally-averaged temperature dissipation rate (nu grad(S)^2)
    
    real, dimension(nxmr) :: vxTr    !! Advective flux of salinity (x)
    real, dimension(nxmr) :: vyTr    !! Advective flux of salinity (y)
    real, dimension(nxmr) :: vzTr    !! Advective flux of salinity (z)

    real :: inyzmr      !! 1.0/nymr/nzmr

    character(30) :: dsetname   !! Dataset name for HDF5 file
    character(30) :: filename   !! HDF5 file name for statistic storage
    character( 5) :: nstat      !! Character string of statistic index

    integer :: i, j, k

    inyzmr = 1.0/nymr/nzmr

    filename = trim("outputdir/means.h5")

    Trbar(:) = 0.0;  Trrms(:) = 0.0;  chiTr(:) = 0.0
    vxTr(:) = 0.0;   vyTr(:) = 0.0;   vzTr(:) = 0.0

    if (IBM) then
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    ! Only record data from fluid phase
                    if (.not. solidr(k,j,i)) then
                        Trbar(k) = Trbar(k) + tempr(k,j,i)
                        Trrms(k) = Trrms(k) + tempr(k,j,i)**2
                    end if
                end do
            end do
        end do
    else
        do i=xstartr(3),xendr(3)
            do j=xstartr(2),xendr(2)
                do k=1,nxmr
                    Trbar(k) = Trbar(k) + tempr(k,j,i)
                    Trrms(k) = Trrms(k) + tempr(k,j,i)**2
                end do
            end do
        end do
    end if

    ! Since velocities are zero in solid, no need to use if statement for IBM here
    do i=xstartr(3),xendr(3)
        do j=xstartr(2),xendr(2)
            do k=1,nxmr
                vxTr(k) = vxTr(k) + 0.5*(vxr(k,j,i)+vxr(k+1,j,i))*tempr(k,j,i)
                vyTr(k) = vyTr(k) + 0.5*(vyr(k,j,i)+vyr(k,j+1,i))*tempr(k,j,i)
                vzTr(k) = vzTr(k) + 0.5*(vzr(k,j,i)+vzr(k,j,i+1))*tempr(k,j,i)
            end do
        end do
    end do

    call CalcDissipationTempr(chiTr)

    call MpiSumReal1D(Trbar, nxmr)
    call MpiSumReal1D(Trrms, nxmr)
    call MpiSumReal1D(vxTr,  nxmr)
    call MpiSumReal1D(vyTr,  nxmr)
    call MpiSumReal1D(vzTr,  nxmr)
    call MpiSumReal1D(chiTr, nxmr)

    ! Turn sums into averages
    ! (and root Trrms & scale chi)
    do k=1,nxmr
        Trbar(k) = Trbar(k)*inyzmr
        Trrms(k) = sqrt(Trrms(k)*inyzmr)
        vxTr(k)  = vxTr(k)*inyzmr
        vyTr(k)  = vyTr(k)*inyzmr
        vzTr(k)  = vzTr(k)*inyzmr
        chiTr(k) = chiTr(k)/pect*inyzmr
    end do

    ! Store index as character string
    write(nstat,"(i5.5)")nint(time/tout)

    if (ismaster) then
        dsetname = trim("Trbar/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, Trbar, nxmr)
        dsetname = trim("Trrms/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, Trrms, nxmr)
        dsetname = trim("vxTr/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vxTr, nxmr)
        dsetname = trim("vyTr/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vyTr, nxmr)
        dsetname = trim("vzTr/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, vzTr, nxmr)
        dsetname = trim("chiTr/"//nstat)
        call HdfSerialWriteReal1D(dsetname, filename, chiTr, nxmr)
    end if

    call MpiBarrier

end subroutine CalcTemprStats

!> Calculate dissipation rate for temperature (on local process, no MPI action here)
subroutine CalcDissipationTempr(chiTr)
    real, dimension(:), intent(out) :: chiTr

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
                            chiTr(k) = chiTr(k) + ((tempr(k,j,i ) - tempr(k,j,im))*0.5*dzr)**2
                        elseif (solidr(k,j,im)) then
                            chiTr(k) = chiTr(k) + ((tempr(k,j,ip) - tempr(k,j,i ))*0.5*dzr)**2
                        else
                            chiTr(k) = chiTr(k) + ((tempr(k,j,ip) - tempr(k,j,im))*0.5*dzr)**2
                        end if
                        if (solidr(k,jp,i)) then
                            chiTr(k) = chiTr(k) + ((tempr(k,j ,i) - tempr(k,jm,i))*0.5*dyr)**2
                        elseif (solidr(k,jm,i)) then
                            chiTr(k) = chiTr(k) + ((tempr(k,jp,i) - tempr(k,j ,i))*0.5*dyr)**2
                        else
                            chiTr(k) = chiTr(k) + ((tempr(k,jp,i) - tempr(k,jm,i))*0.5*dyr)**2
                        end if
                    end if
                end do
                if (.not. solidr(1,j,i)) then
                    chiTr(1) = chiTr(1) + (( &
                                tempr(2,j,i) - tempr(1,j,i) + 2.0*TfixS*(tempr(1,j,i)-temprbp(1,j,i))&
                            )*tdxr(1))**2
                end if
                do k=2,nxmr-1
                    if (.not. solidr(k,j,i)) then
                        chiTr(k) = chiTr(k) + ((tempr(k+1,j,i) - tempr(k-1,j,i))*tdxr(k))**2
                    end if
                end do
                if (.not. solidr(nxmr,j,i)) then
                    chiTr(nxmr) = chiTr(nxmr) + (( &
                                    tempr(nxmr,j,i) - tempr(nxmr-1,j,i) + 2.0*TfixN*(temprtp(1,j,i)-tempr(nxmr,j,i)) &
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
                    chiTr(k) = chiTr(k) + ((tempr(k,j,ip)-tempr(k,j,im))*0.5*dzr)**2
                    chiTr(k) = chiTr(k) + ((tempr(k,jp,i)-tempr(k,jm,i))*0.5*dyr)**2
                end do
                chiTr(1) = chiTr(1) + (( &
                            tempr(2,j,i) - tempr(1,j,i) + 2.0*TfixS*(tempr(1,j,i) - temprbp(1,j,i)) &
                        )*tdxr(1))**2
                do k=2,nxmr-1
                    chiTr(k) = chiTr(k) + ((tempr(k+1,j,i) - tempr(k-1,j,i))*tdxr(k))**2
                end do
                chiTr(nxmr) = chiTr(nxmr) + (( &
                            tempr(nxmr,j,i) - tempr(nxmr-1,j,i) + 2.0*TfixN*(temprtp(1,j,i)-tempr(nxmr,j,i)) &
                        )*tdxr(nxmr))**2
            end do
        end do
    end if
end subroutine CalcDissipationTempr

!> Create the groups in the means.h5 file to store the
!! salinity-related statistics
subroutine CreateTemprH5Groups(filename)
    use HDF5
    
    character(30), intent(in) :: filename
    integer(HID_T) :: file_id, group_id
    integer :: hdf_error

    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, hdf_error)

    call h5gcreate_f(file_id, "Trbar", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "Trrms", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vxTr", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vyTr", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "vzTr", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)
    call h5gcreate_f(file_id, "chiTr", group_id, hdf_error)
    call h5gclose_f(group_id, hdf_error)

    call h5fclose_f(file_id, hdf_error)

end subroutine CreateTemprH5Groups

end module afid_tempr
