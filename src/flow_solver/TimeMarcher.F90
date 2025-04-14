!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: TimeMarcher.F90                                !
!    CONTAINS: subroutine TimeMarcher                     !
!                                                         ! 
!    PURPOSE: Main time integrating routine, which calls  !
!     other subroutines for calculating the Navier-Stokes !
!     equations and advancing velocity and temperature in !
!     time                                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine TimeMarcher
    use param
    use local_arrays
    use mgrd_arrays, only: vxr,vyr,vzr,tempr
    use afid_pressure
    use afid_salinity
    use afid_phasefield
    use mpih
    use decomp_2d
    use ibm_param, only: aldto
    use afid_moisture
    implicit none
    integer :: ns
    integer :: j,k,i

    beta = dt/ren*0.5d0

    do ns=1,nsst                                                 

!RO     Coefficients for time marching integration (alpha, gamma, rho)
        if(ntime.le.1) then
            aldto = alm(1)*dt
        else
            aldto = al*dt
        end if
        al = alm(ns)
        ga = gam(ns)
        ro = rom(ns)

        ! if (melt) call UpdateBCs

        ! iF(ANY(IsNaN(phi))) write(*,*)nrank,'NaN in PHI pre-explicit'
        call ExplicitTermsVX
        call ExplicitTermsVY
        call ExplicitTermsVZ
        if (multires_T) then
                call ExplicitTermsTempr
        else
                call ExplicitTermsTemp
        end if

        if (salinity) then
            call ExplicitSalinity
            
            ! If using salinity as an active scalar, add its buoyancy contribution
            ! to the relevant component of the momentum equation
            if (active_S==1) then
                if (gAxis==1) then
                    call AddSalBuoyancy(qcap)
                elseif (gAxis==2) then
                    call AddSalBuoyancy(dph)
                elseif (gAxis==3) then
                    call AddSalBuoyancy(dq)
                end if
            end if
        end if

        if (phasefield) then
           if (pfield_a) call ExplicitPhase
           call AddVolumePenalty
           if (salinity) then
                call AddSaltFluxInterface
                call AdjustMeltPoint
           end if
           if (pfield_a) call ImplicitPhase
            ! Add the latent heat and salt terms *after* computing the implicit solve for phi
           if (pfield_a .and. multires_T) call AddLatentHeatr
           if (pfield_a .and. .not. multires_T) call AddLatentHeat
           if (salinity) call AddLatentSalt
        end if

        if (moist) call ExplicitHumidity

        if (moist) call AddCondensation

        ! iF(ANY(IsNaN(vx))) write(*,*)nrank,'NaN in VX pre-implicit',ns
        ! iF(ANY(IsNaN(vy))) write(*,*)nrank,'NaN in VY pre-implicit',ns
        ! iF(ANY(IsNaN(vz))) write(*,*)nrank,'NaN in VZ pre-implicit',ns
        ! iF(ANY(IsNaN(temp))) write(*,*)nrank,'NaN in TEMP pre-implicit',ns

        if (phasefield) call update_halo(phi,lvlhalo)

        call ImplicitAndUpdateVX
        call ImplicitAndUpdateVY
        call ImplicitAndUpdateVZ
        if (multires_T) then
                call ImplicitAndUpdateTempr
        else
                call ImplicitAndUpdateTemp
        end if
        if (salinity) call ImplicitSalinity

        if (moist) call ImplicitHumidity

        ! iF(ANY(IsNaN(vx))) write(*,*)nrank,'NaN in VX post-implicit',ns
        ! iF(ANY(IsNaN(vy))) write(*,*)nrank,'NaN in VY post-implicit',ns
        ! iF(ANY(IsNaN(vz))) write(*,*)nrank,'NaN in VZ post-implicit',ns
        ! iF(ANY(IsNaN(temp))) write(*,*)nrank,'NaN in TEMP post-implicit',ns

        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)

        call CalcLocalDivergence
        call SolvePressureCorrection

!EP this copy can be avoided by changing transpose_x_to_y_real and
!transpose_y_to_x_real so these routines can handle arrays with
!halo. This copy is a defacto array temporary. Using inferred size
!arrays in the transpose calls results in 5 more of these, and more
!memory usage.  Time spent on this copy is 0.1% for 65^3 grid.

        do i=xstart(3),xend(3)
            do j=xstart(2),xend(2)
                do k=1,nxm
                    dphhalo(k,j,i) = dph(k,j,i)
                end do
            end do
        end do

        call update_halo(dphhalo,lvlhalo)

        call CorrectVelocity
        call CorrectPressure

        call update_halo(vx,lvlhalo)
        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)
        call update_halo(pr,lvlhalo)
        if (multires_T) then
                call update_halo(tempr,lvlhalo)
        else
                call update_halo(temp,lvlhalo)
        end if
        if (salinity) call update_halo(sal,lvlhalo)
        if (phasefield) call update_halo(phi,lvlhalo)
        if (moist) call update_halo(humid,lvlhalo)

        if (salinity .or. multires_T) then
            call InterpVelMgrd !Vel from base mesh to refined mesh
            call update_halo(vxr,lvlhalo)
            call update_halo(vyr,lvlhalo)
            call update_halo(vzr,lvlhalo)
        end if

        if (salinity) then
            call InterpSalMultigrid !Sal from refined mesh to base mesh
            call update_halo(salc,lvlhalo)
        end if

        if (multires_T) then
            if (phasefield .and. .not. pfield_a) then
               do i = xstartr(3), xendr(3)
                do j = xstartr(2), xendr(2)
                 do k = 1, nxmr
                        if (phi(k,j,i) > 0.1) then
                        tempr(k,j,i) = 0.0
                        end if
                 end do
                end do
               end do
            end if
            call update_halo(tempr,lvlhalo)    
            call interpTempMultigridr
            call update_halo(temp,lvlhalo)
        end if

        if (phasefield) then    
            call InterpPhiMultigrid
            call update_halo(phic,lvlhalo)
            if (.not. pfield_a) then
               do i = xstart(3), xend(3)
                do j = xstart(2), xend(2)
                 do k = 1, nxm
                        if (phic(k,j,i) > 0.1) then
                        temp(k,j,i) = 0.0
                        end if
                 end do
                end do
               end do
            end if
            call update_halo(temp,lvlhalo)
       
            if (.not. multires_T) call InterpTempMultigrid
            if (.not. multires_T) call update_halo(tempr,lvlhalo)
        end if

        if (moist) call UpdateSaturation

        ! iF(ANY(IsNaN(vx))) write(*,*)nrank,'NaN in VX',ns
        ! iF(ANY(IsNaN(vy))) write(*,*)nrank,'NaN in VY',ns
        ! iF(ANY(IsNaN(vz))) write(*,*)nrank,'NaN in VZ',ns
        ! iF(ANY(IsNaN(temp))) write(*,*)nrank,'NaN in PHI',ns

    end do


    return
end subroutine TimeMarcher
