!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateTemp.F90                      !
!    CONTAINS: subroutine ImplicitAndUpdateTemp           !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the temperature and call the implicit solver.       !
!     After this routine, the temperature has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateTempr
    use param
    use mgrd_arrays, only: tempr,hror,rutempr,rhsr
    use decomp_2d, only: xstartr,xendr
    use ibm_param
    implicit none
    integer :: jc,kc,ic
    real    :: alpec,dxxt

    alpec=al/pect

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,temp) &
!$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
!$OMP   SHARED(ga,ro,alpec,dt) &
!$OMP   SHARED(rhs,rutemp,hro) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dxxt)
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr

!   Calculate second derivative of temperature in the x-direction.
!   This is the only term calculated implicitly for temperature.
                if (kc.eq.1) then       !CJH Apply lower BC
                    dxxt = tempr(kc+1,jc,ic)*ap3ssk(kc) &
                        + tempr(kc,jc,ic)*ac3ssk(kc) &
                        - (ap3ssk(kc)+ac3ssk(kc))*tempbp(1,jc,ic)*TfixS
                elseif(kc.eq.nxm) then  !CJH Apply upper BC
                    dxxt = tempr(kc,jc,ic)*ac3ssk(kc) &
                        + tempr(kc-1,jc,ic)*am3ssk(kc) &
                        - (am3ssk(kc)+ac3ssk(kc))*temptp(1,jc,ic)*TfixN
                else
                    dxxt = tempr(kc+1,jc,ic)*ap3ssk(kc) &
                        + tempr(kc  ,jc,ic)*ac3ssk(kc) &
                        + tempr(kc-1,jc,ic)*am3ssk(kc)
                end if


!    Calculate right hand side of Eq. 5 (VO96)

                rhsr(kc,jc,ic)=(ga*hror(kc,jc,ic)+ro*rutempr(kc,jc,ic) &
                        +alpec*dxxt)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

                rutempr(kc,jc,ic)=hror(kc,jc,ic)

            enddo
        enddo
    enddo
!$OMP END PARALLEL DO

    if (IBM .and. .not. phasefield) then
        call SolveImpEqnUpdate_Temp_ibm
    else
!  Solve equation and update temperature

        call SolveImpEqnUpdate_Tempr
    end if

    return
end subroutine ImplicitAndUpdateTempr
