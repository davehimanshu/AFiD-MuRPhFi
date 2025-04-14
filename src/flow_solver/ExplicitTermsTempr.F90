!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsTempr.F90                          !
!    CONTAINS: subroutine ExplicitTermsTempr               !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the temperature when temperature is in refined grid.!
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsTempr
    use param
    use mgrd_arrays, only: vxr, vyr, vzr, tempr, hror
    use decomp_2d, only: xstartr,xendr
    implicit none
    integer :: jc,kc,ic
    integer :: jm,jp,im,ip
    real    :: htx,hty,htz,udyr,udzr
    real    :: udzqr,udyqr
    real    :: dyyt,dzzt
    
    udzr=dzr*0.5d0
    udyr=dyr*0.5d0
    udzqr=dzqr/pect
    udyqr=dyqr/pect
    
    do ic=xstartr(3),xendr(3)
        im=ic-1
        ip=ic+1
        do jc=xstartr(2),xendr(2)
            jm=jc-1
            jp=jc+1
            do kc=1,nxmr

                ! ! x-advection d/dx(vx * T)
                if (kc==1) then
                    htx = ( &
                          vxr(kc+1,jc,ic)*(tempr(kc+1,jc,ic) + tempr(kc,jc,ic)) &
                        - vxr(kc  ,jc,ic)*2.0*tempbp(1,jc,ic) &
                    )*0.5*udx3mr(kc)
                elseif (kc==nxmr) then
                    htx = ( &
                          vxr(kc+1,jc,ic)*2.0*temptp(1,jc,ic) &
                        - vxr(kc  ,jc,ic)*(tempr(kc,jc,ic) + tempr(kc-1,jc,ic)) &
                    )*0.5*udx3mr(kc)
                else
                    htx = ( &
                          vxr(kc+1,jc,ic)*(tempr(kc+1,jc,ic) + tempr(kc  ,jc,ic)) &
                        - vxr(kc  ,jc,ic)*(tempr(kc  ,jc,ic) + tempr(kc-1,jc,ic)) &
                    )*0.5*udx3mr(kc)
                end if

                ! ! z-advection d/dx(vz * T)
                htz = ( &
                      vzr(kc,jc,ip)*(tempr(kc,jc,ip)+tempr(kc,jc,ic)) &
                    - vzr(kc,jc,ic)*(tempr(kc,jc,ic)+tempr(kc,jc,im)) &
                )*udzr
                
                ! ! y-advection d/dx(vy * T)
                hty=( &
                      vyr(kc,jp,ic)*(tempr(kc,jp,ic)+tempr(kc,jc,ic)) &
                    - vyr(kc,jc,ic)*(tempr(kc,jc,ic)+tempr(kc,jm,ic)) &
                )*udyr

                !   zz second derivatives of temp
                dzzt=(tempr(kc,jc,ip) - 2.0*tempr(kc,jc,ic) + tempr(kc,jc,im))*udzqr
                !   yy second derivatives of temp
                dyyt=(tempr(kc,jp,ic) - 2.0*tempr(kc,jc,ic) + tempr(kc,jm,ic))*udyqr

                hror(kc,jc,ic) = -(htx+hty+htz)+dyyt+dzzt
                !   hro(kc,jc,ic) = dyyt+dzzt
            end do
            
        end do
    end do
    !$OMP  END PARALLEL DO
    
    return
end subroutine ExplicitTermsTempr
