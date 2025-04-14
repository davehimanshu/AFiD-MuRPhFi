!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_Temp.F90                     !
!    CONTAINS: subroutine SolveImpEqnUpdate_Temp          !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     temperature, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Tempr
    use param
    use mgrd_arrays, only : tempr,rhsr
    use decomp_2d, only: xstartr,xendr
    implicit none
    real, dimension(nxr) :: amkl,apkl,ackl
    integer :: jc,kc,info,ipkv(nxmr),ic,nrhs
    real :: betadx,ackl_b
    real :: amkT(nxmr-1),ackT(nxmr),apkT(nxmr-1),appk(nxmr-2)

!     Calculate the coefficients of the tridiagonal matrix
!     The coefficients are normalized to prevent floating
!     point errors.

    betadx=0.5d0*al*dt/pect

    do kc=1,nxmr
        ackl_b=1.0d0/(1.0d0-ac3ssk(kc)*betadx)
        amkl(kc)=-am3ssk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ssk(kc)*betadx*ackl_b
    end do

    amkT=amkl(2:nxmr)
    apkT=apkl(1:(nxmr-1))
    ackT=ackl(1:nxmr)

!     Call to LAPACK library to factor tridiagonal matrix.
!     No solving is done in this call.

    call dgttrf(nxmr,amkT,ackT,apkT,appk,ipkv,info)
    
    nrhs=(xendr(3)-xstartr(3)+1)*(xendr(2)-xstartr(2)+1)
    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                ackl_b=1.0/(1.0-ac3ssk(kc)*betadx)
                rhsr(kc,jc,ic)=rhsr(kc,jc,ic)*ackl_b
            end do
        end do
    end do
      
    call dgttrs('N',nxmr,nrhs,amkT,ackT,apkT,appk,ipkv,rhsr(1:nxmr,:,:),nxmr,info)

    do ic=xstartr(3),xendr(3)
        do jc=xstartr(2),xendr(2)
            do kc=1,nxmr
                tempr(kc,jc,ic)=tempr(kc,jc,ic) + rhsr(kc,jc,ic)
            end do
        end do
    end do

    return
end subroutine SolveImpEqnUpdate_Tempr
