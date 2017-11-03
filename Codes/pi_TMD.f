C ***********************************************************************
C  A SIMPLE CODE FOR FOR EVALUATING THE k_perp DEPENDENT INTEGRANDS THAT
C  GO INTO THE MESON DISTRIBUTION FUNCTIONS; THE GOAL IS TO COMPARE
C  NUMERICALLY THE LIGHT-FRONT AND QUASI-PDFs
C
C     LAST EDITED: TJH, JULY 27th, 2017.
C ***********************************************************************
        FUNCTION qpi_int (pz,x,kT,Nq,Lq,mpi,mq,mS,quas,typ,spin)
C
C  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE fpi
C      THIS IS A FUNCTION OF pz, ... AND {x, kT} ARE TO BE INTEGRATED
C      OVER IN THE CALLING PROGRAM
C
C  WRITTEN: T. HOBBS (SUMMER.2017)
C ***********************************************************************
        IMPLICIT NONE
        INTEGER quas,typ,spin
        REAL*8  pz,x,kT,kT2,SInv,WF,WFq
        REAL*8  tCOV,rho,lamb,delta,kbar2    
        REAL*8  qpi_int,Nq,Lq,mq,mS
        REAL*8  pi,mpi,numfa1,numfa2,pbar0,kbar0,qbar0
        REAL*8  pdq,pdk,qdk
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        pi = 4*DATAN(1.D0)
!
!        mpi = 0.139D0                  ! initial state meson MASS; IN GeV!
!
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
           kT2 = kT**2

!  .   .   .   .   .   . FOR THE LIGHT-FRONT CALC :
           SInv = ( kT2 + (1.D0-x)*mq**2 + x*mS**2 )
     &          / (x * (1.D0-x) )

           tCOV = -( kT2 + x * ( mS**2 - (1.D0-x)*mpi**2 ) ) / (1.D0-x)

!  .   .   .   .   .   . FOR THE QUASI-pdf CALC :
           rho = DSQRT( 1.D0 + (kT2 + mS**2)/( (1.D0-x)**2 * pz**2) )
           lamb = (1.D0-x) * pz * rho
           delta = DSQRT( 1.D0 + mpi**2 / pz**2 )
           kbar2 = mpi**2 + mS**2 + 2.D0*(1.D0-x)*pz**2 
     &                                  *( 1.D0 - delta * rho )
!
!   .    .   . ENERGIES .    .    .    .    .     .   .
       pbar0 = pz*delta
       qbar0 = (1.D0-x)*pz*rho
       kbar0 = pbar0 - lamb
!   .    .   . INNER PRODUCTS .   .    .     .    .
       pdq = pbar0 * qbar0 - (1.D0-x) * pz**2 
       pdk = pbar0 * kbar0 - x * pz**2
       qdk = qbar0 * kbar0 + kT2 - x*(1.D0-x) * pz**2 
!
!   .  .  .  .  .  .  EXPLICITLY, THIS CAN BE EXPANDED AS :
!
!       kbar0 = DSQRT( mpi**2 + mS**2 + kT2 + (1.D0+(1.D0-x)**2)*pz**2 
!     &         - 2.D0*(1.D0-x)*pz**2*( 1.D0 - delta * rho ) )

           numfa1 = 2.D0*x * ( mq*mS + kT2 -x*(1.D0-x)*pz**2
     &                         + kbar0*qbar0 )
           numfa2 =  (1.D0-x) * (mq**2 - kbar2)

!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
          IF (typ.EQ.0) THEN
            WF = ( Lq**2 / (Lq**2 - tCOV) )               ! COV. halfPOLE
            WFq = ( Lq**2 / (Lq**2 - kbar2) )              !. . .quasi
          ELSE IF (typ.EQ.1) THEN
            WF = ( Lq**2 / (Lq**2 - tCOV) )**2               ! COV. monoPOLE
            WFq = ( Lq**2 / (Lq**2 - kbar2) )**2              !. . .quasi
          ELSE IF (typ.EQ.2) THEN
            WF = ( Lq**2 / (Lq**2 - tCOV) )**4               ! COV. diPOLE
            WFq = ( Lq**2 / (Lq**2 - kbar2) )**4              !. . .quasi
          ELSE IF (typ.EQ.3) THEN
            WF = ( Lq**2 / (Lq**2 - tCOV) )**6               ! COV. quadPOLE
            WFq = ( Lq**2 / (Lq**2 - kbar2) )**6              !. . .quasi
          ENDIF
!        
!  _    _     _     _   THE TMDS  --  _    _    _    _     _     _    _ 
      IF (spin.EQ.0) THEN            ! COMPUTE PS MESON VAL. PDF

          IF (quas.EQ.0) THEN
!
        qpi_int = Nq * ( 1.D0/(16.D0*pi**2) )
     &             * ( kT2 + ( (1.D0-x)*mq + x*mS )**2 )
     &             * ( 1.D0 / ( x**2 * (1.D0-x)**2 ) ) * 2.D0*kT * WF 
     &             * ( 1.D0 / ( mpi**2 - SInv )**2 )
!   .     .     .     .     .     .     .     .     .     .     .    .
          ELSE IF (quas.EQ.1) THEN   !  .     .    GIVE THE QUASI-PDF:
!
        qpi_int = Nq * ( 1.D0/(8.D0*pi**2) ) * 2.D0*kT * WFq
     &   * ( numfa1 + numfa2 )
     &   * ( 1.D0 / ( 2.D0*(1.D0-x)*rho*(mq**2 - kbar2)**2 ) )
!
          ENDIF
!  -    -    -    -    -    -    -    -   -    -    -    -     -   -
      ELSE IF (spin.EQ.1) THEN   ! COMPUTE VECTOR MESON VAL. PDF

          IF (quas.EQ.0) THEN
!
        qpi_int = Nq * ( 1.D0/(16.D0*pi**2) )
     &             * ( kT2 + 4.D0*x*(1.D0-x)*mq*mS
     &                  + ( (1.D0-x)*mq + x*mS )**2
     &             + (kT2 + mq**2 + x**2*mpi**2)
     &               * ( (1.D0-x)**2 + (kT2 + mS**2)/mpi**2 ) )
     &             * ( 1.D0 / ( x**2 * (1.D0-x)**2 ) ) * 2.D0*kT * WF 
     &             * ( 1.D0 / ( mpi**2 - SInv )**2 )
!   .     .     .     .     .     .     .     .     .     .     .    .
          ELSE IF (quas.EQ.1) THEN   !  .     .    GIVE THE QUASI-PDF:
!
        qpi_int = Nq * ( 1.D0/(8.D0*pi**2) ) * 2.D0*kT * WFq
     &   * ( x * ( 2.D0*pdq*pdk/mpi**2 + 3.D0*mq*mS + qdk )
     &       - ( kbar2 - mq**2 ) * ( pdq/mpi**2 + (1.D0-x)/2.D0 ) )   
     &   * ( 1.D0 / ( (1.D0-x)*rho*(mq**2 - kbar2)**2 ) )
!
          ENDIF
!
      ENDIF
!  .    .    .    .    .    .    .    .    .     .     .     .     .
	RETURN   
	END
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
