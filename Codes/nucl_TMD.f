C ***********************************************************************
C  A SIMPLE CODE FOR FOR EVALUATING THE k_perp DEPENDENT INTEGRANDS THAT
C  GO INTO THE MESON DISTRIBUTION FUNCTIONS; THE GOAL IS TO COMPARE
C  NUMERICALLY THE LIGHT-FRONT AND QUASI-PDFs
C
C     LAST EDITED: TJH, JULY 27th, 2017.
C
C      MODIFIED TO COMPUTE THE NUCLEON DISTRIBUTION(S)
C ***********************************************************************
        FUNCTION qpi_int (pz,x,kT,Nq,Lq,mpi,mq,mS,quas,typ)
C
C  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE fpi
C      THIS IS A FUNCTION OF pz, ... AND {x, kT} ARE TO BE INTEGRATED
C      OVER IN THE CALLING PROGRAM
C
C  WRITTEN: T. HOBBS (SUMMER.2017)
C ***********************************************************************
        IMPLICIT NONE
        INTEGER quas,typ
        REAL*8  pz,x,kT,kT2,SInv,WF,WFq
        REAL*8  tCOV,rho,lamb,delta,kbar2,numfa    
        REAL*8  qpi_int,Nq,Lq,mq,mS
        REAL*8  pi,mpi
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        pi = 4*DATAN(1.D0)
!
!        mpi = 0.938D0                  ! initial state NUCLEON MASS; IN GeV!
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
           numfa = (2.D0*x-1.D0)*mpi**2 - mS**2 + mq**2 + 2.D0*x*mpi*mq
     &           - 2.D0*(1.D0-x)**2*pz**2*( 1.D0 - delta * rho )


!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
          IF (typ.EQ.1) THEN
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
          IF (quas.EQ.0) THEN
!
        qpi_int = Nq * ( 1.D0/(16.D0*pi**2) )
     &             * ( kT2 + ( mq + x*mpi )**2 )
     &             * ( 1.D0 / ( x**2 * (1.D0-x) ) ) * 2.D0*kT * WF 
     &             * ( 1.D0 / ( mpi**2 - SInv )**2 )
!   .     .     .     .     .     .     .     .     .     .     .    .
          ELSE IF (quas.EQ.1) THEN
!
        qpi_int = Nq * ( 1.D0/(8.D0*pi**2) ) * 2.D0*kT * WFq
     &               * numfa
     &   * ( 1.D0 / ( 2.D0*(1.D0-x)*rho*(mq**2 - kbar2)**2 ) )
!
          ENDIF

	RETURN   
	END
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
