C ***********************************************************************
C HERE WE GATHER SEVERAL FUNCTIONS OF (x,Q2,kT) COMPUTED USING THE LFWF
C FORMALISM WITH SCALAR VERTICES; THESE MAY THEN BE CALLED TO COMPUTE,
C E.G., IC PDFs AND THE CHARM SIGMA TERM.
C
C     LAST EDITED: TJH, SEPT. 27th, 2016.
C ***********************************************************************
	FUNCTION F1int (Q2,x,kT,Nq,Lq,mq,mS,typ)
C
C  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE F1(Q2)
C      THIS IS A FUNCTION OF Q2, ... AND {x, kT} ARE TO BE INTEGRATED
C      OVER IN THE CALLING PROGRAM
C
C  WRITTEN: T. HOBBS (OCT. 2014)
C ***********************************************************************
	IMPLICIT NONE
	INTEGER typ
	REAL*8  Q2,x,kT,kT2,SInv,WF,tCOV
	REAL*8  F1int,Nq,Lq,mq,mS
	REAL*8  pi,mN
	REAL*8  kpkm_sum,X1f,X2f,X3f
	REAL*8  Aterm, Bterm

	pi = 4*DATAN(1.D0)
	mN = 0.9382720813D0   ! PROTON MASS; IN GeV!
C
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
           kT2 = kT**2
           SInv = ( kT2 + (1.D0-x)*mq**2 + x*mS**2
     &          + (1.D0-x)**2*Q2/4.D0 ) / (x * (1.D0-x) )
!
!. . . . THE COVARIANT PARAMETER t = k^2 . . . . . . . . . . . . . . . . .
           tCOV = -( kT2 + x * ( mS**2 - (1.D0-x)*mN**2 ) ) / (1.D0-x)

           kpkm_sum = 2.D0 * ( kT2 + (1.D0-x)**2*Q2/4.D0 )

           X1f = (1.D0/(4.D0*Lq**4)) * ( 1.D0/x**2 + 1.D0/(1.D0-x)**2
     &                                       + 2.D0/(x*(1.D0-x)) )

           X2f = (1.D0/(4.D0*Lq**4)) * ( (mq/x)**2 + (mS/(1.D0-x))**2
     &                                + (mq**2+mS**2)/(x*(1.D0-x)) )

           X3f = (1.D0/(4.D0*Lq**4)) * ( mq**4/x**2 + mS**4/(1.D0-x)**2
     &         + 2.D0*(mq**2*mS**2)/(x*(1.D0-x)) )

         Aterm = 1.D0 + SInv/Lq**2 + X1f*(kpkm_sum/2.D0)**2
     &         + X2f*kpkm_sum + X3f

         Bterm = (1.D0-x)**2 * kT2 * Q2 * X1f

!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
          IF (typ.EQ.1) THEN
            WF = DEXP( -SInv / Lq**2 )                                ! s-dep GAUSSIAN
          ELSE IF (typ.EQ.2) THEN
            WF = ( Lq**2 / (Lq**2 - tCOV) )**4                        ! COV. diPOLE
          ELSE IF (typ.EQ.3) THEN
            WF = (1.D0/2.D0) * (2.D0*Aterm - Bterm)/Aterm**3          ! s-dep DIPOLE
     &         * ( DSQRT( Aterm / (Aterm - Bterm) ) )**3
          ELSE IF (typ.EQ.4) THEN
            WF = ( Lq**2 / (Lq**2 - tCOV) )**2                        ! COV. monoPOLE
          ELSE IF (typ.EQ.5) THEN
            WF = ( Lq**2 / (Lq**2 - tCOV) )**6                        ! COV. quadPOLE
          ENDIF

        F1int = Nq * ( 1.D0/(16.D0*pi**2) )
     &             * ( kT2 + (mq + x*mN)**2 - (1.D0-x)**2*Q2/4.D0 )
     &             * ( 1.D0 / ( x**2 * (1.D0-x) ) ) * 2.D0*kT * WF 
     &             * ( 1.D0 / ( mN**2 - SInv )**2 )
	RETURN   
	END
C ***************************************************************************
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	FUNCTION F2int (Q2,x,kT,Nq,Lq,mq,mS,typ)
C
C  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE F2(Q2)
C      THIS IS A FUNCTION OF Q2, ... AND {x, kT} ARE TO BE INTEGRATED
C      OVER IN THE CALLING PROGRAM
C
C  WRITTEN: T. HOBBS (OCT. 2014)
C ***********************************************************************
	IMPLICIT NONE
	INTEGER typ
	REAL*8  Q2,x,kT,kT2,SInv,WF
	REAL*8  F2int,Nq,Lq,mq,mS
	REAL*8  pi,mN
	REAL*8  kpkm_sum,X1f,X2f,X3f
	REAL*8  Aterm, Bterm

	pi = 4*DATAN(1.D0)
	mN = 0.9382720813D0   ! KEEP ALL MASSES IN GeV!
C
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
           kT2 = kT**2
           SInv = ( kT2 + (1.D0-x)*mq**2 + x*mS**2
     &          + (1.D0-x)**2*Q2/4.D0 ) / (x * (1.D0-x) )

           kpkm_sum = 2.D0 * ( kT2 + (1.D0-x)**2*Q2/4.D0 )

           X1f = (1.D0/(4.D0*Lq**4)) * ( 1.D0/x**2 + 1.D0/(1.D0-x)**2
     &                                       + 2.D0/(x*(1.D0-x)) )

           X2f = (1.D0/(4.D0*Lq**4)) * ( (mq/x)**2 + (mS/(1.D0-x))**2
     &                                + (mq**2+mS**2)/(x*(1.D0-x)) )

           X3f = (1.D0/(4.D0*Lq**4)) * ( mq**4/x**2 + mS**4/(1.D0-x)**2
     &         + 2.D0*(mq**2*mS**2)/(x*(1.D0-x)) )

         Aterm = 1.D0 + SInv/Lq**2 + X1f*(kpkm_sum/2.D0)**2
     &         + X2f*kpkm_sum + X3f

         Bterm = (1.D0-x)**2 * kT2 * Q2 * X1f

!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
          IF (typ.EQ.1) THEN
            WF = DEXP( -SInv / Lq**2 )                                ! GAUSSIAN
          ELSE IF (typ.EQ.2) THEN
            WF = 0.D0                                                 ! MONOPOLE...
          ELSE IF (typ.EQ.3) THEN
            WF = (1.D0/2.D0) * (2.D0*Aterm - Bterm)/Aterm**3          ! DIPOLE
     &         * ( DSQRT( Aterm / (Aterm - Bterm) ) )**3
          ENDIF

        F2int = Nq * ( mN/(8.D0*pi**2) )
     &        * ( (mq + x*mN) / x**2) * 2.D0*kT * WF 
     &        * ( 1.D0 / ( mN**2 - SInv )**2 )
	RETURN   
	END
C ***************************************************************************
C ***************************************************************************
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
	FUNCTION GAint (Q2,x,kT,Nq,Lq,mq,mS,typ)
C
C  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE GA(Q2)
C      THIS IS A FUNCTION OF Q2, ... AND {x, kT} ARE TO BE INTEGRATED
C      OVER IN THE CALLING PROGRAM
C
C  WRITTEN: T. HOBBS (Jan. 2015)
C ***********************************************************************
	IMPLICIT NONE
	INTEGER typ
	REAL*8  Q2,x,kT,kT2,SInv,WF
	REAL*8  GAint,Nq,Lq,mq,mS
	REAL*8  pi,mN
	REAL*8  kpkm_sum,kpkm2,Mq4

	pi = 4*DATAN(1.D0)
	mN = 0.9382720813D0   ! KEEP ALL MASSES IN GeV!
C
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
           kT2 = kT**2
           SInv = ( kT2 + (1.D0-x)*mq**2 + x*mS**2
     &          + (1.D0-x)**2*Q2/4.D0 ) / (x * (1.D0-x) )

           kpkm_sum = 2.D0 * ( kT2 + (1.D0-x)**2*Q2/4.D0 )
           kpkm2 = ( kT2 - (1.D0-x)**2*Q2/4.D0 )**2

           Mq4 = kpkm2 * ( 1/x**2 + 1/(1.D0-x)**2 + 2.D0/(x*(1.D0-x)) ) !TERM I
     &         + kpkm_sum                                               !TERM II
     &  * ( (mq/x)**2 + (mS/(1.D0-x))**2 + (mq**2+mS**2)/(x*(1.D0-x)) )
     &         + mq**4/x**2 + mS**4/(1.D0-x)**2                         !TERM III
     &         + 2.D0*(mq**2*mS**2)/(x*(1.D0-x))

!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
          IF (typ.EQ.1) THEN
            WF = DEXP( -SInv / Lq**2 )                                ! GAUSSIAN
          ELSE IF (typ.EQ.2) THEN
            WF = 1.D0 / ( 1.D0 + SInv/Lq**2 + Mq4/(4.D0*Lq**4) )      ! MONOPOLE; DIV.!
          ELSE IF (typ.EQ.3) THEN
            WF = 1.D0 / ( 1.D0 + SInv/Lq**2 + Mq4/(4.D0*Lq**4) )**2   ! DIPOLE
          ENDIF

        GAint = Nq * ( 1.D0/(16.D0*pi**2) )  !A DIMENSIONLESS REDEF.!!!
     &             * ( -kT2 + (mq + x*mN)**2 + (1.D0-x)**2*Q2/4.D0 )
     &             * ( 1.D0 / ( x**2 * (1.D0-x) ) ) * 2.D0*kT * WF 
     &             * ( 1.D0 / ( mN**2 - SInv )**2 )
	RETURN   
	END
C ***************************************************************************
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
C ***********************************************************************
	FUNCTION CCint (Q2,x,kT,Nq,Lq,mq,mS,typ)
C
C  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE THE
C      SCALAR DENSITY <cbar c>
C
C  WRITTEN: T. HOBBS (FEB. 2015)
C  MODIFIED:         (SEP. 2016)
C ***********************************************************************
	IMPLICIT NONE
	INTEGER typ
	REAL*8  Q2,x,kT,kT2,SInv,WF
	REAL*8  CCint,Nq,Lq,mq,mS
	REAL*8  pi,mN
	REAL*8  kpkm_sum,kpkm2,Mq4

	pi = 4*DATAN(1.D0)
	mN = 0.9382720813D0   ! KEEP ALL MASSES IN GeV!
C
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
           kT2 = kT**2
           SInv = ( kT2 + (1.D0-x)*mq**2 + x*mS**2
     &          + (1.D0-x)**2*Q2/4.D0 ) / (x * (1.D0-x) )

           kpkm_sum = 2.D0 * ( kT2 + (1.D0-x)**2*Q2/4.D0 )
           kpkm2 = ( kT2 - (1.D0-x)**2*Q2/4.D0 )**2

           Mq4 = kpkm2 * ( 1/x**2 + 1/(1.D0-x)**2 + 2.D0/(x*(1.D0-x)) ) !TERM I
     &         + kpkm_sum                                               !TERM II
     &  * ( (mq/x)**2 + (mS/(1.D0-x))**2 + (mq**2+mS**2)/(x*(1.D0-x)) )
     &         + mq**4/x**2 + mS**4/(1.D0-x)**2                         !TERM III
     &         + 2.D0*(mq**2*mS**2)/(x*(1.D0-x))

!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
          IF (typ.EQ.1) THEN
            WF = DEXP( -SInv / Lq**2 )                                ! GAUSSIAN
          ELSE IF (typ.EQ.2) THEN
            WF = 1.D0 / ( 1.D0 + SInv/Lq**2 + Mq4/(4.D0*Lq**4) )      ! MONOPOLE; DIV.!
          ELSE IF (typ.EQ.3) THEN
            WF = 1.D0 / ( 1.D0 + SInv/Lq**2 + Mq4/(4.D0*Lq**4) )**2   ! DIPOLE
          ENDIF

        CCint = Nq * ( 1.D0/(16.D0*pi**2) )  !A DIMENSIONLESS REDEF.!!!
     &             * ( kT2 + (mq + x*mN)**2 - (1.D0-x)**2*Q2/4.D0 )
     &             * ( 1.D0 / ( x**3 * (1.D0-x) ) ) * 2.D0*kT * WF
     &             * ( mq / mN )
     &             * ( 1.D0 / ( mN**2 - SInv )**2 )
	RETURN   
	END
C ***************************************************************************
C ***************************************************************************
! . . . NOW FOR THE `ZERO MODE' PART OF THE SIGMA TERM . . . . . . . . . . . 
C ***************************************************************************
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
C ***********************************************************************
	FUNCTION CCZM (Q2,x,kT,Nq,Lq,mq,mS,typ)
C
C  SPECIFIES THE INTEGRAND OF THE x AND kT MOMENTS THAT DETERMINE THE
C      SCALAR DENSITY <cbar c>
C
C  WRITTEN: T. HOBBS (FEB. 2015)
C  MODIFIED:         (JUNE 2017)
C ***********************************************************************
	IMPLICIT NONE
	INTEGER typ
	REAL*8  Q2,x,kT,kT2,SInv,WF
	REAL*8  CCZM,Nq,Lq,mq,mS
	REAL*8  pi,mN
	REAL*8  kpkm_sum,kpkm2,Mq4

	pi = 4*DATAN(1.D0)
	mN = 0.9382720813D0   ! KEEP ALL MASSES IN GeV!
C
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
           kT2 = kT**2
           SInv = ( kT2 + (1.D0-x)*mq**2 + x*mS**2
     &          + (1.D0-x)**2*Q2/4.D0 ) / (x * (1.D0-x) )

           kpkm_sum = 2.D0 * ( kT2 + (1.D0-x)**2*Q2/4.D0 )
           kpkm2 = ( kT2 - (1.D0-x)**2*Q2/4.D0 )**2

           Mq4 = kpkm2 * ( 1/x**2 + 1/(1.D0-x)**2 + 2.D0/(x*(1.D0-x)) ) !TERM I
     &         + kpkm_sum                                               !TERM II
     &  * ( (mq/x)**2 + (mS/(1.D0-x))**2 + (mq**2+mS**2)/(x*(1.D0-x)) )
     &         + mq**4/x**2 + mS**4/(1.D0-x)**2                         !TERM III
     &         + 2.D0*(mq**2*mS**2)/(x*(1.D0-x))

!.....  CHOICES FOR THE FUNCTIONAL FORM OF THE CHARM WAVEFUNCTION....
          IF (typ.EQ.1) THEN
            WF = DEXP( -SInv / Lq**2 )                                ! GAUSSIAN
          ELSE IF (typ.EQ.2) THEN
            WF = 1.D0 / ( 1.D0 + SInv/Lq**2 + Mq4/(4.D0*Lq**4) )      ! MONOPOLE; DIV.!
          ELSE IF (typ.EQ.3) THEN
            WF = 1.D0 / ( 1.D0 + SInv/Lq**2 + Mq4/(4.D0*Lq**4) )**2   ! DIPOLE
          ENDIF

        CCZM = Nq * ( 1.D0/(16.D0*pi**2) )  !A DIMENSIONLESS REDEF.!!!
     &             * ( mq + x*mN ) / mN
     &             * ( 1.D0 / ( x**2 * (1.D0-x) ) ) * 2.D0*kT * WF
     &             * ( 1.D0 / ( mN**2 - SInv ) )
	RETURN   
	END
C ***************************************************************************
C ***************************************************************************

