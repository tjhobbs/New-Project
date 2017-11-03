C **********************************************************************
         FUNCTION qpi_NORM (pz,Lq,mpi,mq,mS,quas,typ,spin)

! .   .  DETERMINE THE PION PDF NORMALIZATION AS AN EXTERNAL FUNCTION
C
C       WE INTEGRATE NUMERICALLY OVER kT2; MOREOVER,
C       THE CODE CALLS THE EXTERNAL FUNCTION qpi_int
C
C  WRITTEN: T. Hobbs (AUG 7, 2017)
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ix,nx,ikT,nkT,typ,quas,spin
	PARAMETER(nx=1000)
	PARAMETER(nkT=1000)
	EXTERNAL qpi_int
	REAL*8  qpi_int, qpi_NORM
	REAL*8  x,xmin,xmax,xint
	REAL*8  kT,kTmin,kTmax,kTint
	REAL*8  pz,Lq,mpi,mq,mS,Nq
	REAL*8  pix,pi_I,pi0
C***********************************************************************
!__________________________________________________________________
           Nq = 1.D0    ! SET TO UNITY HERE TO FIND THE OVERALL NORM.
!---------------------------------------------------------------------
C***********************************************************************
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
!.................PROCEED WITH THE PDF CALCULATION.............
          pi0 = 0.D0
! DEFINE THE BOUNDS IN x:
	  xmin = 0.D0
	     xmax = 0.99999D0                        !INTEGRATE OVER x FULLY!
	     xint = (xmax-xmin)/DBLE(nx)

         DO ix = 1, nx
	      x = xmin + xint * DBLE(ix)
 
! DEFINE THE BOUNDS OF THE kT INTEGRAL --- i.e., kT \in [0, \infty):
	     pix = 0.D0

	  kTmin = 0.D0
	     kTmax = 10.D0
	     kTint = (kTmax-kTmin)/DBLE(nkT)
	   DO ikT = 1, nkT
	      kT = kTmin + kTint * DBLE(ikT)
C________________________________________________________________________________________
!!.... FOR f_pi(x):
              pi_I = qpi_int(pz,x,kT,Nq,Lq,mpi,mq,mS,quas,typ,spin)
C  _    _   _   _    _    _    _    _    _    _    _    _    _    _
!.... FOR f_rho(x):
!              pi_I = qpi_int(pz,x,kT,Nq,Lq,mpi,mq,mS,0,typ,1)
C________________________________________________________________________________________
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
! WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD
!----- FIRST FOR THE kT INTEGRAL:
	     IF (iKT.EQ.0) THEN
	      pix = pix + pi_I
	    ELSE IF (iKT/2*2.NE.iKT) THEN
	      pix  = pix + 4.D0*pi_I
	    ELSE IF (iKT/2*2.EQ.iKT) THEN
	      pix  = pix + 2.D0*pi_I
	     ENDIF

	   ENDDO
	    pix = (kTint/3.D0) * pix
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!----- AND THEN DETERMINE NUMERICAL VALUES FOR THE x MOMENTS ALSO:
             IF (ix.EQ.0) THEN
              pi0 = pi0 + pix
            ELSE IF (ix/2*2.NE.ix) THEN
              pi0 = pi0 + 4.D0*pix
            ELSE IF (ix/2*2.EQ.ix) THEN
              pi0 = pi0 + 2.D0*pix
             ENDIF

          ENDDO
            pi0 = (xint/3.D0) * pi0

        qpi_NORM = 1.D0 / pi0

!________________________________________________________________________________________
        RETURN
	END
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!________________________________________________________________________________________
