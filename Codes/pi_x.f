C **********************************************************************
      FUNCTION qpi_x (pz,x,Lq,mpi,mq,mS,quas,typ,spin)

C
! .   .  PLOTS THE X DEPENDENCE OF THE INTEGRATED PION PDF AND QUASI-PDF
C
C       WE INTEGRATE NUMERICALLY OVER kT2; MOREOVER,
C       THE CODE CALLS THE EXTERNAL FUNCTION qpi_int
C
C  WRITTEN: T. Hobbs (JULY 26, 2017)
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ikT,nkT,typ,quas,spin
	PARAMETER(nkT=1000)
	EXTERNAL qpi_NORM,qpi_int
	REAL*8  qpi_int, qpi_NORM, qpi_x
	REAL*8  kT,kTmin,kTmax,kTint
	REAL*8  x,pz,Lq,mpi,mq,mS,Nq
	REAL*8  pix,pi_I
C***********************************************************************
! .   .   .   .   .   .  .
         Nq = qpi_NORM (pz,Lq,mpi,mq,mS,quas,typ,spin)
!          Nq = 1.D0

!---------------------------------------------------------------------
C***********************************************************************
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
!.................PROCEED WITH THE PDF CALCULATION.............


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
!
!.... FOR f_rho(x):
!              pi_I = qpi_int(pz,x,kT,Nq,Lq,mpi,mq,mS,quas,typ,spin)
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



!_   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _
        RETURN
!________________________________________________________________________________________
	END
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!________________________________________________________________________________________
