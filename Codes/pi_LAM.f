C **********************************************************************
	 PROGRAM pi_LAM
C
! .   .  PLOTS THE LAMBDA DEPENDENCE OF THE INTEGRATED PION PDF MOMENT
C
C       WE INTEGRATE NUMERICALLY OVER kT2; MOREOVER,
C       THE CODE CALLS THE EXTERNAL FUNCTION qpi_int
C
C  WRITTEN: T. Hobbs (JULY 26, 2017)
C **********************************************************************
C VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER il,nl,ix,nx,ikT,nkT,typ
	PARAMETER(nl=300)
	PARAMETER(nx=1000)
	PARAMETER(nkT=1000)
	EXTERNAL qpi_int
!        EXTERNAL GRSPI, qpi_x
        EXTERNAL GRV, qpi_x
	REAL*8  qpi_int, qpi_x, pi_FIT(nx), pmS, pLAM
	REAL*8  lmin,lmax,lint,Lqpt(nl),xpiL(nl)
	REAL*8  x,xmin,xmax,xint,xpt(nx)
	REAL*8  kT,kTmin,kTmax,kTint
	REAL*8  pz,Lq,mpi,mq,mS,Nq
	REAL*8  pi(nx),xpi(nx),piq(nx),disc(nx),pi_I,piq_I
	REAL*8  pi0, xpi0, piq0, GRSv(nx)
!        REAL*8  VAP,QBP,SBP,GLP
        REAL*8  xVpi,xSpi
C***********************************************************************
       pz = 1.D0                ! pion LMET BOOST MOMENTUM
       typ = 0                   ! USE THE COVAR. halfPOLE WF
!       typ = 1                   ! USE THE COVAR. monoPOLE WF
!       typ = 2                   ! USE THE COVAR. DIPOLE WF
!__________________________________________________________________
!_ _ _ _ VALUES OF THE FITTED PARAMETERS WE PLACE HERE _ _ _ _ _ _ _ _
! .   .   .   .   .   .  .
       Nq = 10.D0          ! COUPLING; (NORMED AWAY!)
!  -    -   TUNED TO GIVE LATTICE MOMENTS (BEST ET AL)    -    -    -
!       Lq = 0.0892D0      ! COVARIANT CUTOFF (FOR PION)     
!       mpi = 0.139D0      ! PION MASS
!       mq = 0.33D0        ! LIGHT CONSTIT. QUARK MASS
!       mS = 0.33D0        ! SPECTATOR QUARK MASS
!
!       Lq = 0.275D0        ! COVARIANT CUTOFF (FOR RHO)
       mpi = 0.77D0        ! RHO MASS
       mq = 0.6D0          ! EFF. LIGHT CONSTIT. QUARK MASS
       mS = 0.6D0          ! SPECTATOR QUARK MASS
!  -    -    -    -    -    -    -    -    -    -    -    -    -    -
!
!       mpi = 0.548D0      ! ETA MASS
!       mpi = 0.938D0      ! NUCLEON MASS
!       mS = 1.D0          ! SPECTATOR diQUARK MASS
!__________________________________________________________________

!---------------------------------------------------------------------
C***********************************************************************
!.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  
!.................PROCEED WITH THE PDF CALCULATION.............

	  lmin = 0.05D0
	     lmax = 0.5D0                        !INTEGRATE OVER x FULLY!
	     lint = (lmax-lmin)/DBLE(nl)

         DO il = 1, nl
	      Lq = lmin + lint * DBLE(il)
 
                 Lqpt(il) = Lq

          pi0 = 0.D0
          xpi0 = 0.D0
          piq0 = 0.D0
! DEFINE THE BOUNDS IN x:
	  xmin = 0.D0
	     xmax = 0.99999D0                        !INTEGRATE OVER x FULLY!
	     xint = (xmax-xmin)/DBLE(nx)

         DO ix = 1, nx
	      x = xmin + xint * DBLE(ix)
 
                 xpt(ix) = x

! DEFINE THE BOUNDS OF THE kT INTEGRAL --- i.e., kT \in [0, \infty):
	     pi(ix) = 0.D0
	     piq(ix) = 0.D0

	  kTmin = 0.D0
	     kTmax = 10.D0
	     kTint = (kTmax-kTmin)/DBLE(nkT)
	   DO ikT = 1, nkT
	      kT = kTmin + kTint * DBLE(ikT)
C________________________________________________________________________________________
!!.... FOR f_pi(x):
!              pi_I = qpi_int(pz,x,kT,Nq,Lq,mpi,mq,mS,0,typ,0)
!!
!!.... FOR \tilde{f}_pi(x):
!              piq_I = qpi_int(pz,x,kT,Nq,Lq,mpi,mq,mS,1,typ,0)
C  _    _   _   _    _    _    _    _    _    _    _    _    _    _
!.... FOR f_rho(x):
              pi_I = qpi_int(pz,x,kT,Nq,Lq,mpi,mq,mS,0,typ,1)
!
!.... FOR \tilde{f}_rho(x):
              piq_I = qpi_int(pz,x,kT,Nq,Lq,mpi,mq,mS,1,typ,1)
C________________________________________________________________________________________
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
! WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD
!----- FIRST FOR THE kT INTEGRAL:
	     IF (iKT.EQ.0) THEN
	      pi(ix) = pi(ix) + pi_I
	      piq(ix) = piq(ix) + piq_I
	    ELSE IF (iKT/2*2.NE.iKT) THEN
	      pi(ix)  = pi(ix) + 4.D0*pi_I
	      piq(ix) = piq(ix) + 4.D0*piq_I
	    ELSE IF (iKT/2*2.EQ.iKT) THEN
	      pi(ix)  = pi(ix) + 2.D0*pi_I
	      piq(ix) = piq(ix) + 2.D0*piq_I
	     ENDIF

	   ENDDO
	    pi(ix) = (kTint/3.D0) * pi(ix)
	    piq(ix) = (kTint/3.D0) * piq(ix)

           xpi(ix) = x * pi(ix)
!
!  .  .  . RENDER THIS RESULT AS A FRACTIONAL DISCREPANCY - - - 
!
              disc(ix) = 100.D0 * ( piq(ix) / pi(ix) - 1.D0 )
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .

!              CALL GRSPI(2,x,1.D0,VAP,QBP,SBP,GLP)
!              CALL GRSPI(2,x,0.4D0,VAP,QBP,SBP,GLP)
!                     GRSv(ix) = VAP

!              pmS = 2.1812D0       !'Half-pole'
!              pLAM = 1.5440D0

!              pmS = 0.61D0       !'mono-pole'
!              pLAM = 0.77875D0


!              pmS = 0.76000D0       !'Half-pole'
!              pLAM = 0.74500D0      ! GRV

!              pmS = 0.37150D0       !'mono-pole'
!              pLAM = 1.0230D0       ! GRV

!              pmS = 0.33D0

!         pi_FIT(ix) = x * qpi_x (1.D0,x,pLAM,0.139D0,0.33D0,pmS,0,0,0)
!         pi_FIT(ix) = x * qpi_x (1.D0,x,pLAM,0.139D0,0.33D0,pmS,0,1,0)
!
              CALL GRV(x,0.2501D0,xVpi,xSpi)
                     GRSv(ix) = xVpi

!                     PRINT*, x, GRSv(ix)
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!----- AND THEN DETERMINE NUMERICAL VALUES FOR THE x MOMENTS ALSO:
             IF (ix.EQ.0) THEN
              pi0 = pi0 + pi(ix)
              xpi0 = xpi0 + xpi(ix)
              piq0 = piq0 + piq(ix)
            ELSE IF (ix/2*2.NE.ix) THEN
              pi0 = pi0 + 4.D0*pi(ix)
              xpi0 = xpi0 + 4.D0*xpi(ix)
              piq0 = piq0 + 4.D0*piq(ix)
            ELSE IF (ix/2*2.EQ.ix) THEN
              pi0 = pi0 + 2.D0*pi(ix)
              xpi0 = xpi0 + 2.D0*xpi(ix)
              piq0 = piq0+ 2.D0*piq(ix)
             ENDIF

          ENDDO
            pi0 = (xint/3.D0) * pi0
            xpi0 = (xint/3.D0) * xpi0 / pi0
            piq0 = (xint/3.D0) * piq0


         xpiL(il) = xpi0
         PRINT*, il, Lq, xpiL(il)
         ENDDO

!         PRINT*, "THE TOTAL fpi AMPLITUDE:", pi0
!         PRINT*, "THE TOTAL xfpi 1st MOMENT:", xpi0
!         PRINT*, "THE TOTAL quasi-fpi AMPLITUDE ", piq0
!_   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _   _
C...WRITE DATA TO FILE
!  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
        OPEN (11,FILE='DATA/xrho_Lam.dat',
!        OPEN (11,FILE='DATA/fN_pz=5GeV.dat',
     &                           STATUS='UNKNOWN', FORM='FORMATTED')
          DO il=1,nl
          WRITE (11,*) Lqpt(il), xpiL(il)
        ENDDO
        CLOSE (11)
!
!________________________________________________________________________________________
!________________________________________________________________________________________
	END
!    .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
!________________________________________________________________________________________
