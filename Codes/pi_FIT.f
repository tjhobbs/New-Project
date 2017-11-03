      PROGRAM pi_FIT

C     CREATED:  AUG 4, 2017
C     BY: T. HOBBS

C     THIS IS A SIMPLE CHI_SQUARED MINIMIZATION PROCEDURE TO COMPUTE THE
C     LINE OF BEST FIT TO QUARK DISTRIBUTIONS USING PHENOMENOLOGICAL MODELS

      implicit double precision (a-h,o-z)
      dimension nprm(2),stp(2),arglis(2)
      character*10 pnam(2)
      dimension vstrt(2)
      external fcn

      Nampl = 2 !THIS SPECIFIES THE NUMBER OF INDEPENDENT FITTING PARAMETERS


      OPEN (9,FILE='pi_GRVfit_1.dat',STATUS='UNKNOWN',FORM='FORMATTED')

      DATA NPRM  /   1    ,   2   /
      DATA PNAM  /  'pmS'  , 'pLAM' /
      DATA VSTRT /   0.4   , 0.7  /
      DATA STP   /  0.03  , 0.03  /

!     MINUIT IS INITIALIZED HERE
      call MNINIT(8,9,10)

      zero = 0.

      do 11 i = 1, 2

!     THIS IS THE SUBROUTINE THAT DEFINES/INVOKES THE VARIOUS FITTING
!     PARAMETERS WITH BOUNDS, INITIAL VALUES, ETC.
!     HERE, THE BOUNDS "zero, zero" MEAN UNBOUNDED

      call mnparm(nprm(i),pnam(i),vstrt(i),stp(i),zero,zero,ierflg)

!      call mnparm(nprm(1),pnam(1),vstrt(1),stp(1),0.2D0,1.0D0,ierflg)
!      call mnparm(nprm(2),pnam(2),vstrt(2),stp(2),0.01D0,3.0D0,ierflg)

      IF (ierflg .NE. 0) THEN
!         WRITE (9, '(A,i)') 'UNABLE TO DEFINE PARAMETER NO.', i
          PRINT*, "PROBLEM...."
         STOP
      ENDIF
   11 CONTINUE

!      THE CALLS TO "mnexcm" ARE RESPONSIBLE FOR ISSUING COMMANDS TO 
!      MINUIT.  HERE, 'MINIMIZE' IS THE GRADIENT/SIMPLEX-BASED
!      MINIMIZATION PROCEDURE 

       call mnseti(' PION SF PARAMETERS ')
       arglis(1) = 100000
!       call mnexcm(fcn,'MINIMIZE',arglis,1,ierflg,zero)
       call mnexcm(fcn,'SIMPLEX',arglis,1,ierflg,zero)

       arglis(1) = 9
       call mnexcm(fcn,'SET OUTPUTFILE',arglis,1,ierflg,zero)
       STOP
       end 

      subroutine fcn(npar,grad,fval,xval,iflag,futil)
      implicit double precision (a-h,o-z)
      integer npar,iflag,K
      external qpi_x
      dimension grad(*),xval(*),xp(8),pi(8)

      DATA (xp(K), K=1,8)
     &  /0.2875,0.3750,0.4625,0.5500,0.6375,0.7250,0.8125,0.9000/


      DATA (pi(K), K=1,8)
     &  /0.296299,0.330186,0.354116,0.368448,0.372530,0.364469,
     &   0.339998,0.287869/

!.  .  .  .   .  GRS DATA
!     & /0.28,0.36,0.44,0.52,0.60,0.68,0.76,0.84/    X-POINTS

!     & /0.28095995199473361,0.29517141040044409,0.30072021906187690,
!     &  0.29811420959385798,0.28130000000000005,0.26523642831998134,
!     &  0.23435620005727464,0.19029430190890367/
!     &  0.23435620005727464,0.19029430190890367,
!     &  5.5837236320339524E-002/                    DATA


      pmS = xval(1)
      pLAM = xval(2)

!TRADITIONALLY, 'fval' SHOULD BE A SIMPLE CHI-SQUARE; IN THE ABSENCE OF
!PUBLISHED UNCERTAINTIES, HOWEVER, I HAVE SET THE DENOMINATOR = 1. SO,
!WE JUST MINIMIZE ( theory - data ) ^ 2
      fval = 0.
      DO K = 1, 8
      x = xp(K)
      theory = x * qpi_x (1.D0,x,pLAM,0.139D0,0.33D0,pmS,0,1,0)
!      theory = x * qpi_x (1.D0,x,pLAM,0.139D0,0.33D0,pmS,0,0,0)
      
!      PRINT*, pmS, pLAM, theory

      fval = fval + (theory - pi(K))**2 / pi(K)**2
      ENDDO

      RETURN
      end 
