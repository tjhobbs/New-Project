c      SUBROUTINE MRSPI (X,SCALE,MODE,xuv,xdv,xqs,xst,xch,xbt,xgl)
      SUBROUTINE MRSPI (X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C***************************************************************C
C                                                               C
C    PION DISTRIBUTION FUNCTIONS EXTRACTED FROM DRELL-YAN AND   C
C    PROMPT PHOTON PION BEAM DATA, USING (NUCLEAR TARGET        C
C    CORRECTED) HMRS(B) DISTRIBUTIONS FOR THE NUCLEON.          C
C    MAIN UNCERTAINTY IS IN SEA, HENCE FOLLOWING THREE SETS...  C
C    MODE=2 IS THE "BEST FIT" SET.                              C
C                                                               C
C    MODE 1 CORRESPONDS TO A 10% MOMENTUM SEA    (0.11213)      C
C    MODE 2 CORRESPONDS TO A 15% MOMENTUM SEA    (0.16119)      C
C    MODE 3 CORRESPONDS TO A 20% MOMENTUM SEA    (0.23785)      C
C                                                               C
C    (NUMBER IN BRACKETS IS LAST ENTRY IN FIRST ROW OF GRID)    C
C                                                               C
C                         -*-                                   C
C                                                               C
C    REFERENCE: A.D. MARTIN. R.G. ROBERTS. W.J. STIRLING        C
C    AND P.J. SUTTON, PHYS. REV. D45 (1992) 2349.               C
C                         -*-                                   C
C                                                               C
C    FOR THESE SETS....                                         C
C    * LAMBDA(MSBAR,NF=4) = 190 MEV                             C
C    * UV=DV, STR=USEA=DSEA, BTM=0                              C
C                                                               C
C    NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION         C
C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C
C    IS THE QUARK SEA / 6   I.E. UBAR(X)=DBAR(X)                C
C    = SEA/X FOR A PION. IF IN DOUBT, CHECK THE                 C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV           C
C                                                               C
C                         -*-                                   C
C                                                               C
C     THE RANGE OF APPLICABILITY IS CURRENTLY:                  C
C     10**-5 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C
C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C
C     - PROBLEMS, COMMENTS ETC TO WJS@UK.AC.DUR.HEP             C
C                                                               C
C                                                               C
C***************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
c        real*4  xuv,xdv,xqs,xst,xch,xbt,xgl

	if (x.ge.1.d0) return
	if (scale.lt.dsqrt(5.d0)) scale=dsqrt(5.d0)

      IF(MODE.EQ.1) CALL PION31(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.2) CALL PION32(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.3) CALL PION33(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
c        xuv = SNGL(UPV)
c        xdv = SNGL(DNV)
c        xqs = SNGL(SEA)
c        xst = SNGL(STR)
c        xch = SNGL(CHM)
c        xbt = SNGL(BOT)
c        xgl = SNGL(GLU)
      RETURN
      END
C
      SUBROUTINE PION31(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C
C ::::::::::::  PION STRUCTURE FUNCTION :: 10% SEA :::::::::::::::::
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=50)
      PARAMETER(NQ=19)
      PARAMETER(NTENTH=21)
      DIMENSION F(7,NX,NQ),G(7),XX(NX),N0(7)
      DATA XX/1.D-5,2.D-5,4.D-5,6.D-5,8.D-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.85D0,.9D0,.95D0,.975D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/0,0,3,5,0,5,0/
      DATA INIT/0/
      XSAVE=X
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,NX-1
      DO 20 M=1,NQ
      open(101,file='/u/home/wmelnitc/Work/SF/DATA/mrs101.dat',
     & status='old',form='formatted')
      READ(101,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M)
C 1=UV 2=DV 3=GLUE 4=(UBAR+DBAR)/2 5=CBAR 7=BBAR 6=SBAR
         DO 25 I=1,7
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=2,6
      DO 30 K=1,NQ
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,NTENTH,K)/DLOG(F(I,NTENTH,K))
  50  FORMAT(7F10.5)
      DO 40 I=1,7
      DO 40 M=1,NQ
  40  F(I,NX,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,7
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.NTENTH) GOTO 65
      IF(I.EQ.7.OR.I.EQ.1) GOTO 65
          FAC=(1.D0-B)*F(I,NTENTH,M)+B*F(I,NTENTH,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(2) ! UPBAR DISTRIBUTION = D DISTRIBUTION
      DNV=G(2)
      SEA=G(4) ! THIS SEA IS (UBAR+DBAR)/2
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
      X=XSAVE
      RETURN
      END
C
      SUBROUTINE PION32(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C
C ::::::::::::  PION STRUCTURE FUNCTION :: 15% SEA :::::::::::::::::
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=50)
      PARAMETER(NQ=19)
      PARAMETER(NTENTH=21)
      DIMENSION F(7,NX,NQ),G(7),XX(NX),N0(7)
      DATA XX/1.D-5,2.D-5,4.D-5,6.D-5,8.D-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.85D0,.9D0,.95D0,.975D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/0,0,3,5,0,5,0/
      DATA INIT/0/
      XSAVE=X
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,NX-1
      DO 20 M=1,NQ
      open(102,file='/u/home/wmelnitc/Work/SF/DATA/mrs102.dat',
     & status='old',form='formatted')
      READ(102,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M)
C 1=UV 2=DV 3=GLUE 4=(UBAR+DBAR)/2 5=CBAR 7=BBAR 6=SBAR
         DO 25 I=1,7
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=2,6
      DO 30 K=1,NQ
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,NTENTH,K)/DLOG(F(I,NTENTH,K))
  50  FORMAT(7F10.5)
      DO 40 I=1,7
      DO 40 M=1,NQ
  40  F(I,NX,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,7
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.NTENTH) GOTO 65
      IF(I.EQ.7.OR.I.EQ.1) GOTO 65
          FAC=(1.D0-B)*F(I,NTENTH,M)+B*F(I,NTENTH,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(2) ! UPBAR DISTRIBUTION = D DISTRIBUTION
      DNV=G(2)
      SEA=G(4) ! THIS SEA IS (UBAR+DBAR)/2
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
      X=XSAVE
      RETURN
      END
C
C
      SUBROUTINE PION33(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C
C ::::::::::::  PION STRUCTURE FUNCTION :: 20% SEA :::::::::::::::::
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NX=50)
      PARAMETER(NQ=19)
      PARAMETER(NTENTH=21)
      DIMENSION F(7,NX,NQ),G(7),XX(NX),N0(7)
      DATA XX/1.D-5,2.D-5,4.D-5,6.D-5,8.D-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.85D0,.9D0,.95D0,.975D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/0,0,3,5,0,5,0/
      DATA INIT/0/
      XSAVE=X
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,NX-1
      DO 20 M=1,NQ
      open(103,file='/u/home/wmelnitc/Work/SF/DATA/mrs103.dat',
     & status='old',form='formatted')
      READ(103,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M)
C 1=UV 2=DV 3=GLUE 4=(UBAR+DBAR)/2 5=CBAR 7=BBAR 6=SBAR
         DO 25 I=1,7
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=2,6
      DO 30 K=1,NQ
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,NTENTH,K)/DLOG(F(I,NTENTH,K))
  50  FORMAT(7F10.5)
      DO 40 I=1,7
      DO 40 M=1,NQ
  40  F(I,NX,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,7
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.NTENTH) GOTO 65
      IF(I.EQ.7.OR.I.EQ.1) GOTO 65
          FAC=(1.D0-B)*F(I,NTENTH,M)+B*F(I,NTENTH,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(2) ! UPBAR DISTRIBUTION = D DISTRIBUTION
      DNV=G(2)
      SEA=G(4) ! THIS SEA IS (UBAR+DBAR)/2
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
      X=XSAVE
      RETURN
      END

c From W.J.Stirling@durham.ac.uk Fri Jul  9 17:59:15 1993
c To: WMELNITC@PHYSICS.ADELAIDE.EDU.AU
c Subject: MRSS pions fortran subroutine
