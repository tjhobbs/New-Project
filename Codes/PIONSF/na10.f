C *****************************************************************************
        SUBROUTINE NA10 (MODE,x,Q2in,xuv,xdv,xse,xVpi,xSpi)
C
C  Subroutine giving x-dependence of NA10 parametrizations of proton valence
C  and sea distribution functions xuv, xdv and xse (=xu_sea=xd_sea=2xs_sea),
C  evolved from Q2=5GeV^2,
C  and pion valence and sea distribution functions xVpi, xSpi, evolved from
C  Q2=25GeV^2.
C
C  NA10 Collaboration, B.Betev, et al., Z.Phys. C28 (1985) 15 (table 1).
C *****************************************************************************
	IMPLICIT NONE
        INTEGER MODE
        REAL*4  x,xuv,xdv,xse,xVpi,xSpi
        REAL*4  a,b,C,beta,a0,b0,a1,b1,gam,gpi,Nv,Ns,EULBETA
        REAL*4  Q2in,Q2,Q02,L2,s

	IF (x.LE.1.E-6 .OR. x.GE.0.99) RETURN

C...MODE 1: Evaluate pion distributions
C...MODE 2: Evaluate proton distributions
        IF (MODE.EQ.1) GO TO 100
        IF (MODE.EQ.2) GO TO 200


C...Pion distribution functions................................................
 100    xVpi = 0.0
        xSpi = 0.0
	IF (x.LE.0.0 .OR. x.GE.0.99) RETURN 
        Q02 = 25.0
        Q2 = Q2in
        L2 = 0.09
        IF (Q2.LE.Q02) Q2 = 25.001
        IF (Q2.LE.L2) RETURN
        s = ALOG( ALOG(Q2/L2) / ALOG(Q02/L2) )

C...Valence distribution
        a0 = 0.41
        b0 = 1.11
        a1 = 0.0579 - 0.2370*a0 - 0.0328*b0
        b1 = 0.5150 - 0.0114*a0 + 0.0544*b0
        a = a0 + a1*s
        b = b0 + b1*s
C...Valence normalization: S dx Vpi(x) = 1
        Nv = 1.0 / EULBETA(a,b+1.0)
        xVpi = 0.0
        IF (x.LT.1.0) xVpi = Nv * x**a * (1.0-x)**b
        IF (xVpi.LT.0.0) xVpi = 0.0

C...Sea distribution
        gam = 8.4
        gpi = 0.47
C...Sea normalized s.t. S dx (2 Vpi(x) + 6 Spi(x)) = 1 - <gpi>,
C...where <gpi> = gluon fractional momentum
        Ns = 1.0 / (6.0 * EULBETA(1.0,gam+1.0))
     &     * (1.0 - gpi - 2.0*EULBETA(a+1.0,b+1.0) / EULBETA(a,b+1.0) )
        xSpi = 0.0
        IF (x.LT.1.0) xSpi = Ns * (1.0-x)**gam
        IF (xSpi.LT.0.0) xSpi = 0.0
        RETURN



C...Proton distribution functions..............................................
 200    Q02 = 5.0
        Q2 = Q2in
        IF (Q2.LE.Q02) Q2 = 5.001
        L2 = 0.09
        s = ALOG( ALOG(Q2/L2) / ALOG(Q02/L2) )

C...u and d valence distributions
        a = 0.6190 - 0.1678*s
        b = 2.8670 + 0.6687*s
        xuv = 2.0 / EULBETA(a,b+1.0) * x**a * (1.0-x)**b
        xdv = 1.0 / EULBETA(a,b+2.0) * x**a * (1.0-x)**(b+1.0)
C...Sea distribution
        C = (0.50758 + 0.23006*s + 0.067345*s**2) / 2.8
        beta = 7.417 - 1.138*s + 13.22*s**2 - 4.966*s**3 - 1.86*s**4
        xse = C * (1.0-x)**beta
        RETURN

        END


C *****************************************************************************
        FUNCTION EULBETA (X,Y)
C  Euler's Beta function.
C *****************************************************************************
c        IMPLICIT UNDEFINED (A-Z)
	IMPLICIT NONE
        REAL*4  X,Y,GAMA,EULBETA

        EULBETA = GAMA(X) * GAMA(Y) / GAMA(X+Y)
        RETURN
        END


C *****************************************************************************
        FUNCTION GAMA (X)
C  Gamma function.
C *****************************************************************************
c        IMPLICIT UNDEFINED (A-Z)
	IMPLICIT NONE
        REAL*4  C(13),Z,X,F,GAMA
        DATA    C
     & /0.00053 96989 58808, 0.00261 93072 82746, 0.02044 96308 23590,
     &  0.07309 48364 14370, 0.27964 36915 78538, 0.55338 76923 85769,
     &  0.99999 99999 99998,-0.00083 27247 08684, 0.00469 86580 79622,
     &  0.02252 38347 47260,-0.17044 79328 74746,-0.05681 03350 86194,
     &  1.13060 33572 86556/

        Z=X
        IF (X .GT. 0.0) GO TO 1
        IF (X .EQ. AINT(X)) GO TO 5
        Z=1.0-Z
 1      F=1.0/Z
        IF (Z .LE. 1.0) GO TO 4
        F=1.0
 2      IF (Z .LT. 2.0) GO TO 3
        Z=Z-1.0
        F=F*Z
        GO TO 2
 3      Z=Z-1.0
 4      GAMA=
     &   F*((((((C(1)*Z+C(2))*Z+C(3))*Z+C(4))*Z+C(5))*Z+C(6))*Z+C(7))/
     &   ((((((C(8)*Z+C(9))*Z+C(10))*Z+C(11))*Z+C(12))*Z+C(13))*Z+1.0)
        IF (X .GT. 0.0) RETURN
        GAMA=3.141592653589793/(SIN(3.141592653589793*X)*GAMA)
        RETURN
 5      GAMA=0.0
        RETURN
c 10    FORMAT(1X,45HGAMA ... ARGUMENT IS NON-POSITIVE INTEGER = ,E20.5)
        END
