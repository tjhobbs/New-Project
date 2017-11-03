C ******************************************************************************
        SUBROUTINE GRV (x,Q2,xVpi,xSpi)
C
C  Subroutine giving x-dependence of parametrization of LO pion valence
C  and sea distribution functions xVpi, xSpi, valid for 0.25 < Q2 < 10^8 GeV^2,
C  and 10^-5 < x < 1.
C
C  Gluck, Reya, Vogt: Z.Phys. C53 (1992) 651 (appendix 1).
C ******************************************************************************
!        IMPLICIT UNDEFINED (A-Z)
        REAL*8  x,xVpi,a,AA,D,Nv,
     &          xSpi,alpha,as,AAs,Bs,Ds,E,Epr,beta
        REAL*8  Q2,Q02,L,s

        xVpi = 0.D0
        xSpi = 0.D0

	IF (x.LT.1.D-5) x=1.01D-5
	IF (Q2.LT.0.25D0) Q2=0.2501D0

        Q02 = 0.25D0
        L = 0.232D0
        IF (Q2.LE.Q02) RETURN
        s = DLOG( DLOG(Q2/L**2) / DLOG(Q02/L**2) )

C...Valence distribution
        Nv = 0.519D0 + 0.180D0*s - 0.011D0*s**2
        a = 0.499D0 - 0.027D0*s
        AA = 0.381D0 - 0.419D0*s
        D = 0.367D0 + 0.563D0*s
        xVpi = 0.D0
        IF (x.LT.1.D0)
     &	  xVpi = Nv * x**a * (1.D0+AA*DSQRT(x)) * (1.D0-x)**D

C...Sea distribution (SU(3) symmetric)
        alpha = 0.55D0
        as = 2.538D0 - 0.763D0*s
        AAs = -0.748D0
        Bs = 0.313D0 + 0.935D0*s
        Ds = 3.359D0
        E = 4.433D0 + 1.301D0*s
        Epr = 9.30D0 - 0.887D0*s
        beta = 0.56D0
        xSpi = 0.D0
        IF (x.LT.1.D0)
     &	  xSpi = s**alpha / (DLOG(1.D0/x))**as
     &         * (1.D0 + AAs*DSQRT(x) + Bs*x) * (1.D0 - x)**Ds
     &         * DEXP(-E + DSQRT(Epr*s**beta*DLOG(1.D0/x)))

        RETURN
        END
