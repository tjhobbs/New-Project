C **************************************************************************
	SUBROUTINE E615 (x,xVpi,xSpi)
C
C  Subroutine giving x-dependence of E615 parametrizations of pion valence
C  and sea distribution functions xVpi, xSpi, at Q2 = 16.4 to 73.1 GeV^2.
C
C  E615 Collaboration, J.S.Conway, et al., P.R.D39 (1989) 92 (4 of table 5).
C *************************************************************************** 
        IMPLICIT UNDEFINED (A-Z)
        REAL*8  x,xVpi,xSpi
        REAL*8  a,b,gamma,delta,Av,As, Q2
        LOGICAL rem
        DATA    rem /.TRUE./
 
C...Reminder about Q2
        IF (rem) THEN
          WRITE (6,*) '*** E615 PION DISTRIBUTIONS AT Q2=16-73GeV^2 ***'
          rem = .FALSE.
        ENDIF

	Q2 = 16.D0 ! GeV^2 average
 
C...Valence distribution function, with normalization: S dx Vpi(x) = 1
C... =>  Av = 1 / Beta[a,b+1] = 0.929026 @ Q2 ~ 25 Gev^2
        a = 0.56D0
        b = 1.21D0
        gamma = 0.83D0	!0.63D0
        Av = 0.929026D0
        xVpi = 0.D0
        IF (x.LE.1.D0) xVpi = Av * x**a * (1.D0-x)**b
     &                      + gamma * 2.D0 * x**2 / (9.D0 * Q2)
 
C...Sea distribution, normalized s.t. S dx (2 Vpi(x) + 6 Spi(x)) = 1 - <gpi>,
C                                                    where <gpi> = 0.47
C      => Ns = 1 / (6 Beta[1,gamma+1])
C		 (1 - <gpi> - 2 Beta[a+1,b+1]/Beta[a,b+1])
C...         = 0.19688 @ Q2 ~ 25 GeV^2
        As = 0.19688D0
        delta = 8.4D0
        xSpi = 0.D0
        IF (x.LT.1.D0) xSpi = As * (1.D0-x)**delta
 
        END
