Chd|====================================================================
Chd|  JAC98_UNI                   source/materials/mat/mat098/lossfun_98.F
Chd|-- called by -----------
Chd|  LAW98_IDENT                     source/materials/mat/mat098/law98_upd.F
Chd|-- calls ---------------
Chd|        FINTER58                      source/materials/mat/mat098/lossfun_98.F
CChd|====================================================================
      SUBROUTINE JAC98_UNI(NVAR   ,NPTU   ,NPTB   ,XUNIAX ,YUNIAX ,
     .                     XFIB   ,YFIB   ,STRETCH,FLEX1  ,FLEX2  ,
     .                     WEIGHT ,FSQR   ,JAC    ,IVAR   )
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
c----------------------------------------------- 
      INTEGER ,INTENT(IN)  :: NPTB,NPTU,NVAR,IVAR
      my_real ,INTENT(IN)  :: STRETCH,FLEX1,FLEX2
      my_real ,DIMENSION(NPTU) ,INTENT(IN)  :: XUNIAX,YUNIAX
      my_real ,DIMENSION(NPTB) ,INTENT(IN)  :: XFIB,YFIB
c
      my_real ,INTENT(OUT) :: WEIGHT
      my_real ,DIMENSION(NPTU) ,INTENT(OUT) :: FSQR
      my_real ,DIMENSION(NVAR,NPTU) ,INTENT(OUT) :: JAC
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c---------------------------- 
      INTEGER :: IPT,ITER,NITER
      my_real :: DC0,DC,HC0,HC,YC,DCC,UDC,HDC,LC,
     .           FLEX,FLEXP,FLEX1P,FSQRP,FSQR2,DFLEX,DSTR,FC,FPC,
     .           AA,PHI,DPHI,SIG
      my_real :: FINTER58
      EXTERNAL FINTER58
c----------------------------------------------- 
c     sum of squared differences between experimental input curve 
c     and calculated response in uniaxial tension test
c     SUM ( (Fexp(eps_i) - Fsim(eps_i)) / Fexp(eps_i))**2
c=======================================================================
      NITER = 5
      DC0   = ONE + STRETCH
      HC0   = SQRT(DC0*DC0 - ONE)
      FLEX  = FLEX1 + FLEX2
      YC    = ZERO
c      
      DO IPT = 1,NPTU
        LC = ONE + XUNIAX(IPT)
c
        DO ITER = 1, NITER
          HC   = HC0 + YC 
          DC   = SQRT(LC*LC + HC*HC) 
          DCC  = DC - DC0
          UDC  = ONE / DC
          HDC  = HC * UDC
          FC   = FINTER58(DCC,XFIB,YFIB,FPC,NPTB)
          FPC  = FPC*HDC
          PHI  = FLEX*YC + FC*HDC
          DPHI = FLEX + FPC*HDC + FC*UDC*(ONE - HDC*HDC)
          YC   = YC - PHI / MAX(DPHI,EM20)
        ENDDO ! iter
        YC  = MAX(YC,-HC0)
        HC  = HC0 + YC  
        DC  = SQRT(LC*LC + HC*HC) 
        DCC = DC - DC0    
        FC  = FINTER58(DCC,XFIB,YFIB,FPC,NPTB)
c
        SIG = FC * LC / DC   ! ENGINEERING STRESS (initial state => EC2=0)
        AA  = (YUNIAX(IPT) - SIG) / MAX(EM20,YUNIAX(IPT))
c        FSQR(IPT) = WEIGHT * AA**2  ! sum of squared differences : target func to minimize
        FSQR(IPT) = AA
      ENDDO
c------------------------------------------------------
c     Calculate target function derivatives  (jacobian)
c------------------------------------------------------
      DSTR  = MAX(EM03 * STRETCH  , EM10)
      DFLEX = MAX(EM03 * FLEX1 , EM10)

c     over stretch
c
      FLEX  = FLEX1 + FLEX2
      DC0   = ONE + STRETCH + DSTR
      HC0   = SQRT(DC0*DC0 - ONE)
      YC    = ZERO
      JAC(:,:) = ZERO
c
      DO IPT = 1,NPTU
        LC = ONE + XUNIAX(IPT)
c
        DO ITER = 1, NITER
          HC   = HC0 + YC 
          DC   = SQRT(LC*LC + HC*HC) 
          DCC  = DC - DC0
          UDC  = ONE / DC
          HDC  = HC * UDC
          FC   = FINTER58(DCC,XFIB,YFIB,FPC,NPTB)
          FPC  = FPC*HDC
          PHI  = FLEX*YC + FC*HDC
          DPHI = FLEX + FPC*HDC + FC*UDC*(ONE - HDC*HDC)
          YC   = YC - PHI / MAX(DPHI,EM20)
        ENDDO !iter
        YC  = MAX(YC,-HC0)
        HC  = HC0 + YC  
        DC  = SQRT(LC*LC + HC*HC) 
        DCC = DC - DC0    
        FC  = FINTER58(DCC,XFIB,YFIB,FPC,NPTB)
c
        SIG = FC * LC / DC   ! ENGINEERING STRESS (initial state => EC2=0)
        AA  = (YUNIAX(IPT) - SIG) / MAX(EM20,YUNIAX(IPT))
c        FSQRP = WEIGHT * AA**2  ! sum of squared differences : target func to minimize
        FSQRP = AA
        JAC(IVAR,IPT) = (FSQR(IPT) - FSQRP) / DSTR
      ENDDO 
c
c     over flex
      FLEX  = FLEX1 + FLEX2 + DFLEX
      DC0   = ONE + STRETCH
      HC0   = SQRT(DC0*DC0 - ONE)
      FLEX  = FLEX1P + FLEX2
      YC    = ZERO
c
      DO IPT = 1,NPTU
        LC = ONE + XUNIAX(IPT)
c
        DO ITER = 1, NITER
          HC   = HC0 + YC 
          DC   = SQRT(LC*LC + HC*HC) 
          DCC  = DC - DC0
          UDC  = ONE / DC
          HDC  = HC * UDC
          FC   = FINTER58(DCC,XFIB,YFIB,FPC,NPTB)
          FPC  = FPC*HDC
          PHI  = FLEX*YC + FC*HDC
          DPHI = FLEX + FPC*HDC + FC*UDC*(ONE - HDC*HDC)
          YC   = YC - PHI / MAX(DPHI,EM20)
        ENDDO !iter
        YC  = MAX(YC,-HC0)
        HC  = HC0 + YC  
        DC  = SQRT(LC*LC + HC*HC) 
        DCC = DC - DC0    
        FC  = FINTER58(DCC,XFIB,YFIB,FPC,NPTB)
c
        SIG = FC * LC / DC   ! ENGINEERING STRESS (initial state => EC2=0)
        AA  = (YUNIAX(IPT) - SIG) / MAX(EM20,YUNIAX(IPT))
c        FSQRP = WEIGHT * AA**2  ! sum of squared differences : target func to minimize
        FSQRP = AA
        JAC(3,IPT) = (FSQR(IPT) - FSQRP) / DFLEX
        JAC(4,IPT) = JAC(3,IPT)
      ENDDO
c-----------
      RETURN
      END 
