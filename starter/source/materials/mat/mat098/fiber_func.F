Chd|====================================================================
Chd|  FIBER_SFUNC                     source/materials/mat/mat098/lossfun_98.F
Chd|-- called by -----------
Chd|-- calls ---------------
Chd|        MESSAGE_MOD                   share/message_module/message_mod.F
Chd|        TABLE_MOD                     share/modules1/table_mod.F    
Chd|====================================================================
      SUBROUTINE FIBER_FUNC(NPTB  ,STR1   ,STR2   ,FLEX1   ,FLEX2   ,
     .                      XBIAX ,YBIAX1 ,YBIAX2 ,XFIB1   ,YFIB1   ,
     .                      XFIB2 ,YFIB2  )              
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "param_c.inc"
#include      "com04_c.inc"
#include      "scr17_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER ,INTENT(IN) :: NPTB
      my_real ,INTENT(IN) :: STR1,STR2,FLEX1,FLEX2
      my_real ,DIMENSION(NPTB) ,INTENT(IN)  :: XBIAX,YBIAX1,YBIAX2
      my_real ,DIMENSION(NPTB) ,INTENT(OUT) :: XFIB1,YFIB1,XFIB2,YFIB2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I,ITER,NITER
      my_real :: FLEX,DC0,HC0,DT0,HC,HT0,YC,LC,DC,HT,YT,LT,DT,Y,FC,FT,HDC,HDT,
     .           FPC,FPT,PHI,DPHI,STISSUC,STISSUT,UDC,UDT
c=======================================================================
      NITER = 5
      DC0 = ONE + STR1
      DT0 = ONE + STR2 
      HC0 = SQRT(DC0*DC0 - ONE)
      HT0 = SQRT(DT0*DT0 - ONE)
      FLEX  = FLEX1 + FLEX2
      Y   = ZERO
c
      DO I = 1, NPTB   
        LC = ONE + XBIAX(I)  
        LT = ONE + XBIAX(I) 
        STISSUC = YBIAX1(I)
        STISSUT = YBIAX2(I)

c         Y = ZERO
         DO ITER = 1, NITER
            HC = HC0 + Y 
            HT = HT0 - Y 
            DC = SQRT(LC *LC  + HC*HC) 
            DT = SQRT(LT *LT  + HT*HT)
            UDC= ONE / DC
            UDT= ONE / DT
            HDC= HC * UDC
            HDT= HT * UDT
            FC  = STISSUC * LC/DC
            FT  = STISSUT * LT/DT
            FPC = -FC * HDC/DC
            !FPC = STISSUC * HDC
            FPT = -FT * HDT/DT
            !FPT = STISSUT * HDT
            PHI  = FLEX*Y + FC * HDC - FT * HDT
            DPHI = FLEX + FPC*HDC + FC*UDC*(ONE - HDC*HDC)
     .                 + FPT*HDT + FT*UDT*(ONE - HDT*HDT)
            Y = Y - PHI / DPHI
C
            IF (Y > 0) THEN
              Y = MIN(Y, HT0)
            ELSE
              Y = MAX(Y,-HC0)
            ENDIF
          ENDDO !iter

          YC = Y
          YT =-Y
          HC = HC0 + YC   
          HT = HT0 + YT   
          DC = SQRT(LC *LC  + HC*HC) 
          DT = SQRT(LT *LT  + HT*HT) 
          XFIB1(I) = DC - DC0
          XFIB2(I) = DT - DT0
          YFIB1(I) = STISSUC * LC/DC
          YFIB2(I) = STISSUT * LT/DT


          PRINT*, 'C,XFIB, YFIB ', I, XFIB1(I), YFIB1(I)
          PRINT*, 'T,XFIB, YFIB ', I, XFIB2(I), YFIB2(I)
      ENDDO
c
c-----------
      RETURN
      END 
