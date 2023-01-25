Chd|====================================================================
Chd|  UNIFY_2TAB_X                       source/tools/curve/func_slope.F
Chd|-- called by -----------
Chd|-- calls ---------------
Chd|        MESSAGE_MOD                   share/message_module/message_mod.F
Chd|        TABLE_MOD                     share/modules1/table_mod.F    
Chd|====================================================================
      SUBROUTINE UNIFY_2TAB_X(FUNC1,FUNC2,TABLE,LEN1,LEN2,NPOINT,XUNI)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
      USE TABLE_MOD
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
      INTEGER :: FUNC1,FUNC2,NPOINT,LEN1,LEN2
      my_real :: XUNI(NPOINT) 
      TYPE (TTABLE), DIMENSION(NTABLE) ,TARGET ::  TABLE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I,J,K
      my_real :: X1,X2
c=======================================================================
      I = 1
      J = 1
      K = 0
      DO While (I <= LEN1 .and. J <= LEN2)
        X1 = TABLE(FUNC1)%X(1)%VALUES(I)
        X2 = TABLE(FUNC2)%X(1)%VALUES(J)
        K = K + 1
        
        IF (X1 == X2) THEN
          I = I + 1
          J = J + 1
          XUNI(K) = X2 
        ELSE IF (X1 < X2) THEN
          I  = I + 1
          XUNI(K) = X1 
        ELSE IF (X1 > X2) THEN  
          J = J + 1
          XUNI(K) = X2 
        ENDIF
      ENDDO
c
      NPOINT = K   ! updated common abscissa length
c-----------
      RETURN
      END
