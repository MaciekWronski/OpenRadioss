Chd|====================================================================
Chd|  UNIFY_X                       source/tools/curve/func_slope.F
Chd|-- called by -----------
Chd|-- calls ---------------
Chd|        MESSAGE_MOD                   share/message_module/message_mod.F
Chd|        TABLE_MOD                     share/modules1/table_mod.F    
Chd|====================================================================
      SUBROUTINE UNIFY_X(IDN1,IDN2,NPC,PLD,NPOINT,LEN1,LEN2,XUNI,NPTNEW)
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
      INTEGER     IDN1,IDN2,IDN,NPOINT,LEN1,LEN2,
     .        NPTNEW,NPC(*) 
      my_real 
     .        PLD(*),XUNI(NPOINT) 
C-----------------------------------------------
      INTENT(IN)    :: NPC,PLD
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I,J,K,ID,NP1,NP2,K1,PN1,PN2
      my_real 
     .      EC , ET  
c=======================================================================
       !IDN1 = IFUNC(1) ! uni C
       !IDN2 = IFUNC(2) ! uni T
       I = 0
       J = 0
       EC = PLD(NPC(IDN1) )
       ET = PLD(NPC(IDN2) )
       DO K = 1,NPOINT
        IF(I == 2*LEN1 .AND. J == 2*LEN2 )THEN 
         EXIT
        ELSE
         IF     ((EC < ET.AND.I<2*LEN1) .OR. J >= 2*LEN2)THEN
           XUNI(K) = EC 
           I = I + 2 
           EC = PLD(NPC(IDN1)+ I )
         ELSEIF ((EC > ET.AND.J<2*LEN2)   .OR. I >= 2*LEN1)THEN  
           XUNI(K) = ET 
           J = J + 2
           ET = PLD(NPC(IDN2)+ J )
         ELSEIF (EC == ET)THEN 
           XUNI(K) = ET 
           I = I + 2
           J = J + 2
           EC = PLD(NPC(IDN1)+ I )
           ET = PLD(NPC(IDN2)+ J )
         ENDIF
        ENDIF
       ENDDO
       NPTNEW = K
c-----------
      RETURN
      END
