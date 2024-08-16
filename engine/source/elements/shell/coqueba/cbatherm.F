Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      !||====================================================================
      !||    cbatherm   ../engine/source/elements/shell/coqueba/cbatherm.F
      !||--- called by ------------------------------------------------------
      !||    cbaforc3   ../engine/source/elements/shell/coqueba/cbaforc3.F
      !||====================================================================
      SUBROUTINE CBATHERM(NEL    ,THERM  ,THK   ,
     &                    BM     ,AREA   ,DTIME ,TEMPN  ,TEMPEL ,DHEAT  ,
     &                    NPLAT  ,IPLAT  ,FPHI  ,THEACCFACT,
     &                    IX     ,NIX    ,NUMEL ,NUMNOD ,MVSIZ  )
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use therm_param_mod
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ,intent(in)  :: NEL
      integer ,intent(in)  :: NPLAT
      integer ,intent(in)  :: NUMNOD
      integer ,intent(in)  :: NIX
      integer ,intent(in)  :: NUMEL
      integer ,intent(in)  :: MVSIZ
      integer ,intent(in)  :: IPLAT(NEL)
      integer ,intent(in)  :: IX(NIX,NUMEL)
      my_real ,intent(in)  :: THEACCFACT
      my_real ,intent(in)  :: dtime
      my_real ,intent(in)  :: BM(MVSIZ,7)
      my_real ,intent(out) :: FPHI(MVSIZ,4)
      my_real ,DIMENSION(NUMNOD) ,intent(in) :: TEMPN
      my_real ,DIMENSION(NEL)    ,intent(in) :: TEMPEL
      my_real ,DIMENSION(NEL)    ,intent(in) :: AREA
      my_real ,DIMENSION(NEL)    ,intent(in) :: DHEAT
      my_real ,DIMENSION(NEL)    ,intent(in) :: THK
      type (therm_param_)        ,intent(in) :: therm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I,EP,N1,N2,N3,N4
      my_real :: CA,CB,A,KC,PHIX,PHIY
      my_real :: TEMP1,TEMP2,TEMP3,TEMP4,TEMP13,TEMP24,TEMPH
C------------------------------------------
      CA = THERM%AS
      CB = THERM%BS
!
      DO EP=1,NPLAT 
        I=IPLAT(EP)

        KC = (CA + CB*TEMPEL(I))*DTIME*THEACCFACT
C
        N1 = IX(2,I)
        N2 = IX(3,I)
        N3 = IX(4,I)
        N4 = IX(5,I)

        TEMP13 = TEMPN(N1)-TEMPN(N3)
        TEMP24 = TEMPN(N2)-TEMPN(N4)
        TEMPH  = TEMPN(N1)-TEMPN(N2)+TEMPN(N3)-TEMPN(N4)
      !     
      ! thermal flux : K gradT = B * T
      !     
        PHIX =  TEMP13*BM(I,1) + TEMP24*BM(I,2)+BM(I,3)*TEMPH        
        PHIY =  TEMP13*BM(I,5) + TEMP24*BM(I,6)+BM(I,7)*TEMPH
!     
        PHIX = KC*PHIX*THK(I)*AREA(I)
        PHIY = KC*PHIY*THK(I)*AREA(I)
        !     
        ! nodal thermal force
        !     
        A = FOURTH *FOURTH * DHEAT(I)

        FPHI(I,1) = A - (PHIX*(BM(I,1)+BM(I,3)) + PHIY*(BM(I,5)+BM(I,7)))
        FPHI(I,2) = A - (PHIX*(BM(I,2)-BM(I,3)) + PHIY*(BM(I,6)-BM(I,7)))
        FPHI(I,3) = A - (PHIX*(BM(I,3)-BM(I,1)) + PHIY*(BM(I,7)-BM(I,5)))
        FPHI(I,4) = A + (PHIX*(BM(I,2)+BM(I,3)) + PHIY*(BM(I,6)+BM(I,7)))

      END DO
!--------------------------------------------------------------------------
      DO EP=NPLAT+1,NEL
        I=IPLAT(EP)

        KC = (CA + CB*TEMPEL(I)) * DTIME
C
        N1 = IX(2,I)
        N2 = IX(3,I)
        N3 = IX(4,I)
        N4 = IX(5,I)

        TEMP1 = TEMPN(N1)
        TEMP2 = TEMPN(N2)
        TEMP3 = TEMPN(N3)
        TEMP4 = TEMPN(N4)
        !     
        ! thermal flux : K gradT = B * T
        !     
        PHIX =  (BM(I,1)+BM(I,4))*TEMP1 + (BM(I,10)+BM(I,13))*TEMP2
     .       +  (BM(I,19)+BM(I,22))*TEMP3 + (BM(I,28)+BM(I,31))*TEMP2  
        PHIY =  (BM(I,2)+BM(I,5))*TEMP1 + (BM(I,11)+BM(I,14))*TEMP2
     .       +  (BM(I,20)+BM(I,23))*TEMP3 + (BM(I,29)+BM(I,32))*TEMP2  

        PHIX = KC*PHIX*THK(I)*AREA(I)
        PHIY = KC*PHIY*THK(I)*AREA(I)    
        !     
        ! nodal thermal force
        !     
        A = FOURTH *FOURTH * DHEAT(I)

        FPHI(I,1) = A - (PHIX*(BM(I,1)+BM(I,4))   + PHIY*(BM(I,2)+BM(I,5)))
        FPHI(I,2) = A - (PHIX*(BM(I,10)+BM(I,13)) + PHIY*(BM(I,11)+BM(I,14)))
        FPHI(I,3) = A - (PHIX*(BM(I,19)+BM(I,22)) + PHIY*(BM(I,20)+BM(I,23)))
        FPHI(I,4) = A - (PHIX*(BM(I,28)+BM(I,31))  + PHIY*(BM(I,29)+BM(I,32)))
      END DO
!----------- 
      RETURN
      END
