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
      !||    thermc        ../engine/source/materials/mat_share/thermc.F
      !||--- called by ------------------------------------------------------
      !||    cforc3        ../engine/source/elements/shell/coque/cforc3.F
      !||    cforc3_crk    ../engine/source/elements/xfem/cforc3_crk.F
      !||    czforc3       ../engine/source/elements/shell/coquez/czforc3.F
      !||    czforc3_crk   ../engine/source/elements/xfem/czforc3_crk.F
      !||====================================================================
      SUBROUTINE THERMC(NEL    ,THERM ,THK   , 
     &                  PX1    ,PX2   ,PY1   ,PY2   ,AREA      ,timestep  ,
     &                  TEMPN  ,TEMPEL,DHEAT ,FPHI  ,THEACCFACT,
     &                  IX     ,NIX   ,NUMEL ,NUMNOD,MVSIZ )
C-----------------------------------------------
C    calculates nodal thermic force from heat energy increment
C-----------------------------------------------
C    M o d u l e s
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
      integer ,intent(in)  :: NUMNOD
      integer ,intent(in)  :: NIX
      integer ,intent(in)  :: NUMEL
      integer ,intent(in)  :: MVSIZ
      integer ,intent(in)  :: IX(NIX,NUMEL)
      my_real ,intent(in)  :: THEACCFACT
      my_real ,intent(in)  :: timestep
      my_real ,intent(out) :: FPHI(MVSIZ,4)
      my_real ,DIMENSION(NUMNOD) ,intent(in) :: TEMPN
      my_real ,DIMENSION(NEL)    ,intent(in) :: TEMPEL
      my_real ,DIMENSION(NEL)    ,intent(in) :: AREA
      my_real ,DIMENSION(NEL)    ,intent(in) :: PX1
      my_real ,DIMENSION(NEL)    ,intent(in) :: PX2
      my_real ,DIMENSION(NEL)    ,intent(in) :: PY1
      my_real ,DIMENSION(NEL)    ,intent(in) :: PY2
      my_real ,DIMENSION(NEL)    ,intent(in) :: DHEAT
      my_real ,DIMENSION(NEL)    ,intent(in) :: THK
      type (therm_param_)        ,intent(in) :: therm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I,N1,N2,N3,N4
      my_real :: AS,BS,A,KC,PHIX,PHIY
C-----------------------------------------------
      AS = THERM%AS
      BS = THERM%BS
!
      DO I=1,NEL
        KC = (AS + BS*TEMPEL(I))*timestep*THEACCFACT / MAX(EM20,AREA(I))
        N1 = IX(2,I)
        N2 = IX(3,I)
        N3 = IX(4,I)
        N4 = IX(5,I)
        !     
        ! thermal flux : K gradT = B * T
        !     
        PHIX = TEMPN(N1)*PX1(I)+TEMPN(N2)*PX2(I)-TEMPN(N3)*PX1(I)-TEMPN(N4)*PX2(I)
        PHIY = TEMPN(N1)*PY1(I)+TEMPN(N2)*PY2(I)-TEMPN(N3)*PY1(I)-TEMPN(N4)*PY2(I)
C     
        PHIX = KC*PHIX*THK(I)
        PHIY = KC*PHIY*THK(I)
        !     
        ! nodal thermal force
        !     
        A = FOURTH * DHEAT(I)
        FPHI(I,1) = A - (PHIX*PX1(I) + PHIY*PY1(I))
        FPHI(I,2) = A - (PHIX*PX2(I) + PHIY*PY2(I))
        FPHI(I,3) = A + (PHIX*PX1(I) + PHIY*PY1(I))
        FPHI(I,4) = A + (PHIX*PX2(I) + PHIY*PY2(I))
      ENDDO
!-----------
      RETURN
      END
