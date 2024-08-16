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
      !||    s20therm   ../engine/source/elements/solid/solide20/s20therm.F
      !||--- called by ------------------------------------------------------
      !||    s16forc3   ../engine/source/elements/thickshell/solide16/s16forc3.F
      !||    s20forc3   ../engine/source/elements/solid/solide20/s20forc3.F
      !||====================================================================
      SUBROUTINE S20THERM(
     1   THERM   ,NEL     ,NPE,     NC,
     2   VOL     ,PX,      PY,      PZ,
     3   NI      ,DT1,     TEMPNC,  TEMPEL,
     4   HEAT    ,FPHI,    OFFG,    OFF,
     5   THEACCFACT)
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
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "mvsiz_p.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: NEL
      INTEGER NPE,NC(MVSIZ,NPE)
      my_real, INTENT(IN) :: THEACCFACT
      my_real :: VOL(*), PX(MVSIZ,NPE), PY(MVSIZ,NPE),PZ(MVSIZ,NPE),
     .    NI(NPE),TEMPNC(*), FPHI(MVSIZ,NPE), HEAT(*),
     .    DT1,TEMPEL(*),OFF(*),OFFG(*)
      type (therm_param_) ,intent(in) :: therm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I,N
      my_real :: AS, BS, KC, PHIX(MVSIZ), PHIY(MVSIZ), PHIZ(MVSIZ)
!==============================================================================
      AS = THERM%AS    !< thermal conductivity coefficient A for solid phase
      BS = THERM%BS    !< thermal conductivity coefficient B for solid phase
!
      PHIX(1:MVSIZ) = ZERO
      PHIY(1:MVSIZ) = ZERO
      PHIZ(1:MVSIZ) = ZERO
C 
C - flux
C
      DO N=1,NPE
        DO I= 1,NEL
         IF(OFF(I)==ZERO.OR.OFFG(I)<=ZERO) CYCLE
          PHIX(I) = PHIX(I) +  TEMPNC(NC(I,N))*PX(I,N)
          PHIY(I) = PHIY(I) +  TEMPNC(NC(I,N))*PY(I,N)
          PHIZ(I) = PHIZ(I) +  TEMPNC(NC(I,N))*PZ(I,N)
        ENDDO
      ENDDO
      DO I= 1,NEL
        IF(OFF(I)==ZERO.OR.OFFG(I)<=ZERO) CYCLE
        KC = (AS + BS*TEMPEL(I))*DT1*VOL(I)*THEACCFACT          
        PHIX(I) = KC*PHIX(I)
        PHIY(I) = KC*PHIY(I)
        PHIZ(I) = KC*PHIZ(I)
      ENDDO
C
C force thermique nodale
C
      DO N=1,NPE
       DO I= 1,NEL
         IF (OFF(I)==ZERO.OR.OFFG(I)<=ZERO) CYCLE        
         FPHI(I,N) = FPHI(I,N) + HEAT(I)*NI(N)
     .             - (PHIX(I)*PX(I,N) + PHIY(I)*PY(I,N) + PHIZ(I)*PZ(I,N))
       ENDDO
      ENDDO
!----------- 
      RETURN
      END
