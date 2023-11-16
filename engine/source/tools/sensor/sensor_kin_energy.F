Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2023 Altair Engineering Inc.
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
Chd|====================================================================
Chd|  SENSOR_KIN_ENERGY               source/tools/sensor/SENSOR_KIN_ENERGY.F
Chd|-- called by -----------
Chd|        RESOL                         source/engine/resol.F         
Chd|-- calls ---------------
Chd|        SENSOR_MOD                    share/modules/sensor_mod.F    
Chd|====================================================================
      SUBROUTINE SENSOR_KIN_ENERGY(
     .           NUMNOD   ,NSENSOR  ,SENSOR_TAB ,ENCIN   ,MS       ,
     .           V        ,A        ,WEIGHT     ,DT1     ,IPARIT   ,
     .           NSPMD    )
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------  
      USE SENSOR_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "task_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER ,INTENT(IN) :: IPARIT
      INTEGER ,INTENT(IN) :: NSPMD
      INTEGER ,INTENT(IN) :: NUMNOD
      INTEGER ,INTENT(IN) :: NSENSOR
      my_real ,INTENT(IN) :: DT1
      my_real ,INTENT(IN) :: ENCIN
      INTEGER ,DIMENSION(NUMNOD)   ,INTENT(IN) :: WEIGHT
      my_real ,DIMENSION(NUMNOD)   ,INTENT(IN) :: MS
      my_real ,DIMENSION(3,NUMNOD) ,INTENT(IN) :: V,A
      TYPE (SENSOR_STR_) ,DIMENSION(NSENSOR) ,INTENT(INOUT) :: SENSOR_TAB
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I,ISENS,ISELECT,IKINE
      my_real :: VX,VY,VZ,MAS,EKIN,DT05
      my_real :: ENTMP(1)
C-----------------------------------------------
C   D e s c r i p t i o n
C-----------------------------------------------
C     this routine calculates the total kinematic system energy for 
C     /SENSOR/ENERGY
c=======================================================================
      IKINE = 0
      DO ISENS=1,NSENSOR
        IF (SENSOR_TAB(ISENS)%TYPE == 14) THEN   ! /sensor/energy
          ISELECT = SENSOR_TAB(ISENS)%IPARAM(4)  ! Iselect == 2 => total kinetic system energy
          IF (ISELECT == 2) THEN
            IKINE = 1
            EXIT
          END IF
        ENDIF
      ENDDO
c
      IF (IKINE == 1) THEN
        EKIN = ZERO
        DT05 = HALF*DT1
        IF (IPARIT == 0) THEN
          DO I = 1, NUMNOD
            IF (WEIGHT(I) == 1) THEN
              MAS  = MS(I) * WEIGHT(I)
              VX   = V(1,I) + DT05*A(1,I)
              VY   = V(2,I) + DT05*A(2,I)
              VZ   = V(3,I) + DT05*A(3,I)
              EKIN = EKIN  + (VX*VX + VY*VY + VZ*VZ)*HALF*MAS
            END IF 
          ENDDO
        ELSE
          DO I = 1, NUMNOD
            IF (WEIGHT(I) == 1) THEN
              MAS  = MS(I) * WEIGHT(I)
              VX   = V(1,I) + DT05*A(1,I)
              VY   = V(2,I) + DT05*A(2,I)
              VZ   = V(3,I) + DT05*A(3,I)
              EKIN = EKIN  + (VX*VX + VY*VY + VZ*VZ)*HALF*MAS
            END IF 
          ENDDO
        END IF
c
        IF (NSPMD > 1) THEN
          ENTMP(1) = ENCIN
          CALL SPMD_GLOB_DSUM9(ENTMP,1)
          IF (ISPMD == 0) THEN
            EKIN = ENTMP(1)
          END IF
        END IF
c
        IF (ISPMD == 0) THEN
          DO ISENS=1,NSENSOR
            IF (SENSOR_TAB(ISENS)%TYPE == 14) THEN 
              ISELECT = SENSOR_TAB(ISENS)%IPARAM(4)
              IF (ISELECT == 2) THEN
                SENSOR_TAB(ISENS)%VAR(9) = EKIN + ENCIN
              END IF
            ENDIF
          ENDDO
        END IF
      ENDIF
C-----------
      RETURN
      END
