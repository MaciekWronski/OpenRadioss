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
Chd|  SENSOR_OUTPUT               source/tools/sensor/SENSOR_OUTPUT.F
Chd|-- called by -----------
Chd|        RESOL                         source/engine/resol.F         
Chd|-- calls ---------------
Chd|        SENSOR_MOD                    share/modules/sensor_mod.F    
Chd|====================================================================
      SUBROUTINE SENSOR_OUTPUT(SENSORS,ISPMD,
     .           MX_OUTP2 ,MX_STAT2 ,MSTATT,MOUTPT,
     .           IWARN,MSTOP,MREST,NSANIM,MDESS,
     .           MSTAT,MOUTP)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------  
      USE SENSOR_MOD
      USE MESSAGE_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER ,INTENT(IN)    :: ISPMD
      INTEGER ,INTENT(IN)    :: MX_OUTP2
      INTEGER ,INTENT(IN)    :: MX_STAT2
      INTEGER ,INTENT(IN)    :: NSANIM
      INTEGER ,INTENT(INOUT) :: MSTATT
      INTEGER ,INTENT(INOUT) :: MOUTPT
      INTEGER ,INTENT(INOUT) :: IWARN
      INTEGER ,INTENT(INOUT) :: MSTOP
      INTEGER ,INTENT(INOUT) :: MREST
      INTEGER ,INTENT(INOUT) :: MDESS
      INTEGER ,INTENT(INOUT) :: MSTAT(MX_STAT2)
      INTEGER ,INTENT(INOUT) :: MOUTP(MX_OUTP2)
      TYPE (SENSORS_) ,INTENT(IN) :: SENSORS
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I,ISENS
      my_real :: TS,TT
C-----------------------------------------------
C   D e s c r i p t i o n
C-----------------------------------------------
C     this routine sets flags for radioss output files or run stop 
C     depending on sensor status
c=======================================================================
C     /STATE/LSENSOR
C---------------------------------------------
      IF (SENSORS%NSTAT  > 0) THEN
        MSTATT = 0
        DO I=1,SENSORS%NSTAT
           ISENS = SENSORS%STAT(I)
           TS = SENSORS%SENSOR_TAB(ISENS)%TSTART    
           IF (TT >= TS) THEN
             MSTAT(I) = MSTAT(I)+1    
           ENDIF
           IF (MSTAT(I)==1) MSTATT=1
        ENDDO
      ENDIF
C---------------------------------------------
C     /OUTP/LSENSOR
C---------------------------------------------
      IF (SENSORS%NOUTP > 0) THEN
        MOUTPT = 0
        DO I=1,SENSORS%NOUTP
           ISENS = SENSORS%OUTP(I)
           TS = SENSORS%SENSOR_TAB(ISENS)%TSTART    
           IF (TT >= TS) THEN
              MOUTP(I) = MOUTP(I)+1    
           ENDIF
           IF(MOUTP(I)==1)  MOUTPT=1
        ENDDO
      ENDIF
C---------------------------------------------
C     /STOP/LSENSOR
C---------------------------------------------
      IF (SENSORS%NSTOP > 0) THEN
        DO I=1,SENSORS%NSTOP
          ISENS = SENSORS%STOP(I)
          IF (ISPMD == 0 .and. SENSORS%SENSOR_TAB(ISENS)%STATUS == 1) THEN
            CALL ANCMSG(MSGID=234,ANMODE=ANINFO,
     .                  I1 = SENSORS%SENSOR_TAB(ISENS)%SENS_ID)
            IWARN = IWARN+1
            MSTOP = 1
            MREST = 1
            IF (NSANIM/=0) MDESS = 1
          ENDIF
        ENDDO
      ENDIF
C-----------
      RETURN
      END
