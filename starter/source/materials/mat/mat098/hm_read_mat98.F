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
Chd|  HM_READ_MAT98                 source/materials/mat/mat088/hm_read_mat98.F
Chd|-- called by -----------
Chd|        HM_READ_MAT                   source/materials/mat/hm_read_mat.F
Chd|-- calls ---------------
Chd|        ANCMSG                        source/output/message/message.F
Chd|        HM_GET_FLOATV                 source/devtools/hm_reader/hm_get_floatv.F
Chd|        HM_GET_FLOATV_DIM             source/devtools/hm_reader/hm_get_floatv_dim.F
Chd|        HM_GET_INTV                   source/devtools/hm_reader/hm_get_intv.F
Chd|        HM_OPTION_IS_ENCRYPTED        source/devtools/hm_reader/hm_option_is_encrypted.F
Chd|        ELBUFTAG_MOD                  share/modules1/elbuftag_mod.F 
Chd|        MESSAGE_MOD                   share/message_module/message_mod.F
Chd|        SUBMODEL_MOD                  share/modules1/submodel_mod.F 
Chd|====================================================================
      SUBROUTINE HM_READ_MAT98(UPARAM   ,MATPARAM ,MAXUPARAM,NUPARAM  ,NUVAR    ,
     .                         MAXTABL  ,NUMTABL  ,ITABLE   ,MTAG     ,PARMAT   ,
     .                         PM       ,UNITAB   ,LSUBMODEL,MAT_ID   ,TITR     )
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE UNITAB_MOD
      USE ELBUFTAG_MOD            
      USE MESSAGE_MOD 
      USE SUBMODEL_MOD
      USE MATPARAM_DEF_MOD          
      USE TABLE_MOD
C-----------------------------------------------
C   ROUTINE DESCRIPTION :
C   ===================
C   READ MAT LAW98 WITH HM READER
C-----------------------------------------------
C   DUMMY ARGUMENTS DESCRIPTION:
C   ===================
C     UNITAB          UNITS ARRAY
C     MAT_ID          MATERIAL ID(INTEGER)
C     TITR            MATERIAL TITLE
C     LSUBMODEL       SUBMODEL STRUCTURE    
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "scr17_c.inc"
#include      "com04_c.inc"
#include      "units_c.inc"
#include      "submod_c.inc"
#include      "sysunit.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN)  :: MAT_ID,MAXTABL,MAXUPARAM
      INTEGER, INTENT(OUT) :: NUMTABL,NUPARAM,NUVAR
      my_real, DIMENSION(NPROPM) ,INTENT(INOUT)    :: PM     
      CHARACTER*nchartitle ,INTENT(IN)             :: TITR
      INTEGER, DIMENSION(MAXTABL)   ,INTENT(OUT)   :: ITABLE
      my_real, DIMENSION(MAXUPARAM) ,INTENT(OUT)   :: UPARAM
      my_real, DIMENSION(100),INTENT(OUT) :: PARMAT
      TYPE (UNIT_TYPE_),INTENT(IN) ::UNITAB 
      TYPE (SUBMODEL_DATA), DIMENSION(*),INTENT(IN) :: LSUBMODEL
      TYPE (MATPARAM_STRUCT_) ,INTENT(INOUT) :: MATPARAM
      TYPE (MLAW_TAG_), INTENT(INOUT) :: MTAG
      TYPE (TTABLE) TABLE(NTABLE)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      LOGICAL :: IS_AVAILABLE,IS_ENCRYPTED
      INTEGER :: I,ILAW,ISENS,ILOAD,ULOAD
      my_real :: RHO0,RHOR,YOUNG,EC,ET,BC,BT,G,G0,GT,GB,GSH,GFROT,
     .   KC,KT,KKC,KKT,KXC,KXT,KFC,KFT,FLEX,FLEX1,FLEX2,EMBC,EMBT,
     .   LC0,LT0,DC0,DT0,HC0,HT0,COSIN,TAN_LOCK,PHI_LOCK,
     .   VISCE,VISCG,AREAMIN1,AREAMIN2,ZEROSTRESS,STRESS_UNIT,FAC_M,FAC_T,FAC_L
      my_real ,DIMENSION(8) ::  YFAC
C=======================================================================
      IS_ENCRYPTED = .FALSE.
      IS_AVAILABLE = .FALSE.
      ILAW  = 98
      ILOAD = 0
      AREAMIN1 = ZERO

      NUMTABL = 8
      ITABLE(1:NUMTABL) = 0
      YFAC (1:NUMTABL) = ZERO
c
      CALL HM_OPTION_IS_ENCRYPTED(IS_ENCRYPTED)
c---------------------------------------------
card1
      CALL HM_GET_FLOATV('MAT_RHO'       ,RHO0      ,IS_AVAILABLE, LSUBMODEL, UNITAB)
      CALL HM_GET_FLOATV('Refer_Rho'     ,RHOR      ,IS_AVAILABLE, LSUBMODEL, UNITAB)
card2
c      
      CALL HM_GET_FLOATV('MAT_E1'        ,EC        ,IS_AVAILABLE, LSUBMODEL, UNITAB)
      CALL HM_GET_FLOATV('MAT_E2'        ,ET        ,IS_AVAILABLE, LSUBMODEL, UNITAB)
      CALL HM_GET_INTV  ('SENS_ID'       ,ISENS     ,IS_AVAILABLE,LSUBMODEL)
      CALL HM_GET_FLOATV('Zerostress'    ,ZEROSTRESS,IS_AVAILABLE, LSUBMODEL, UNITAB)
card3
      CALL HM_GET_INTV  ('FUN_BL1'       ,ITABLE(1)  ,IS_AVAILABLE,LSUBMODEL)
      CALL HM_GET_INTV  ('FUN_BL2'       ,ITABLE(2)  ,IS_AVAILABLE,LSUBMODEL)
      CALL HM_GET_INTV  ('FUN_UL1'       ,ITABLE(3)  ,IS_AVAILABLE,LSUBMODEL)
      CALL HM_GET_INTV  ('FUN_UL2'       ,ITABLE(4)  ,IS_AVAILABLE,LSUBMODEL)
      CALL HM_GET_INTV  ('FUN_SL'        ,ITABLE(5)  ,IS_AVAILABLE,LSUBMODEL)
card4
      CALL HM_GET_FLOATV('MAT_BL1'       ,YFAC(1)   ,IS_AVAILABLE, LSUBMODEL, UNITAB)
      CALL HM_GET_FLOATV('MAT_BL2'       ,YFAC(2)   ,IS_AVAILABLE, LSUBMODEL, UNITAB)
      CALL HM_GET_FLOATV('MAT_UL1'       ,YFAC(3)   ,IS_AVAILABLE, LSUBMODEL, UNITAB)
      CALL HM_GET_FLOATV('MAT_UL2'       ,YFAC(4)   ,IS_AVAILABLE, LSUBMODEL, UNITAB)
      CALL HM_GET_FLOATV('MAT_SL'        ,YFAC(5)   ,IS_AVAILABLE, LSUBMODEL, UNITAB)
c-----------------------------------------------------------------------
c     Check consistency of tabulated input data (loading and unloading)
c     a) there's no unloading functions => loading curves are optional
c        analytic and tabulated loading may be mixed
c     b) at least one unloading curve is defined => all loading corves must be defined
c        missing unloading curves may be created by Radioss by copying the loading ones
c-----------------------------------------------------------------------
c     Default values
c-----------------------------------------------------------------------
      CALL HM_GET_FLOATV_DIM('MAT_BL1',STRESS_UNIT ,IS_AVAILABLE, LSUBMODEL, UNITAB)
c
      DO I=1,8       
        IF (YFAC(I) == ZERO) YFAC(I) = ONE * STRESS_UNIT
      ENDDO
c-----------------------------------------------------------------------
c
      !LC0 = UN 
      !LT0 = UN 
      !DC0 = UN+EMBC
      !DT0 = UN+EMBT
      !HC0 = SQRT(DC0*DC0 - LC0*LC0)
      !HT0 = SQRT(DT0*DT0 - LT0*LT0)
c
c---  rigidite fil
      KC = EC
      KT = ET
c
      UPARAM( 1) = EC
      UPARAM( 2) = ET
          !UPARAM( 3) = HC0
          !UPARAM( 4) = HT0
          !UPARAM( 5) = KFC
          !UPARAM( 6) = KFT


      UPARAM(7) = ZEROSTRESS
      UPARAM(8) = ISENS
      UPARAM(9:8+NUMTABL) = YFAC(1:NUMTABL)

c
      NUVAR   = 39     
      NUPARAM = 46 !4pts pour l intersection +flag+PR SHEAR
C-----------------------
      YOUNG = MAX(EC,ET)
      PARMAT(1) = YOUNG*THIRD
      PARMAT(2) = YOUNG
      PARMAT(3) = ZERO
      PARMAT(4) = ZERO
      PARMAT(5) = ZERO
c--------------------------
      PM(1)  = RHOR
      PM(89) = RHO0
      PM(23) = YOUNG
c--------------------------
      MTAG%L_ANG  = 1
c--------------------------------------------------
c     Starter output
c--------------------------------------------------
      WRITE(IOUT,1000) TRIM(TITR),MAT_ID,ILAW
      WRITE(IOUT,1100)
      IF (IS_ENCRYPTED) THEN
        WRITE(IOUT,'(5X,A,//)')'CONFIDENTIAL DATA'
      ELSE
        WRITE(IOUT,1200) RHO0
        WRITE(IOUT,1300) EC,ET,ISENS,ZEROSTRESS,
     .                   ITABLE (1),ITABLE (2),ITABLE (3),ITABLE(4),
     .                   ITABLE(5),ITABLE(6),ITABLE(7),ITABLE(8),
     .                   YFAC(1),YFAC(2),YFAC(3),YFAC(4),
     .                   YFAC(5),YFAC(6),YFAC(7),YFAC(8)




c        WRITE(IOUT,1300) EC,ET,ISENS,ZEROSTRESS
c        WRITE(IOUT,1400) ITABLE(1),ITABLE(2),ITABLE(3),ITABLE(4)
c        WRITE(IOUT,1500) ITABLE(1),ITABLE(2),ITABLE(3),YFAC(1),YFAC(2),YFAC(3)




c        IF (ILOAD == 0) THEN
c          WRITE(IOUT,1300) BC,BT,G0,GT,PHI_LOCK*HUNDRED80/PI
c        ELSE
c          WRITE(IOUT,1400) ITABLE(1),ITABLE(2),ITABLE(3),YFAC(1),YFAC(2),YFAC(3)
c          IF (ILOAD == 2)
c     .    WRITE(IOUT,1500) ITABLE(4),ITABLE(5),ITABLE(6),YFAC(4),YFAC(5),YFAC(6)
c        ENDIF
c        WRITE(IOUT,1600) VISCE,VISCG,GFROT,GSH,ZEROSTRESS,
c     .        EMBC,EMBT,NC,NT,ISENS,FLEX1,FLEX2
      ENDIF     
c-----------------------------------------------------------------------
 1000 FORMAT(/
     & 5X,A,/,
     & 5X,'MATERIAL NUMBER. . . . . . . . . . . . . . . . .  .=',I10/,
     & 5X,'MATERIAL LAW . . . . . . . . . . . . . . . . . .  .=',I10/)
 1100 FORMAT
     &(5X,'MATERIAL MODEL : ANISOTROPIC FABRIC (LAW98)       ',/,
     & 5X,'--------------------------------------------      ',/)
 1200 FORMAT(
     & 5X,'INITIAL DENSITY . . . .  .  . . . . . . . . . . . .=',1PG20.13)  
 1300 FORMAT(
     & 5X,'YOUNG MODULUS E1 (WARP DIRECTION) . . . . . . . . .=',1PG20.13/
     & 5X,'YOUNG MODULUS E2 (WEFT DIRECTION) . . . . . . . . .=',1PG20.13/
     & 5X,'SENSOR ID . . . . . . . . . . . . . . . . . . . . .=',I10/
     & 5X,'REF-STATE STRESS RELAXATION FACTOR. . . . . . . . .=',1PG20.13/
     & 5X,'UNIAXIAL LOADING FUNCTION ID IN WARP DIRECTION. . . .=',I10/
     & 5X,'UNIAXIAL LOADING FUNCTION ID IN WEFT DIRECTION. . . .=',I10/ 
     & 5X,'SHEAR LOADING FUNCTION ID IN SHEAR. . . . . . . . . .=',I10/
     & 5X,'UNIAXIAL UNLOADING FUNCTION ID IN WARP DIRECTION. . .=',I10/
     & 5X,'UNIAXIAL UNLOADING FUNCTION ID IN WEFT DIRECTION. . .=',I10/ 
     & 5X,'SHEAR UNLOADING FUNCTION ID IN SHEAR. . . . . . . . .=',I10/
     & 5X,'BIAXIAL LOADING FUNCTION ID IN WARP DIRECTION . . . .=',I10/
     & 5X,'BIAXIAL LOADING FUNCTION ID IN WEFT DIRECTION . . . .=',I10/ 
     & 5X,'SCALE FACTOR UNIAXIAL LOADING(WARP) . . . . . . . . .=',1PG20.13/
     & 5X,'SCALE FACTOR UNIAXIAL LOADING(WEFT) . . . . . . . . .=',1PG20.13/
     & 5X,'SCALE FACTOR LOADING(SHEAR) . . . . . . . . . . . . .=',1PG20.13/ 
     & 5X,'SCALE FACTOR UNIAXIAL UNLOADING(WARP) . . . . . . . .=',1PG20.13/
     & 5X,'SCALE FACTOR UNIAXIAL UNLOADING(WEFT) . . . . . . . .=',1PG20.13/
     & 5X,'SCALE FACTOR UNLOADING(SHEAR) . . . . . . . . . . . .=',1PG20.13/ 
     & 5X,'SCALE FACTOR BIAXIAL LOADING(WARP). . . . . . . . . .=',1PG20.13/
     & 5X,'SCALE FACTOR BIAXIAL LOADING(WEFT). . . . . . . . . .=',1PG20.13/)


c 1400 FORMAT(
c     & 5X,'UNIAXIAL LOADING FUNCTION ID IN WARP DIRECTION. . .=',I10/
c     & 5X,'UNIAXIAL LOADING FUNCTION ID IN WEFT DIRECTION. . .=',I10/ 
c     & 5X,'BIAXIAL  LOADING FUNCTION ID IN WARP DIRECTION. . .=',I10/
c     & 5X,'BIAXIAL  LOADING FUNCTION ID IN WEFT DIRECTION. . .=',I10/)
c 1500 FORMAT(
c     & 5X,'UNLOADING STRESS FUNCTION ID IN WARP DIRECTION. . .=',I10/
c     & 5X,'UNLOADING STRESS FUNCTION ID IN WEFT DIRECTION. . .=',I10/ 
c     & 5X,'UNLOADING STRESS FUNCTION ID IN SHEAR DIRECTION . .=',I10/
c     & 5X,'UNIAXIAL LOADING SCALE FACTOR (WARP). . . . . . . .=',1PG20.13/
c     & 5X,'UNIAXIAL LOADING SCALE FACTOR (WEFT). . . . . . . .=',1PG20.13/
c     & 5X,'BIAXIAL  LOADING SCALE FACTOR (WARP). . . . . . . .=',1PG20.13/
c     & 5X,'BIAXIAL  LOADING SCALE FACTOR (WEFT). . . . . . . .=',1PG20.13/)
c-----------------------------------------------------------------------
      RETURN
      END
