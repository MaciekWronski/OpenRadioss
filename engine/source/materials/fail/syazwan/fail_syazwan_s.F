Copyright>        OpenRadioss
Copyright>        Copyright (C) 2026 Siemens
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
Copyright>        Commercial Alternative: Simcenter Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Siemens also offers Simcenter(TM) Radioss(R)
Copyright>        software under a commercial license.  Contact Siemens to discuss further if the
Copyright>        commercial version may interest you: 
Copyright>        https://www.siemens.com/en-us/products/simcenter/mechanical-simulation/radioss/.
!||====================================================================
!||    fail_syazwan_s   ../engine/source/materials/fail/syazwan/fail_syazwan_s.F
!||--- called by ------------------------------------------------------
!||    mmain            ../engine/source/materials/mat_share/mmain.F90
!||    mulaw            ../engine/source/materials/mat_share/mulaw.F90
!||    usermat_solid    ../engine/source/materials/mat_share/usermat_solid.F
!||--- calls      -----------------------------------------------------
!||    finter           ../engine/source/tools/curve/finter.F
!||====================================================================
      SUBROUTINE FAIL_SYAZWAN_S(FAIL  ,
     1     NEL     ,NUVAR    ,NVARTMP ,UVAR    ,VARTMP  ,
     2     TIME    ,NGL      ,IPG     ,DPLA    ,TDELE   ,
     3     SIGNXX  ,SIGNYY   ,SIGNZZ  ,SIGNXY  ,SIGNYZ  ,SIGNZX  ,
     4     DFMAX   ,ALDT     ,OFF     ,UELR    ,NPG     ,LOFF    )
C-----------------------------------------------
C    Modules
C-----------------------------------------------
      use fail_param_mod
      use table_mat_vinterp_mod
C!-----------------------------------------------
C!   I m p l i c i t   T y p e s
C!-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   G l o b a l   P a r a m e t e r s
C-----------------------------------------------
#include      "units_c.inc"
#include      "comlock.inc"
C-----------------------------------------------
C   I N P U T   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: NEL,NVARTMP,NUVAR,NPG,IPG
      INTEGER :: NGL(NEL)
      my_real, INTENT(IN) :: TIME
      my_real, INTENT(IN) :: DPLA(NEL)  ,ALDT(NEL)
      my_real, INTENT(IN) :: SIGNXX(NEL),SIGNYY(NEL),SIGNZZ(NEL),SIGNXY(NEL),SIGNYZ(NEL),SIGNZX(NEL)
      type (fail_param_) ,intent(inout) :: fail    !< failure model data structure
C-----------------------------------------------
C   I N P U T   O U T P U T   A r g u m e n t s 
C-----------------------------------------------
      integer ,dimension(nel,nvartmp) ,intent(inout) :: vartmp
      my_real ,dimension(nel,nuvar)   ,intent(inout) :: uvar
      my_real, INTENT(INOUT) :: LOFF(NEL),OFF(NEL),DFMAX(NEL),TDELE(NEL),UELR(NEL)
C!-----------------------------------------------
C!   L o c a l   V a r i a b l e s
C!-----------------------------------------------
      INTEGER I,J,NINDX,NINDX2,FAILIP
      INTEGER, DIMENSION(NEL) :: INDX,INDX2
      my_real :: C1,C2,C3,C4,C5,C6,REF_LEN  
      my_real :: P,SVM,TRIAX,COS3THETA,LODEP,EPSFAIL,DET,SXX,SYY,SZZ,EPFMIN
      my_real ,dimension(nel) :: fac,lambda, dydx
      !=======================================================================
      ! - INITIALISATION OF COMPUTATION ON TIME STEP
      !=======================================================================
      ! Recovering failure criterion parameters
      !
      C1        = FAIL%UPARAM(1) 
      C2        = FAIL%UPARAM(2) 
      C3        = FAIL%UPARAM(3)
      C4        = FAIL%UPARAM(4) 
      C5        = FAIL%UPARAM(5)
      C6        = FAIL%UPARAM(6)
      REF_LEN   = FAIL%UPARAM(11) 
      EPFMIN    = FAIL%UPARAM(12)

      FAILIP    = MIN(FAIL%IPARAM(4),NPG)
c
      ! At first timestep, initialization of the critical damage and 
      ! the element size scaling factor
      IF (UVAR(1,1) == ZERO) THEN
        IF (fail%table4d(1)%notable > 0) THEN 
          LAMBDA(1:nel) = ALDT(1:nel)/REF_LEN
          call table_mat_vinterp(fail%table4d(1),nel,nel,vartmp,lambda,uvar(1:nel,1:1),dydx)
        ELSE 
          UVAR(1:NEL,1) = ONE
        ENDIF
      ENDIF
c
      DO I=1,NEL
        ! Recover element size scaling
        FAC(I) = UVAR(I,1)
        IF (OFF(I) < EM01) OFF(I) = ZERO
        IF (OFF(I) < ONE .AND. OFF(I) > ZERO) OFF(I) = OFF(I)*FOUR_OVER_5
      ENDDO
c
      !====================================================================
      ! - COMPUTATION OF THE DAMAGE VARIABLE EVOLUTION
      !==================================================================== 
      ! Initialization of element failure index
      NINDX = 0  
      INDX(1:NEL) = 0
      NINDX2 = 0
      INDX2(1:NEL) = 0
c
      ! Loop over the elements 
      DO I=1,NEL
c
        IF (LOFF(I) == ONE .AND. DPLA(I) > ZERO .AND. OFF(I) == ONE) THEN
c
          ! Computation of hydrostatic stress, Von Mises stress, and stress triaxiality
          P   = THIRD*(SIGNXX(I) + SIGNYY(I) + SIGNZZ(I))     
          SXX = SIGNXX(I) - P                                 
          SYY = SIGNYY(I) - P                                 
          SZZ = SIGNZZ(I) - P                                 
          SVM = HALF*(SXX**2 + SYY**2 + SZZ**2)             
     .         + SIGNXY(I)**2+ SIGNZX(I)**2 + SIGNYZ(I)**2    
          SVM = SQRT(THREE*SVM) 
          TRIAX = P/MAX(EM20,SVM)
          IF (TRIAX < -TWO_THIRD) TRIAX = -TWO_THIRD
          IF (TRIAX >  TWO_THIRD) TRIAX =  TWO_THIRD
c
          ! Computation of Lode parameter
          DET = SXX*SYY*SZZ + TWO*SIGNXY(I)*SIGNZX(I)*SIGNYZ(I)- 
     .          SXX*SIGNYZ(I)**2-SZZ*SIGNXY(I)**2-SYY*SIGNZX(I)**2
          COS3THETA = HALF*TWENTY7*DET/MAX(EM20,SVM**3)
          IF(COS3THETA < -ONE) COS3THETA = -ONE
          IF(COS3THETA > ONE)  COS3THETA = ONE
          LODEP = ONE - TWO*ACOS(COS3THETA)/PI
c
          ! Computation of the plastic strain at failure
          EPSFAIL = C1 + C2*TRIAX + C3*LODEP + C4*(TRIAX**2) + 
     .              C5*(LODEP**2) + C6*TRIAX*LODEP
          EPSFAIL = EPSFAIL*FAC(I)
          EPSFAIL = MAX(EPFMIN,EPSFAIL)
c
          ! Computation of the damage variable update 
          DFMAX(I) = DFMAX(I) + DPLA(I)/MAX(EPSFAIL,EM20)
          DFMAX(I) = MIN(DFMAX(I),ONE) 
c              
          ! Check element failure    
          IF (DFMAX(I) >= ONE) THEN
            LOFF(I)     = ZERO
            NINDX       = NINDX + 1
            INDX(NINDX) = I
            UELR(I)     = UELR(I) + ONE
            IF (NINT(UELR(I)) >= FAILIP) THEN 
              OFF(I)    = FOUR_OVER_5
              TDELE(I)  = TIME
              NINDX2    = NINDX2 + 1 
              INDX2(NINDX2) = I
            ENDIF
          ENDIF
        ENDIF
      ENDDO
c
      !====================================================================
      ! - PRINTOUT DATA ABOUT FAILED ELEMENTS
      !====================================================================
      IF (NINDX>0) THEN
        DO J=1,NINDX
          I = INDX(J)     
#include "lockon.inc"
          WRITE(IOUT, 1000) NGL(I),IPG,TIME
          WRITE(ISTDO,1000) NGL(I),IPG,TIME
#include "lockoff.inc"
        END DO
      ENDIF    
      IF (NINDX2>0) THEN
        DO J=1,NINDX2
          I = INDX2(J)     
#include "lockon.inc"
          WRITE(IOUT, 2000) NGL(I),TIME
          WRITE(ISTDO,2000) NGL(I),TIME
#include "lockoff.inc"
        END DO
      ENDIF       
C------------------
 1000 FORMAT(1X,'FOR SOLID ELEMENT NUMBER el#',I10,
     .          ' FAILURE (SYAZWAN) AT GAUSS POINT ',I5,
     .          ' AT TIME :',1PE12.4)
 2000 FORMAT(1X,'-- RUPTURE OF SOLID ELEMENT :',I10,
     .          ' AT TIME :',1PE12.4)     
      END
