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
Chd|  LAW98_IDENTPAR                      source/materials/mat/mat098/law98_upd.F
Chd|-- called by -----------
Chd|        UPDMAT                        source/materials/updmat.F     
Chd|-- calls ---------------
Chd|        ANCMSG                        source/output/message/message.F
Chd|        FUNC_INTERS                   source/tools/curve/func_inters.F
Chd|        FUNC_INTERS_SHEAR             source/tools/curve/func_inters.F
Chd|        FUNC_SLOPE                    source/tools/curve/func_slope.F
Chd|        MATFUN_USR2SYS                source/materials/tools/matfun_usr2sys.F
Chd|        MESSAGE_MOD                   share/message_module/message_mod.F
Chd|        SENSOR_MOD                    share/modules1/sensor_mod.F   
Chd|        TABLE_MOD                     share/modules1/table_mod.F    
Chd|====================================================================
      SUBROUTINE LAW98_IDENTPAR(IOUT      ,TITR      ,MAT_ID   ,NUPARAM  ,MATPARAM ,
     .                       UPARAM    ,NUMTABL   ,ITABLE   ,TABLE    ,PM       ,
     .                       SENSOR_TAB)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
      USE MATPARAM_DEF_MOD
      USE TABLE_MOD
      USE SENSOR_MOD
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
#include      "tabsiz_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      CHARACTER*nchartitle  :: TITR
      INTEGER :: MAT_ID,IOUT,NUMTABL,NUPARAM,NPTS
      INTEGER ,DIMENSION(NUMTABL) :: ITABLE
      my_real ,DIMENSION(NUPARAM) :: UPARAM
      my_real ,DIMENSION(NPROPM)  :: PM
      TYPE (SENSOR_STR_) ,DIMENSION(NSENSOR) :: SENSOR_TAB
      TYPE (MATPARAM_STRUCT_)  ,TARGET :: MATPARAM
      TYPE (TTABLE), DIMENSION(NTABLE) ,TARGET ::  TABLE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I,K,FUNC,FUND,UNLOAD,PN,IOK,ISENS,SENS_ID,TEST,NVAR,
     .           NITER,NDIM,IVAR
      INTEGER :: FUNC_B1,FUNC_B2,FUNC_U1,FUNC_U2,FUNC_SH,INFO
      my_real KC,KT,KCMAX,KTMAX,KFC,KFT,GMAX,STIFFMIN,STIFFMAX,STIFFINI,STIFFAVG
c
      INTEGER ID_F1,ID_F2
      my_real :: XINT,YINT,XINT1,YINT1,XINT2 ,YINT2,NU,LAMBDA,
     .           YFAC_B1,YFAC_B2,YFAC_U1,YFAC_U2,YFAC_SH,
     .           DC0,DT0,HC0,HT0,XINTC,YINTC
      my_real :: STIFF,STIFFC,STIFFT,STFAVG,WEIGHT1,WEIGHT2,RI2,RI2_OLD,EC,ET,
     .           DENOM,RHO,JACTDY_N2
      my_real :: STR1,STR2,FLEX1,FLEX2
      my_real, DIMENSION(:)   ,ALLOCATABLE :: DERI,FAC,DERIINI,XUNIAX,XBIAX,
     .   YBIAX1,YBIAX2,YUNIAX1,YUNIAX2,SLOPE1,SLOPE2
      my_real, DIMENSION(:)   ,ALLOCATABLE :: XX,XXTR,TTMP
      my_real, DIMENSION(:,:) ,ALLOCATABLE :: JAC_U1,JAC_U2,JAC_U1TR,JAC_U2TR,
     .                                        JAC_U1T,JAC_U2T,JTJ,WW

      my_real, DIMENSION(:) ,ALLOCATABLE :: XFIB1,XFIB2,YFIB1,YFIB2,FSQR1,FSQR2,
     .                                      FSQR1TR,FSQR2TR,JACTDY,VECT,MU,XV
      INTEGER, DIMENSION(:)   ,ALLOCATABLE :: IPIV
      INTEGER, DIMENSION(:,:) ,ALLOCATABLE :: IPOS

      INTEGER ICONTR,LWORK,MAXFUN,MAXLIN,NCON,NEQ,NDC,IGOTO,ICONV,ISYM,
     .   ITER,PN1,PN2,SIZEPN,IDN,IDN1,IDN2,NPOINT,NPTU,NPTB,LEN1,LEN2

      TYPE(TTABLE), DIMENSION(:) ,POINTER ::  TABLE_MAT
      INTEGER, PARAMETER :: MXITER = 10000
      LOGICAL :: CONV   
      my_real :: EPS1 
      my_real :: EPS2 
      my_real :: ALPHA
C=======================================================================
C     SENSOR NUMBERING CHECK
C----------------------------
      SENS_ID = UPARAM(8)
      ISENS   = 0
      IF (SENS_ID > 0 ) THEN
        DO I=1,NSENSOR
          IF (SENS_ID == SENSOR_TAB(I)%SENS_ID) THEN
            ISENS = I
            UPARAM(8) = ISENS
            EXIT
          END IF
        ENDDO
        IF (ISENS == 0)  
     .    CALL ANCMSG(MSGID=1240,ANMODE=ANINFO,MSGTYPE=MSGWARNING,
     .  	      I1=MAT_ID,C1=TITR,I2=ISENS) 
      ENDIF
c---------------------------------------------------------------
      FUNC_B1 = ITABLE(1)
      FUNC_B2 = ITABLE(2)
      FUNC_U1 = ITABLE(3)
      FUNC_U2 = ITABLE(4)
      FUNC_SH = ITABLE(5)
c
      YFAC_B1 = UPARAM(9)
      YFAC_B2 = UPARAM(10)
      YFAC_U1 = UPARAM(11)
      YFAC_U2 = UPARAM(12)
      YFAC_SH = UPARAM(13)
c---------------------------------------------------------------

!      !DIRECTION CHAINE
!      ID_F1 = IFUNC(1)
!      ID_F2 = IFUNC(4)
!
!      FAC1 = UPARAM(8+1)
!      FAC2 = UPARAM(8+4)
!      CALL FUNC_INTERS(TITR,MAT_ID,ID_F1,ID_F2,FAC1,FAC2 ,NPC,PLD ,XINTC,YINTC )
!      UPARAM(36) = XINTC
!      UPARAM(37) = YINTC
!      !PRINT*, 'DIRECTION CHAINE XINT YINT',XINTC, YINTC
!c
!      !DIRECTION TRAME
!      ID_F1 = IFUNC(2)
!      ID_F2 = IFUNC(5)
!      FAC1 = UPARAM(8+2)
!      FAC2 = UPARAM(8+5)   
!      CALL FUNC_INTERS(TITR,MAT_ID,ID_F1,ID_F2,FAC1,FAC2 ,NPC,PLD ,XINT,YINT )
!      UPARAM(38) = XINT
!      UPARAM(39) = YINT
!      !PRINT*, 'DIRECTION TRAME  XINT YINT',XINT, YINT
!
!      !DIRECTION SHEAR
!      ID_F1 = IFUNC(3)
!      ID_F2 = IFUNC(6)
!      FAC1 = UPARAM(8+3)
!      FAC2 = UPARAM(8+6)   
!      CALL FUNC_INTERS_SHEAR(TITR ,MAT_ID,ID_F1,ID_F2 ,FAC1  ,FAC2   ,
!     .                       NPC  ,PLD   ,XINT1 ,YINT1  ,XINT2 ,YINT2)
!      UPARAM(43) = XINT1
!      UPARAM(44) = YINT1
!
!      UPARAM(45) = XINT2
!      UPARAM(46) = YINT2
      !PRINT*, 'DIRECTION SHEAR  XINT1 YINT1',XINT1, YINT1
      !PRINT*, '                 XINT2 YINT2',XINT2, YINT2


       
       
C====================================================================
c       MAXIMUM FUNCTION SLOPE -> STIFFNESS for time step computation


        CALL TABLE_SLOPE1(TABLE(FUNC_B1),YFAC_B1,STIFFMIN,STIFFMAX,STIFFINI,STIFFAVG)
        
        EC = STIFFMAX

        CALL TABLE_SLOPE1(TABLE(FUNC_B2),YFAC_B2,STIFFMIN,STIFFMAX,STIFFINI,STIFFAVG)
        
        ET = STIFFMAX

        UPARAM(1) = EC
        UPARAM(2) = ET



!       STIFFC  = MAX(DERI(1),DERI(7))
!       STIFFT  = MAX(DERI(2),DERI(8))
!       STIFF   = MAX(STIFFC,STIFFT)
!       !PRINT*, ' STIFFC, STIFFT ', STIFFC, STIFFT
!       PM(20)= STIFF
!       PM(22)= STIFF*HALF
!       PM(24)= STIFF
!       PM(32)= STIFF/THREE
!       !UPARAM( 1) = STIFFC 
!       !UPARAM( 2) = STIFFT 
!
!       UPARAM(40)  = DERIINI(1)
!       UPARAM(41)  = DERIINI(2)
!C====================================================================

c     unify abscissa vectors for both UNIAXIAL curves ! not needed but to uniformize data
      LEN1 = SIZE(TABLE(FUNC_U1)%X(1)%VALUES)
      LEN2 = SIZE(TABLE(FUNC_U2)%X(1)%VALUES)
      NPOINT = LEN1 +  LEN2
      ALLOCATE(XUNIAX(NPOINT))
      ALLOCATE(XBIAX(NPOINT))
      ALLOCATE(SLOPE1(NPOINT))
      ALLOCATE(SLOPE2(NPOINT))

      CALL UNIFY_2TAB_X(FUNC_U1,FUNC_U2,TABLE,LEN1,LEN2,NPOINT,XUNIAX)
      
      NPTU = NPOINT
      NDIM = TABLE(FUNC_U1)%NDIM

      ALLOCATE(YUNIAX1(NPTU))
      ALLOCATE(YUNIAX2(NPTU))
      
c   
c     Interpolate data with new abscissa vector

      DO I=1,NPTU
        CALL TABLE_INTERP1(TABLE(FUNC_U1),NDIM,XUNIAX(I),YUNIAX1(I))
        CALL TABLE_INTERP1(TABLE(FUNC_U2),NDIM,XUNIAX(I),YUNIAX2(I))
      END DO
c
c-----------------------------------------------------------------
c     unify abscissa vectors for both BIAXIAL curves
      LEN1 = SIZE(TABLE(FUNC_B1)%X(1)%VALUES)
      LEN2 = SIZE(TABLE(FUNC_B2)%X(1)%VALUES)
      NPOINT = LEN1 +  LEN2

      CALL UNIFY_2TAB_X(FUNC_B1,FUNC_B1,TABLE,LEN1,LEN2,NPOINT,XBIAX)

      NPTB = NPOINT
      ALLOCATE(YBIAX1(NPTB))
      ALLOCATE(YBIAX2(NPTB))
c   
c     Interpolate data with new abscissa vector

      DO I=1,NPTB
        CALL TABLE_INTERP1(TABLE(FUNC_B1),NDIM,XBIAX(I),YBIAX1(I))
        CALL TABLE_INTERP1(TABLE(FUNC_B2),NDIM,XBIAX(I),YBIAX2(I))
      END DO
c     
      YUNIAX1(1:NPTB) = YUNIAX1(1:NPTU) * YFAC_U1
      YUNIAX2(1:NPTB) = YUNIAX2(1:NPTU) * YFAC_U2
      YBIAX1 (1:NPTB) = YBIAX1 (1:NPTB) * YFAC_B1
      YBIAX2 (1:NPTB) = YBIAX2 (1:NPTB) * YFAC_B2
c-----------------------------------------------------------------
      
      NVAR = 4
      ALLOCATE(XFIB1(NPTB))
      ALLOCATE(XFIB2(NPTB))
      ALLOCATE(YFIB1(NPTB))
      ALLOCATE(YFIB2(NPTB))

      ALLOCATE(FSQR1(NPTU))
      ALLOCATE(FSQR2(NPTU))
      ALLOCATE(FSQR1TR(NPTU))
      ALLOCATE(FSQR2TR(NPTU))
      ALLOCATE(JAC_U1(NPTU,NVAR))
      ALLOCATE(JAC_U2(NPTU,NVAR))
      ALLOCATE(JAC_U1T(NVAR,NPTU))
      ALLOCATE(JAC_U2T(NVAR,NPTU))
      ALLOCATE(JAC_U1TR(NPTU,NVAR))
      ALLOCATE(JAC_U2TR(NPTU,NVAR))
      ALLOCATE(JTJ(NVAR,NVAR))
      ALLOCATE(JACTDY(NVAR))  ! gradient of residues
      ALLOCATE(VECT(NVAR))
      ALLOCATE(MU(NVAR))
      ALLOCATE(IPIV(NVAR))
      ALLOCATE(XX(NVAR))
      ALLOCATE(XXTR(NVAR))
      ALLOCATE(WW(NPTU,NPTU))
      ALLOCATE(TTMP(NPTU))
      
c     design variables 
c          XX(1) = stretch1  -> HC0,DC0
c          XX(2) = stretch2  -> HT0,DT0
c          XX(3) = flex1
c          XX(4) = flex2

      XX(1) = 0.1
      XX(2) = 0.1
      XX(3) = 0.2*EC
      XX(4) = 0.2*ET

      STR1  = XX(1)
      STR2  = XX(2)
      FLEX1 = XX(3)
      FLEX2 = XX(4)
      WEIGHT1 = HALF
      WEIGHT2 = HALF
      EPS1  = EM08
      EPS2  = EM10
      ALPHA = HALF
c---------------------------------------------
c     Initialize target function + jacobian
c---------------------------------------------
         
c     Calculate fiber stress-strain from biaxial tension 

      CALL FIBER_FUNC(NPTB  ,STR1   ,STR2   ,FLEX1   ,FLEX2   ,
     .                XBIAX ,YBIAX1 ,YBIAX2 ,XFIB1   ,YFIB1   ,
     .                XFIB2 ,YFIB2  )              
      
c     uniaxial tension fiber 1
      IVAR = 1
      CALL JAC98_UNI(NVAR   ,NPTU  ,NPTB  ,XUNIAX,YUNIAX1,
     .               XFIB1  ,YFIB1 ,STR1  ,FLEX1   ,FLEX2   ,
     .               WEIGHT1,FSQR1 ,JAC_U1,IVAR    )
c     uniaxial tension fiber 2
      IVAR = 2
      CALL JAC98_UNI(NVAR   ,NPTU  ,NPTB  ,XUNIAX,YUNIAX2,
     .               XFIB2  ,YFIB2 ,STR2  ,FLEX1   ,FLEX2   ,
     .               WEIGHT2,FSQR2 ,JAC_U2,IVAR    )
                     
c     Computation of the squared sum of residues
      RI2 = DOT_PRODUCT(FSQR1,FSQR1) + DOT_PRODUCT(FSQR2,FSQR2)

      ! Initialization of the damping parameter and the convergence criterion
      LAMBDA = MAXVAL(JAC_U1) + MAXVAL(JAC_U2)
      NU     = TWO
      CONV   = .FALSE.
      NITER  = 0

      WW(:,:) = ZERO
      DO I = 1,NVAR
        WW(I,I) = ONE
      ENDDO      
c----------------------------------------------------------------------
c       Levenberg-Marquardt Algorithm + Line-search
c----------------------------------------------------------------------
ccc#ifndef WITHOUT_LINALG

      DO WHILE ((.NOT.CONV).and.(NITER < MXITER))
C
        ! Transpose of Jacobian matrix
        JAC_U1T = TRANSPOSE(JAC_U1)
        JAC_U2T = TRANSPOSE(JAC_U2)
        JTJ     = MATMUL(JAC_U1T,JAC_U1) + MATMUL(JAC_U2T,JAC_U2)
        !  -> Adding the damping parameter
        DO I = 1,NVAR
          JTJ(I,I) = JTJ(I,I) + LAMBDA
        ENDDO
                     
        !  -> Compute the gradient of residues 
        
c        TTMP = MATMUL(WW,FSQR1)
        JACTDY = -MATMUL(JAC_U1T,FSQR1) - MATMUL(JAC_U2T,FSQR2)

c        JACTDY = -MATMUL(JAC_U1T,MATMUL(WW,FSQR1)) - MATMUL(JAC_U2T,MATMUL(WW,FSQR2))
      
c
        ! b) Compute the trial values of XI
        CALL DGESV(NVAR, 1, JTJ, NVAR, IPIV, JACTDY , NVAR, INFO)
        XXTR(1:NVAR) = XX(1:NVAR) + JACTDY(1:NVAR)

        ! Checking convergence
        IF (SQRT(DOT_PRODUCT(JACTDY,JACTDY)) <= EPS2*(SQRT(DOT_PRODUCT(XX,XX))+EPS2)) THEN
          CONV = .TRUE.

        ELSE
        
          ! c) Compute the trial values of the function
          STR1  = MAX(ZERO,XXTR(1))
          STR2  = MAX(ZERO,XXTR(2))
          FLEX1 = MAX(ZERO,XXTR(3))
          FLEX2 = MAX(ZERO,XXTR(4))
          
          CALL FIBER_FUNC(NPTB  ,STR1   ,STR2   ,FLEX1   ,FLEX2   ,
     .                    XBIAX ,YBIAX1 ,YBIAX2 ,XFIB1   ,YFIB1   ,
     .                    XFIB2 ,YFIB2  )              
      
c         uniaxial tension fiber 1
          IVAR = 1
          CALL JAC98_UNI(NVAR   ,NPTU  ,NPTB  ,XUNIAX,YUNIAX1,
     .                   XFIB1  ,YFIB1 ,STR1  ,FLEX1   ,FLEX2   ,
     .                   WEIGHT1,FSQR1TR  ,JAC_U1TR,IVAR    )
c         uniaxial tension fiber 2
          IVAR = 2
          CALL JAC98_UNI(NVAR   ,NPTU  ,NPTB  ,XUNIAX,YUNIAX2,
     .                   XFIB2  ,YFIB2 ,STR2  ,FLEX1   ,FLEX2   ,
     .                   WEIGHT2,FSQR2TR  ,JAC_U2TR,IVAR    )
          
          
          ! d) Computation of the squared sum of residues
          RI2_OLD = RI2
          RI2     = DOT_PRODUCT(FSQR1TR,FSQR1TR) + DOT_PRODUCT(FSQR2TR,FSQR2TR)
C
          ! Compute gain ratio
          VECT  = MATMUL(JAC_U1T,FSQR1) + MATMUL(JAC_U2T,FSQR2)
          VECT  = LAMBDA*JACTDY - VECT
          DENOM = DOT_PRODUCT(JACTDY,VECT)
          RHO   = (RI2_OLD - RI2) / MAX(DENOM,EM20)
C          
          ! e) Update parameters
          !  -> Step accepted
          IF (RHO > ZERO .and. INFO == 0) THEN
            ! Update damping parameter
            LAMBDA = LAMBDA*MAX(THIRD,ONE-(TWO*RHO-ONE)**3)
            NU  = TWO
            ! Update XI vector
            !  -> Compute the MU vector for line search
            DO I = 1,NVAR
              IF (JACTDY(I) >= ZERO) THEN
                MU(I) = ONE
              ELSE
                MU(I) = MIN(ONE,((ONE-ALPHA)/(-JACTDY(I)))*XX(I))
              ENDIF
              XX(I) = XX(I) + MU(I)*JACTDY(I)
            ENDDO  
            ! Update FXI vector and Jacobian matrix
          
            STR1  = MAX(ZERO,XX(1))
            STR2  = MAX(ZERO,XX(2))
            FLEX1 = MAX(ZERO,XX(3))
            FLEX2 = MAX(ZERO,XX(4))
      
          
            CALL FIBER_FUNC(NPTB  ,STR1   ,STR2   ,FLEX1   ,FLEX2   ,
     .                      XBIAX ,YBIAX1 ,YBIAX2 ,XFIB1   ,YFIB1   ,
     .                      XFIB2 ,YFIB2  )              
      
c           uniaxial tension fiber 1
            IVAR = 1
            CALL JAC98_UNI(NVAR    ,NPTU     ,NPTB    ,XUNIAX  ,YUNIAX1 ,
     .                     XFIB1   ,YFIB1    ,STR1    ,FLEX1   ,FLEX2   ,
     .                     WEIGHT1 ,FSQR1    ,JAC_U1  ,IVAR    )
c           uniaxial tension fiber 2
            IVAR = 2
            CALL JAC98_UNI(NVAR    ,NPTU     ,NPTB    ,XUNIAX  ,YUNIAX2 ,
     .                     XFIB2   ,YFIB2    ,STR2    ,FLEX1   ,FLEX2   ,
     .                     WEIGHT2 ,FSQR2    ,JAC_U2  ,IVAR    )

            ! Update residue vector
            RI2 = DOT_PRODUCT(FSQR1,FSQR1) + DOT_PRODUCT(FSQR2,FSQR2)

            JAC_U1T = TRANSPOSE(JAC_U1)
            JAC_U2T = TRANSPOSE(JAC_U2)
          
          ELSE   !  -> Step not accepted
                   
            LAMBDA = LAMBDA*NU  ! Update damping parameter    
            NU     = NU*TWO

          ENDIF          
C          
          ! f) Computing the convergence criterion
          !  -> Compute the gradient of residues 
          JACTDY = TWO*MATMUL(JAC_U1T,FSQR1) 
     .           + TWO*MATMUL(JAC_U2T,FSQR2)
          JACTDY_N2 = SQRT(DOT_PRODUCT(JACTDY,JACTDY))
          IF (JACTDY_N2 < EPS1) CONV = .TRUE.
c
          ! g) Increasing the iteration number
          NITER = NITER + 1
c
        ENDIF ! CONVERGED
C
      ENDDO  ! WHILE NOT CONVERGED
C
c------------------------------
c     end of optimization loop and check of results
c------------------------------

      ! Algorithm has not converged
      IF (NITER == MXITER .OR. .NOT.CONV) THEN
        CALL ANCMSG(MSGID=1926,MSGTYPE=MSGERROR,
     .              ANMODE=ANINFO_BLIND_1,I1=MAT_ID)
      ! Solutions values must be checked
      ELSEIF (ABS(JACTDY_N2)>EPS1 .OR. MAXVAL(XX)>EP10 .OR. MINVAL(XX)<EM10) THEN
        CALL ANCMSG(MSGID=1927,MSGTYPE=MSGWARNING,
     .              ANMODE=ANINFO_BLIND_1,I1=MAT_ID)   
      ENDIF
c------------------------------
      ! Storage of the results
c
c-----------------------------------------------------------------------
 1000 FORMAT(
     & 5X,40H  MATERIAL CALCULATED DATA LAW98                  ,/,
     & 5X,40H  -------------------------------                 ,//)
 1100 FORMAT(
     & 5X,'NOMINAL STRETCH . (WARP DIRECTION). . . . .=',1PG20.13/
     & 5X,'FIBER BENDING MODULUS (WARP DIRECTION). . .=',1PG20.13/
     & 5X,'NOMINAL STRETCH . (WEFT DIRECTION). . . . .=',1PG20.13/
     & 5X,'FIBER BENDING MODULUS (WEFT DIRECTION). . .=',1PG20.13/)
C
 1200 FORMAT(
     & 5X,'NOMINAL STRETCH . (WARP DIRECTION). . . . .=',1PG20.13/
     & 5X,'FIBER BENDING MODULUS (WARP DIRECTION). . .=',1PG20.13/
     & 5X,'NOMINAL STRETCH . (WEFT DIRECTION). . . . .=',1PG20.13/
     & 5X,'FIBER BENDING MODULUS (WEFT DIRECTION). . .=',1PG20.13/)

 1300 FORMAT(
     & 5X,'INTERSECTION ABSCISSA (WARP DIRECTION). . .=',1PG20.13/
     & 5X,'INTERSECTION ORDINATE (WARP DIRECTION). . .=',1PG20.13/
     & 5X,'INTERSECTION ABSCISSA (WEFT DIRECTION). . .=',1PG20.13/
     & 5X,'INTERSECTION ORDINATE (WEFT DIRECTION). . .=',1PG20.13/
     & 5X,'INTERSECTION ABSCISSA (SHEAR POINT1  ). . .=',1PG20.13/
     & 5X,'INTERSECTION ORDINATE (SHEAR POINT1  ). . .=',1PG20.13/
     & 5X,'INTERSECTION ABSCISSA (SHEAR POINT2  ). . .=',1PG20.13/
     & 5X,'INTERSECTION ORDINATE (SHEAR POINT2  ). . .=',1PG20.13/)
c-----------
      RETURN
      END
