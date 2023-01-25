c
c=======================================================================
Chd|====================================================================
Chd|  MINCF                         src/mater1/mincf.F            
Chd|-- called by -----------
Chd|        LAW98_UPD                     source/materials/mat/mat098/law98_upd.F
Chd|-- calls ---------------
Chd|        VF13AD                        src/mater1/mincf.F            
Chd|====================================================================
      SUBROUTINE MINCF(dvarl,dvaru,N,L,NDC,LEQ,X,f,df,c,dc,
     .                 EPS,MAXFUN,MAXLIN,W,IW,iiw,ICONTR,
     .                 igoto,ICONV,step_mincf)
c=======================================================================
C   MINIMIZATION OF A NON-LINEAR FUNCTION SUBJECT TO NON-LINEAR
C   CONSTRAINTS (MATHEMATICAL PROGRAMMING).
C   THIS SUBROUTINE CALLS SUBROUTINE VF13AD FROM THE HARWELL LIBRARY.
C
C   FOR A PROGRAM DESCRIPTION SEE:
C   KAJ MADSEN, OLE TINGLEFF, PER CHRISTIAN HANSEN & WOJCIECH OWCZARZ:
C   'ROBUST SUBROUTINES FOR NON-LINEAR OPTIMIZATION', REPORT NI 90-06,
C   INSTITUTE FOR NUMERICAL ANALYSIS, TECHNICAL UNIVERSITY OF DENMARK,
C   APTIL 1990
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,L,NDC,LEQ,MAXFUN,ICONTR,IGOTO,IW,MAXLIN,ICONV
      DOUBLE PRECISION F,EPS,step_mincf,GMAX1
      LOGICAL lopened
      INTEGER IIW(N+1)
      DOUBLE PRECISION X(N),W(IW),df(n),c(NDC),dc(n+1,NDC)
      DOUBLE PRECISION DVARL(N),DVARU(N)
c
      INTEGER icall, lwow, np1,lp, iprint, ilofcn, i, k
      DOUBLE PRECISION FSCAL
c
      COMMON /SQPSAV/ILOFCN
      COMMON /MINCFPRT/IPRINT, LP
      COMMON /SCAL_SQP/ FSCAL
c=======================================================================
c     N      nb of design variables N > 0
c     NDC    nb of variable constraints
c     C(NDC)      table of constraints
c     DC(N+1,NDC) constraint derivatives
c     F           Target function value
c     DF(N)       Function derivatives
c     DVARL(N)  lower bounds of design variables
c     DVARU(N)  upper bounds of design variables
c     L    Total number of constraints         L > 0
c     LEQ  Number of equality constraints      0 <= LEQ <= MIN(L,N)
c     X(N) Input table, depending on ICONTR
c        ICONTR > 0 :  on entry : X(N) = inital approximation to x
c                      on exit  : X(N) = computed solution <=> table of fited variables
c        ICONTR < 0 :  X(N) = point at which the gradient should be checked (input only)
c     DX      only used with ICONTR <=0 on entry
c     EPS     only used with ICONTR > 0, desired accuracy of the results, EPS > 0
c     MAXFUN  only used with ICONTR > 0
c                    on entry :  Upper bound on the number of calls of FDF (user objective function)  
c                    on exit  :  effective number of calls of FDF
c     MAXLIN  upper bound of calls of VF13AD in a line search
c     W(IW)   workspace table
c     IW      length of workspace table
c     IIW     ???
c     ICONTR  on entry :  control flag
c                 ICONTR = 1 : Start minimization.
c                 ICONTR = 2 : Start minimization and print information during the iteration.
c                 ICONTR < 0 : Check gradient. No iteration

c            On exit: Information about performance
c            ICONTR = 1 : Successful call.
c            ICONTR = 2 : Iteration stopped because too many iterations were needed, see MAXFUN.
c            ICONTR = 3 : Iteration stopped because more than 5 calls of FDFCDC 
c                         was needed in one line search. Check your gradients.
c            ICONTR = 4 : Iteration stopped because an uphill search direction was suggested. 
c                         Check your gradients.
c            ICONTR = 5 : Iteration failed because it was not possible to find 
c                         a starting point satisfying all constraints.
c            ICONTR < 0 : Computation did not start for the following reason,
c                   ICONTR = -2 : N <= 0
c                   ICONTR = -3 : L <= 0
c                   ICONTR = -4 : LEQ < 0 or LEQ > min(L,N)
c                   ICONTR = -6 : |DX| is too small
c                   ICONTR = -7 : EPS <= 0
c                   ICONTR = -8 : MAXFUN <= 0
c                   ICONTR = -9 : IW < 5/2 n(n+9) + (n+8)L + 15
c=======================================================================
c ... Set parameters
c----------------
      GMAX1  = 1.0D-3
      if (igoto == 8)  icall=0
      if (igoto == 5)  icall=1
      if (igoto == 3)  icall=2
      if (igoto /= 0)  goto 22
c
      icall  = 0
      ILOFCN = 0
      FSCAL  = 1.D0
      LWOW   = 5*L+12+(N*(5*N+57))/2+NDC
      IF (ICONTR < 0)       ICONTR=0
      IF (N <= 0)           ICONTR=-2
      IF (L < 0)            ICONTR=-3
      IF (LEQ < 0 .OR. LEQ > L) ICONTR=-4
      IF (EPS < 0D0)        ICONTR=-7
      IF (MAXFUN <= 0)      ICONTR=-8
      IF (IW < LWOW)        ICONTR=-10
C     ****                                         ***
C
      IF (ICONTR < 0) THEN
        IGOTO = 0
        RETURN
      ENDIF
C     ****    EXIT IF ERROR IN INPUT PARAMETERS    ***
c-----------------------------------------------------------------------
      if (IPRINT.GT.0) then
        lopened = .false.
        inquire(unit=LP, opened=lopened)
        if (.not. lopened) then
!         the file has not been opened yet, so no prints are allowed
          IPRINT = -1
        endif
      endif
c
21    continue
      if (igoto == 0) igoto = 8
      return

22    continue
      igoto = 0

      ILOFCN=ILOFCN+1
c
c        icall=0 is full function and gradient evaluation for all functions
c        icall=1 is function evaluation only for all functions
c        icall=2 is gradient evaluation only (by FD steps) for all functions 
c                   (functions at a point are assumed to be already known)
c
      NP1 = N+1
      CALL VF13AD(dvarl,dvaru,N,L,NDC,LEQ,X,f,df,c,dc,NP1,
     .           MAXFUN,MAXLIN,EPS,ICONTR,
     .           W,IW,IIW,icall,step_mincf,GMAX1)
c
      if (icall == 0) igoto = 8
      if (icall == 1) igoto = 5
      if (icall == 2) igoto = 3
c 
      IF (ICONTR == 0) GO TO 21
c
c     converged
c
      igoto = 0
      ICONV = 1
c
      MAXFUN=ILOFCN
C
      IF (FSCAL /= 1.D0) THEN
        k = L+2*N
        DO I=1,L
          W(k+I) = W(k+I)/FSCAL
        ENDDO
      ENDIF
c-----------
      RETURN
      END

C
Chd|====================================================================
Chd|  DIFERR                        src/mater1/mincf.F            
Chd|-- called by -----------
Chd|-- calls ---------------
Chd|====================================================================
      SUBROUTINE DIFERR(YM,YF,YB,DIFF,STEPL,ERM,ERF,ERB,CHM,CHF,CHB)
      IMPLICIT NONE
      DOUBLE PRECISION YM,YF,YB,DIFF,STEPL,ERM,ERF,ERB,S
      LOGICAL CHM,CHF,CHB
c-----------
        S=(YF-YM)/STEPL-DIFF
        CHF=(DABS(S) .GT. ERF)
        IF (CHF) ERF=DABS(S)
        S=(YM-YB)*2D0/STEPL-DIFF
        CHB=(DABS(S) .GT. ERB)
        IF (CHB) ERB=DABS(S)
        S=(YM+YF/3D0-4D0*YB/3D0)/STEPL-DIFF
        CHM=(DABS(S) .GT. ERM)
        IF (CHM) ERM=DABS(S)
c-----------
      RETURN
      END
C
Chd|====================================================================
Chd|  VF13BD                        src/mater1/mincf.F            
Chd|-- called by -----------
Chd|        VF13AD                        src/mater1/mincf.F            
Chd|-- calls ---------------
Chd|        TRYVMU                        src/mater1/mincf.F            
Chd|        VF13CD                        src/mater1/mincf.F            
Chd|        ZEROD                         src/mater1/mincf.F            
Chd|====================================================================
      SUBROUTINE VF13BD (N, M, NDC, MEQ, X, FI, G, C, CN,
     1 LCN, MAXFUN, MAXLIN, ACC, INF, W, IWDIM, IW, NP,
     2  VMU, dvarl, dvaru, VLAM, XOPT, COPT, B,
     3   DELTA, PREDPC, GRADCD, GAMMA, GLFIN,
     4    XBEG, BDELTA, GRAD, XBEST, icall,
     5  step_mincf,GMAX1)
C
      IMPLICIT NONE
      INTEGER N,M,NDC,MEQ,LCN,MAXFUN,INF,IWDIM,NP,ICALL
      DIMENSION X(N),G(N),C(NDC),CN(LCN,NDC),W(IWDIM)
      DOUBLE PRECISION X,FI,G,C,CN,ACC,W
      DOUBLE PRECISION dvarl(N),dvaru(N),XBEST(N)
      INTEGER IW(NP)
      DOUBLE PRECISION VLAM(M+2*N),VMU(M+2*N),PREDPC(M+2*N),
     1 XOPT(N), COPT(NDC), XBEG(N), GRADCD(M+2*N), DELTA(NP), GLFIN(N),
     2  B(NP,NP), GAMMA(N), BDELTA(N), GRAD(NP)
C
C     VLAM IS THE VECTOR OF LAGRANGE MULTIPLIERS
C     dvarl AND dvaru ARE PASSED TO VF13CD
C     VMU HOLDS THE PARAMETERS FOR THE LINE SEARCH FUNCTION
C     PREDPC GIVES THE PREDICTED REDUCTIONS IN CONSTRAINT VIOLATIONS
C     XOPT IS THE OPTIMAL CALCULATED VECTOR OF VARIABLES
C     COPT HOLDS THE CONSTRAINT FUNCTIONS FOR XOPT
C     XBEG GIVES THE VARIABLES AT THE START OF A LINE SEARCH
C     GRADCD GIVES DIRECTIONAL DERIVATIVES OF THE CONSTRAINTS
C     DELTA IS THE SEARCH DIRECTION OF AN ITERATION
C       PART OF DELTA IS USED AS WORKSPACE IN VF13CD
C     GAMMA IS THE SCALED CHANGE IN GRADIENT IN THE UPDATING OF B
C     GLFIN IS THE FINAL LAGRANGE FUNCTION GRADIENT OF AN ITERATION
C     B IS THE VARIABLE METRIC MATRIX
C     BDELTA IS B TIMES DELTA
C     GRAD IS THE GRADIENT OF THE QP OBJECTIVE FUNCTION
C     CRES GIVES THE CONSTRAINT RESIDUALS OF THE QP CALCULATION
C     W IS A WORKSPACE ARRAY OF LENGTH IWDIM FOR THE QP CALCULATION
C
C
C     VF13ED IS USED TO HOLD PRIVATE VARIABLES WHICH MUST HAVE
C       THEIR VALUES PRESERVED FOR SUBSEQUENT CALLS
C
      DOUBLE PRECISION  RATBD, RATNEW, SLOPE, STEPBD,
     1 DETB, TINY, FOPT, WOPT, WVAL, GRADFD, SUM, PREDMC,
     2  DBD, TARG, RATIO, AUX, TSWDOG, FBEG, VLBEG, WBEG,
     3   WPRED, VL, DG, THETA, THCOMP, step_mincf
C
      DOUBLE PRECISION CMAX(10), CL, CU, CMAXBEST, FBEST,
     *                   STEP0, STEPL, STEPU
      INTEGER ISCREEN,NMB,JVMU,LCMAX
C
      INTEGER MNN,MAXLIN,NWDOGT,NWDXXX,INFO,ITER,ITERP,IPRINT,ITERXX,
     *         NPR,MPR,ISWDOG,NFUN,I,J,K,NWDOGL,NFBEG,LP,IC
      DOUBLE PRECISION CT,CTEMP,TEMP,T1,T2,SMAL,GMAX1,FPRED,FBK,CBK,
     *         F,FSCAL,rScreen,rScreen0
C
      INTEGER IndexDesign0,INDEXD1,INDEXD2
C
      COMMON /MINCFPRT/IPRINT, LP
      COMMON /VF13ED/ FOPT, WOPT, WVAL, GRADFD, PREDMC, DBD,
     1 TARG, RATIO, TSWDOG, FBEG, FPRED, VLBEG, WPRED, WBEG,
     2  INFO, ITER, ITERP, ITERXX, NPR, MPR,
     3   ISWDOG, NFUN, NWDOGL, NFBEG
      COMMON /SCREEN_SQP/ CMAX,CL,CU,STEP0,STEPL,STEPU,FBK,CBK,
     *                CMAXBEST,FBEST,F,rScreen,
     *                ISCREEN,NMB,JVMU
      COMMON /SCREEN_RATE/ rScreen0
      COMMON /SCAL_SQP/ FSCAL
      COMMON /INDEX_SQP2/ INDEXD1,INDEXD2
      COMMON /INDEX_SQP/IndexDesign0 !PMO
C
C     SET THE FOLLOWING PARAMETERS
C     MAXLIN BOUNDS THE CALLS OF VF13AD IN A LINE SEARCH
C     NWDOGT IS THE T-PARAMETER OF THE WATCHDOG TECHNIQUE
C     NWDXXX RESTRICTS THE FREQUENCY OF UNSUCCESSFUL WATCHDOG ATTEMPTS
C     RATBD IS THE SMALLEST ACCEPTABLE VALUE AND RATNEW IS THE VALUE
C       AFTER REVISION OF THE LINE SEARCH GRADIENT RATIO IN THE
C       ADJUSTMENT OF VMU
C     SLOPE IS THE ARMIJO PARAMETER OF THE LINE SEARCH (ITS VALUE
C       MUST BE LESS THAN 0.5)
C     STEPBD IS A LOWER BOUND ON THE STEP LENGTH REDUCTION
C     DETB IS A LOWER BOUND ON THE DETERMINANT REDUCTION
C     TINY IS THE NORMAL INITIAL VALUE OF EACH COMPONENT OF VMU
C
      MNN=M+N+N
      NWDOGT=2
      NWDXXX=10
      LCMAX=8
      RATBD=1.5D0
      RATNEW=2.0D0
      SLOPE=0.1D0
      STEPBD=0.1D0
      DETB=0.1D0
      TINY=1.0D-20
      SMAL=1.0D-10
C
      IF ((ICALL == 0).OR.(FSCAL == 1.0D0)) THEN
        F = FI
      ELSE
        IF (ICALL == 2) THEN
          DO I=1,N
            G(I) = FSCAL*G(I)
          ENDDO
        ELSE IF (ICALL == 1) THEN
          F = FSCAL*FI
        ENDIF
      ENDIF

c  if coming back with derivatives (icall=2) go straight to evaluation of new search direction
       if (icall == 2) go to 351
c
c
C     INITIALISE ONLY IF THE FIRST ENTRY
c
c       write (*,*)'when inf =0 , skip setting init. conds, goto 305' 
C
      IF (INF == 0) GO TO 305
c
c       write(*,*)'when inf ne 0 , set init. conds, inf will be 0'
C
C     SET THE INITIAL CONDITIONS
C     ITER IS THE ITERATION NUMBER
C     ITERP IS THE NEXT ITERATION ON WHICH PRINTING IS GIVEN
C     ITERXX IS USUALLY THE ITERATION OF THE MOST RECENT BACK-TRACKING
C     NPR AND MPR ARE USED IN THE CONTROL OF FORMAT STATEMENTS
C     ISWDOG=1 GIVES THE STANDARD LINE SEARCH CRITERION
C     NFUN IS THE NUMBER OF CALLS OF VF13AD
C
      FOPT    = 0.D0
      WOPT    = 0.D0
      WVAL    = 0.D0
      GRADFD  = 0.D0
      PREDMC  = 0.D0
      DBD     = 0.D0
      TARG    = 0.D0
      RATIO   = 0.D0
      TSWDOG  = 0.D0
      FBEG    = 0.D0
      FPRED   = 0.D0
      VLBEG   = 0.D0
      WPRED   = 0.D0
      WBEG    = 0.D0
      CALL ZEROD(GLFIN,N)
      CALL ZEROD(GAMMA,N)
      CALL ZEROD(CMAX,LCMAX)
      CL      = 0.D0
      CU      = 0.D0
      STEP0   = 0.D0
      STEPL   = 0.D0
      STEPU   = 0.D0
      FBK     = 0.D0
      CBK     = 0.D0
      NWDOGL  = 0
      NFBEG   = 0
      INFO=INF
      INF=0
      ITER=0
      ITERP=MIN0(IPRINT-1,0)
      ITERXX=-NWDXXX
      NPR=MIN0(N,5)
      MPR=MIN0(M,5)
      ISWDOG=1
      NFUN=1
      ISCREEN=0
      NMB=0
      JVMU=0
      CMAXBEST=1.0D20
      FBEST=1.0D20
      FSCAL=1.0D0
      IF (rScreen0 == 0.8D0) THEN
        rScreen=0.8D0
      ELSE
        rScreen=0.95D0
      ENDIF
      IF (GMAX1.LT.0.0D0 .OR. GMAX1.GT.1.0D0) GMAX1 = 5.0D-3
      INDEXD1 = 1
      INDEXD2 = 1
C
C     SET THE INITIAL ELEMENTS OF VMU AND B
C
      IF (INFO  ==  -111) GO TO 20
      IF (INFO  ==  -110) GO TO 20
      DO 10 K=1,MNN
   10 VMU(K)=TINY
      IF (INFO .NE. -101) INFO=-2
   20 DO 40 I=1,N
      DO 30 J=1,N
   30 B(I,J)=0.0D0
   40 B(I,I)=1.0D0
C
C     INITIALISATION COMPLETE
C
      GO TO 305
C
C     STORE THE OPTIMAL X AND C
C     NWDOGL IS THE ITERATION NUMBER THAT GIVES THE OPTIMAL X
C     FOPT AND WOPT ARE THE OPTIMAL VALUES SO FAR OF THE OBJECTIVE
C       FUNCTION AND THE LINE SEARCH OBJECTIVE FUNCTION
C     WVAL IS THE CURRENT VALUE OF THE LINE SEARCH OBJECTIVE FUNCTION
C
   50 NWDOGL=ITER
   60 DO 70 I=1,N
   70 XOPT(I)=X(I)
      FOPT=F
      WOPT=WVAL
      INDEXD1 = IndexDesign0
      DO 80 K=1,M
   80 COPT(K)=C(K)
      IF (NWDOGL .LT. ITER) GO TO 220
C
      CMAX(1) = 0.0D0
      DO I=1,M
        IF (I.LE.MEQ .OR. C(I).LT.0.0D0) THEN
          IF (CMAX(1).LT.ABS(C(I))) CMAX(1) = ABS(C(I))
        ENDIF
      ENDDO
C
85    CONTINUE
      IF (ISCREEN.GT.0) THEN
C       do constraint screening
        DO I=1,M
          IF (I.LE.MEQ) THEN
            CT = ABS(COPT(I))
          ELSE
            CT = -COPT(I)
          ENDIF
          IF (CT.LT.(CMAX(1)*rScreen)) THEN
            C(I) = 0.0D0
            DO J=1,N
              CN(J,I) = 0.0D0
            ENDDO
          ELSE
            C(I) = COPT(I)
          ENDIF
          CN(NP,I) = 0.0D0
        ENDDO
C
        DO I=1,N
          G(I) = 0.0D0
        ENDDO
      ENDIF
C
C     CALL THE QUADRATIC PROGRAMMING SUBROUTINE AND TEST WHETHER THE
C       VALUE OF INFO INDICATES AN ERROR RETURN
C
   90 CALL VF13CD (N, M, NDC, MEQ, G, C, CN, LCN,
     1  NP, VLAM, B, DELTA, dvarl, dvaru, X,
     2  GRAD, W, IWDIM, IW,INFO)
      IF (INFO  ==  111) GO TO 523
      IF (INFO  ==  5) GO TO 510
C
C     TEST FOR CONVERGENCE
C     GRADFD IS THE DIRECTIONAL DERIVATIVE OF F
C     INFO=1 INDICATES THAT THE CURRENT X SATISFIES THE TERMINATION TEST
C
      GRADFD=0.0D0
      DO 100 I=1,N
  100 GRADFD=GRADFD+DELTA(I)*G(I)
      SUM=DABS(GRADFD)
      DO 110 K=1,M
  110 SUM=SUM+DABS(VLAM(K)*C(K))
      DO I=1,N
        SUM=SUM+DABS(VLAM(M+2*I-1)*(X(I)-dvarl(I)))
        SUM=SUM+DABS(VLAM(M+2*I)*(dvaru(I)-X(I)))
      ENDDO
      IF (ISCREEN.NE.0) GO TO 120
      IF (SUM .GT. ACC) GO TO 120
C
C     dispose violated constraints
C
      TEMP = 0.0D0
      DO I=1,M
        K = 0
        IF (I.LE.MEQ) THEN
          IF (DABS(C(I)).GT.GMAX1) K = 1
        ELSE
          IF (C(I).LT.-GMAX1) K = 1
        ENDIF
        IF (K == 1) THEN
          T1 = -2.0D0*ACC/C(I)
          IF (VLAM(I) == 0.0D0) THEN
            VLAM(I) = T1
          ELSE IF (DABS(VLAM(I)).LT.DABS(T1)) THEN
            VLAM(I) = SIGN(T1,VLAM(I))
          ENDIF
          IF (TEMP.LT.DABS(C(I))) TEMP=DABS(C(I))
        ENDIF
      ENDDO
      IF (TEMP.GT.GMAX1) GOTO 120
C
C     SQP is terminated at the initial design because of small gradient
C     scale objective gradient to have a try
      IF ((ICALL == 0).AND.(FSCAL == 1.D0)) THEN
        TEMP = 0.0D0
        DO I=1,N
          TEMP = TEMP+G(I)*G(I)
        ENDDO
        TEMP = DSQRT(TEMP)
        IF ((TEMP.LT.0.99D0).AND.(TEMP.GE.1.0D-10)) THEN
          FSCAL = 1.0D0/TEMP
          F = FSCAL*F
          DO I=1,N
            G(I) = FSCAL*G(I)
            X(I) = XOPT(I)
          ENDDO
          IndexDesign0 = INDEXD1
          INFO=-2
          GOTO 305
        ENDIF
      ENDIF
C
      INFO=1
      IF (WVAL .LE. WOPT) GO TO 550
C
C     SET THE ELEMENTS OF GRADCD AND PREDPC
C     SET PREDMC TO THE PREDICTED REDUCTION IN THE PENALTY TERM OF
C       THE LINE SEARCH OBJECTIVE FUNCTION
C
  120 PREDMC=0.0D0
      DO 140 K=1,MNN
      GRADCD(K)=0.0D0
      IF(K.GT.M) GO TO 131
      DO 130 I=1,N
  130 GRADCD(K)=GRADCD(K)+DELTA(I)*CN(I,K)
      PREDPC(K)=DMAX1(0.0D0,-C(K)+DMIN1(0.0D0,C(K)+GRADCD(K)))
      IF (K .GT. MEQ) GO TO 140
      PREDPC(K)=DMAX1(0.0D0,DABS(C(K))-DABS(C(K)+GRADCD(K)))
      GO TO 140
  131 I=MOD(K-M,2)
      IF(I == 0) GO TO 132
      I=INT((K-M)/2)+1
      GRADCD(K)=DELTA(I)
      CTEMP=X(I)-dvarl(I)
      GO TO 135
  132 I=INT((K-M)/2)
      GRADCD(K)=-DELTA(I)
      CTEMP=-X(I)+dvaru(I)
  135 PREDPC(K)=DMAX1(0.0D0,-CTEMP+DMIN1(0.0D0,CTEMP+GRADCD(K)))
  140 PREDMC=PREDMC+VMU(K)*PREDPC(K)
C
C     BRANCH IF THERE IS NO NEED TO INCREASE VMU
C     THE VECTOR BDELTA IS CALCULATED FOR THIS TEST, AND DBD IS SET TO
C       THE SCALAR PRODUCT (DELTA,BDELTA).  BOTH ARE USED AGAIN LATER
C
      DBD=0.0D0
      DO 160 I=1,N
      BDELTA(I)=0.0D0
      DO 150 J=1,N
  150 BDELTA(I)=BDELTA(I)+B(I,J)*DELTA(J)
  160 DBD=DBD+DELTA(I)*BDELTA(I)
C
C     TEST FOR LOSS OF POSITIVE DEFINITENESS.
C     THIS CAN ONLY OCCUR BECAUSE OF COMPUTER
C       ROUNDING ERRORS
C
      IF (DBD .LE. 0.0D0) GO TO 523
      IF ((GRADFD+0.5D0*DBD) .LE. 0.0D0) GO TO 230
      IF (PREDMC .GE. (RATBD*DABS(GRADFD)+SMAL)) GO TO 230
C
C     REVISE THE COMPONENTS OF VMU AND REVISE WOPT AND WVAL
C     THE INCREASE IN VMU IS SUCH THAT TARG IS THE PREDICTED REDUCTION
C       IN THE PENALTY TERM OF THE NEW LINE SEARCH FUNCTION
C
      JVMU=1
      TARG=RATNEW*DABS(GRADFD)+SMAL
      SUM=SMAL
      DO 170 K=1,MNN
      IF (VLAM(K) .NE. 0.0D0) GO TO 170
      TARG=TARG-VMU(K)*PREDPC(K)
      PREDPC(K)=0.0D0
  170 SUM=SUM+DABS(VLAM(K)*PREDPC(K))
  180 RATIO=TARG/SUM
      AUX=-1.0D0
      DO 190 K=1,MNN
      IF (PREDPC(K)  ==  0.0D0) GO TO 190
      IF (VMU(K) .LE. RATIO*DABS(VLAM(K))) GO TO 190
      AUX=1.0D0
      TARG=TARG-VMU(K)*PREDPC(K)
      SUM=SUM-DABS(VLAM(K)*PREDPC(K))
      PREDPC(K)=0.0D0
  190 CONTINUE
      IF (AUX .GE. 0.0D0) GO TO 180
      WOPT=FOPT
      WVAL=F
      DO 210 K=1,M
      IF (K.LE.MEQ) THEN
        IF (PREDPC(K) .NE. 0.0D0) VMU(K)=RATIO*DABS(VLAM(K))
      ELSE
        IF (PREDPC(K) .GT. 0.0D0) VMU(K)=RATIO*DABS(VLAM(K))
      ENDIF
      IF (K .GT. MEQ) GO TO 200
      WOPT=WOPT+VMU(K)*DABS(COPT(K))
      WVAL=WVAL+VMU(K)*DABS(C(K))
      GO TO 210
  200 WOPT=WOPT+VMU(K)*DMAX1(0.0D0,-COPT(K))
      WVAL=WVAL+VMU(K)*DMAX1(0.0D0,-C(K))
  210 CONTINUE
      DO I=1,N
        K=M+2*I-1
        IF (PREDPC(K) .GT. 0.0D0) VMU(K)=RATIO*DABS(VLAM(K))
        WOPT=WOPT+VMU(K)*DMAX1(0.0D0,-XOPT(I)+dvarl(I))
        WVAL=WVAL+VMU(K)*DMAX1(0.0D0,-X(I)+dvarl(I))
        IF (PREDPC(K+1) .GT. 0.0D0) VMU(K+1)=RATIO*DABS(VLAM(K+1))
        WOPT=WOPT+VMU(K+1)*DMAX1(0.0D0,XOPT(I)-dvaru(I))
        WVAL=WVAL+VMU(K+1)*DMAX1(0.0D0,X(I)-dvaru(I))
      ENDDO
C
C     SET TSWDOG, WHICH IS USED TO DETERMINE THE LINE SEARCH CRITERION
C       OF THE NEXT ITERATION, TO THE NEW VALUE OF WOPT
C     REVISE THE OPTIMAL X AND C IF NECESSARY
C
      IF (NWDOGL  ==  ITER) GO TO 230
      TSWDOG=WOPT
      IF (WOPT .LT. WVAL) GO TO 230
      IF (INFO  ==  1) GO TO 550
      GO TO 60
  220 NWDOGL=ITER
C
C     PROVIDE PRINTING IF REQUIRED
C
  230 INFO=0
      IndexDesign0 = 0
      IF (ITER  ==  ITERP) GO TO 570
C
C     NFBEG IS THE NUMBER OF FUNCTION CALLS AT THE START OF THE LINE
C       SEARCH
C     STORE THE ELEMENTS OF XBEG
C     STORE THE CONSTRAINT FUNCTIONS FOR XBEG IN THE VECTOR PREDPC
C     STORE THE INITIAL LAGRANGE FUNCTION GRADIENT OF AN ITERATION
C       IN THE VECTOR GAMMA
C     FBEG IS THE VALUE OF F AT THE START OF THE LINE SEARCH
C     VLBEG IS THE VALUE OF THE LAGRANGIAN APPROXIMATION AT THE START
C       OF THE LINE SEARCH
C     STEP IS THE STEP LENGTH OF THE LINE SEARCH
C
  240 NFBEG=NFUN
      GO TO 350
  250 FBEG=F
      DO 260 I=1,N
      GAMMA(I)=GLFIN(I)
  260 XBEG(I)=X(I)
C
C     BACKUP THE BEST DESIGN
      IF (M.GT.0) THEN
        CTEMP = 0.0D0
        DO I=1,M
          IF (I.LE.MEQ .OR. C(I).LT.0.0D0) THEN
            IF (CTEMP.LT.ABS(C(I))) CTEMP = ABS(C(I))
          ENDIF
        ENDDO
        K = 0
        IF (CMAXBEST.GT.GMAX1) THEN
          IF (CTEMP.LT.CMAXBEST) K=1
        ELSE IF (CTEMP.LE.GMAX1) THEN
          IF (F.LT.FBEST) K=1
        ENDIF
        IF (K == 1) THEN
          FBEST=F
          CMAXBEST=CTEMP
          INDEXD2 = IndexDesign0
          DO I=1,N
            XBEST(I)=X(I)
          ENDDO
        ENDIF
      ENDIF
C
      VLBEG=F
      DO 270 K=1,M
      PREDPC(K)=C(K)
  270 VLBEG=VLBEG-VLAM(K)*C(K)
      DO I=1,N
        K=M+2*I-1
        PREDPC(K)=X(I)-dvarl(I)
        VLBEG=VLBEG-VLAM(K)*(X(I)-dvarl(I))
        PREDPC(K+1)=-X(I)+dvaru(I)
        VLBEG=VLBEG-VLAM(K+1)*(-X(I)+dvaru(I))
      ENDDO
      WBEG=WVAL
!      step_mincf = 1.0D0
C
C     CALCULATE THE NEW X UNLESS NFUN IS TOO LARGE
C
  280 IF (NFUN .GE. MAXFUN) GO TO 450
      DO 290 I=1,N
  290 X(I)=XBEG(I)+step_mincf*DELTA(I)
  300 NFUN=NFUN+1
C
C     RETURN FOR MORE VALUES OF F, G, C, CN
c
      ICALL = 1
c
c       write(*,*) 'SQP1 Line search, NFUN, STEP', nfun, step
c      write(*,*) 'SQP2 x out', x(1),x(2)
C
      RETURN
  305 CONTINUE
c
c      write(*,*) 'SQP3 F in', F
c       write(*,*) 'SQP4 C in', c(1),c(2),c(3),c(4)
c
      IF (INFO .GT. 0)  GO TO 550
      IF (INFO  ==  -1) GO TO 440
      IF (INFO  ==  -3) THEN
C       NOW SWITCHED TO THE BEST DESIGN, TRY TO GET GRADIENTS
        INF=0
        CALL TRYVMU(3,FOPT,FBEG,FPRED,F,WOPT,WBEG,WPRED,WVAL,
     &              COPT,C,VMU,M,MEQ)
        FBEG=F
        DO I=1,N
          XBEG(I)=X(I)
          GAMMA(I)=G(I)
        ENDDO
        DO I=1,NP
          DO J=1,NP
            B(I,J)=0.0D0
          ENDDO
          B(I,I)=1.0D0
        ENDDO
        B(NP,NP)=0.0D0
        call zerod(VLAM,M+2*N)
        call zerod(W,IWDIM)
        call zeroi(IW,NP)
        GOTO 350
      ENDIF
C
      IF (ISCREEN.GT.0) THEN
C       dispose screening branch
        IC = MIN0(NFUN-NFBEG+1,LCMAX)
        CMAX(IC) = 0.0D0
        DO I=1,M
          IF (I.LE.MEQ .OR. C(I).LT.0.0D0) THEN
            IF (CMAX(IC).LT.ABS(C(I))) CMAX(IC) = ABS(C(I))
          ENDIF
        ENDDO

        TEMP = (CMAX(1)-CMAX(IC))/CMAX(1)
        IF (TEMP.GT.1.0D-2) THEN
C         better solution was found
          IF (CMAX(IC).LT.GMAX1 .OR. TEMP.LT.1.0D-1) THEN
            ISCREEN = 0
            CALL TRYVMU(3,FOPT,FBEG,FPRED,F,WOPT,WBEG,WPRED,WVAL,
     &                  COPT,C,VMU,M,MEQ)
            WOPT=WVAL
          ENDIF
          IndexDesign0 = 0
          GOTO 350
        ELSE
          IF (IC.GE.LCMAX) THEN
C           cannot find better solution, return
            CALL TRYVMU(3,FOPT,FBEG,FPRED,F,WOPT,WBEG,WPRED,WVAL,
     &                  COPT,C,VMU,M,MEQ)
            CTEMP = 0.0D0
            DO I=1,M
              IF (I.LE.MEQ .OR. COPT(I).LT.0.0D0) THEN
                IF (CTEMP.LT.ABS(COPT(I))) CTEMP = ABS(COPT(I))
              ENDIF
            ENDDO
            IF (CMAX(LCMAX).LT.CTEMP) THEN
              WVAL=WOPT-1.0D0 ! ACCEPT SMALLER VIOLATION
            ELSE
              WVAL=WOPT+1.0D0 ! REVERSE BACK TO THE OPTIMA
            ENDIF
            GOTO 470
          ELSE
            IF (IC == 2) THEN
              STEPL = 0.0D0
              STEPU = 1.0D0
              CL    = CMAX(1)
              CU    = CMAX(2)
              IF (CU.GE.CL) THEN
                STEP0 = STEPL
              ELSE
                STEP0 = STEPU
              ENDIF
              step_mincf = 0.1D0
            ELSE
              IF (CMAX(IC).GE.CMAX(1) .AND. CMAX(IC-1).GE.CMAX(1)) THEN
                STEPU = step_mincf
                CU    = CMAX(IC)
                step_mincf = 0.3D0*step_mincf
              ELSE
                IF (STEP0 == 1.0D0 .AND. CMAX(IC).GE.CMAX(IC-1)) THEN
                  STEPL = step_mincf
                  CL    = CMAX(IC)
                  step_mincf = STEPL+0.7D0*(STEPU-STEPL)
                  CMAX(IC) = CMAX(IC-1)
                ELSE
                  IF (STEP0.LT.1.0D-10) THEN
                    STEP0 = step_mincf
                  ELSE
                    IF (CMAX(IC).LE.CMAX(IC-1)) THEN
                      IF (STEP0.GT.step_mincf) THEN
                        STEPU = STEP0
                        CU    = CMAX(IC-1)
                      ELSE
                        STEPL = STEP0
                        CL    = CMAX(IC-1)
                      ENDIF
                      STEP0 = step_mincf
                    ELSE
                      IF (STEP0.GT.step_mincf) THEN
                        STEPL = step_mincf
                        CL    = CMAX(IC)
                      ELSE
                        STEPU = step_mincf
                        CU    = CMAX(IC)
                      ENDIF
                      CMAX(IC) = CMAX(IC-1)
                    ENDIF
                  ENDIF
                  T1 = (CL-CMAX(IC))*(STEP0-STEPU)
                  T2 = (CMAX(IC)-CU)*(STEPL-STEP0)
                  IF (ABS(T1-T2).GT.1.0D-10) THEN
                    step_mincf = T1*(STEP0+STEPU)-T2*(STEPL+STEP0)
                    step_mincf = step_mincf/(T1-T2)/2.0D0
                  ELSE
                    IF (ABS(CL+CU).GT.1.0D-6) THEN
                      TEMP  = CL/(CL+CU)
                      step_mincf = STEPL+(STEPU-STEPL)*TEMP
                    ELSE
                      step_mincf = (STEPL+STEPU)/2.0D0
                    ENDIF
                  ENDIF
                  T1 = STEPL+(STEPU-STEPL)*0.1D0
                  T2 = STEPL+(STEPU-STEPL)*0.9D0
                  IF (step_mincf.LT.T1) THEN
                    step_mincf = T1
                  ELSE IF (step_mincf.GT.T2) THEN
                    step_mincf = T2
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            GOTO 280
          ENDIF
        ENDIF
      ENDIF
c
      WVAL=F
      DO 310 K=1,M
      AUX=DMAX1(0.0D0,-C(K))
      IF (K .LE. MEQ) AUX=DABS(C(K))
  310 WVAL=WVAL+AUX*VMU(K)
      DO I=1,N
        K=M+2*I-1
        AUX=DMAX1(0.0D0,-X(I)+dvarl(I))
        WVAL=WVAL+AUX*VMU(K)
        AUX=DMAX1(0.0D0,X(I)-dvaru(I))
        WVAL=WVAL+AUX*VMU(K+1)
      ENDDO
      IF (INFO .LE. -2) GO TO 50
C
C     SET WPRED TO THE PREDICTED VALUE OF WVAL THAT IS OBTAINED FROM
C       WBEG AND THE CALCULATED FIRST DERIVATIVES
C
      FPRED=FBEG+step_mincf*GRADFD
      WPRED=FPRED
      DO 320 K=1,MNN
      AUX=PREDPC(K)+step_mincf*GRADCD(K)
      IF (K .LE. MEQ) AUX=DABS(AUX)
      IF (K .GT. MEQ) AUX=DMAX1(0.0D0,-AUX)
  320 WPRED=WPRED+AUX*VMU(K)
C
C     TEST FOR AN ERROR RETURN AND FOR THE END OF THE LINE SEARCH
C     VL IS THE NEW VALUE OF THE LAGRANGIAN FUNCTION
C
      IndexDesign0 = 0
      IF (WPRED .GE. WBEG) GO TO 490
      AUX=(WBEG-WVAL)/(WBEG-WPRED)
      IF (AUX .GE. SLOPE) GO TO 350
      IF (ISWDOG  ==  1) GO TO 340
      IF (WVAL .LE. WBEG) GO TO 350
      VL=F
      DO 330 K=1,M
  330 VL=VL-VLAM(K)*C(K)
      DO I=1,N
        K=M+2*I-1
        VL=VL-VLAM(K)*(X(I)-dvarl(I))
        VL=VL-VLAM(K+1)*(-X(I)+dvaru(I))
      ENDDO
      IF (VL .LE. VLBEG) GO TO 350
C
C     REVISE THE STEP LENGTH TO CONTINUE THE LINE SEARCH
C
  340 IF (NFUN .GE. NFBEG+MAXLIN) GO TO 470
      step_mincf=DMAX1(STEPBD,0.5D0/(1.0D0-AUX))*step_mincf
      GO TO 280
C
C     CALCULATE THE GRADIENT OF THE LAGRANGIAN FUNCTION
C
  350   continue
c  We are here when the line search was successfully completed.
c  Get the derivatives only:
c
c XXX
c      if (icall.ne.0) then
       if ((icall.ne.0).and.(icall.ne.2)) then
c go out for derivatives
          icall=2
          return
       endif
c
  351  continue
c      came back with derivatives
c
      IF (IPRINT.GT.1) THEN
        WRITE(LP,'(A30)') 'COME BACK WITH DERIVATIVES'
        WRITE(LP,'(A9,I5)') 'INFO = ',INFO
        WRITE(LP,'(A9)') 'VLAM ='
        WRITE(LP,'(10E13.5)') (VLAM(I),I=1,MNN)
        WRITE(LP,'(A9)') 'G ='
        WRITE(LP,'(10E13.5)') (G(I),I=1,N)
        WRITE(LP,'(A9)') 'CN ='
        DO J=1,M
          WRITE(LP,'(10E13.5)') (CN(I,J),I=1,N)
        ENDDO
        WRITE(LP,'(A9,I5)') 'NFUN = ',NFUN
        WRITE(LP,'(A9,I5)') 'NFBEG = ',NFBEG
        WRITE(LP,'(A9)') 'DELTA ='
        WRITE(LP,'(10E13.5)') (DELTA(I),I=1,N)
        WRITE(LP,'(A9,E13.5)') 'STEP = ',step_mincf
      ENDIF
c
      IF (INFO == -3) THEN
        INFO=0
        GOTO 50
      ENDIF
c
      DO I=1,N
        GLFIN(I)=G(I)
      ENDDO
      DO K=1,M
        IF (VLAM(K).NE.0.D0) THEN
          DO I=1,N
            GLFIN(I)=GLFIN(I)-VLAM(K)*CN(I,K)
          ENDDO
        ENDIF
      ENDDO
      DO K=M+1,MNN
        IF (VLAM(K).NE.0.D0) THEN
          J=K-M
          I=MOD(J,2)
          IF (I == 0) THEN
            I=INT(J/2)
            GLFIN(I)=GLFIN(I)+VLAM(K)
          ELSE
            I=INT(J/2)+1
            GLFIN(I)=GLFIN(I)-VLAM(K)
          ENDIF
        ENDIF
      ENDDO
      IF (IPRINT.GT.1) THEN
        WRITE(LP,'(A9)') 'GLFIN ='
        WRITE(LP,'(10E13.5)') (GLFIN(I),I=1,N)
      ENDIF
      IF (NFUN  ==  NFBEG) GO TO 250
C
C     CALCULATE GAMMA IN ORDER TO REVISE B
C     NOTE THAT BDELTA AND DBD ARE ALREADY AVAILABLE
C     SET DG TO THE SCALAR PRODUCT OF DELTA AND GAMMA
C
      DG=0.0D0
      DO 390 I=1,N
      GAMMA(I)=(GLFIN(I)-GAMMA(I))/step_mincf
  390 DG=DG+DELTA(I)*GAMMA(I)
C
C     OVERWRITE GAMMA AND DG IF NECESSARY IN ORDER TO SATISFY THE
C       DETERMINANT CONDITION
C
      AUX=DETB*DBD
      IF (IPRINT.GT.1) THEN
        WRITE(LP,'(A20)') 'IF (DG .GE. AUX)'
        WRITE(LP,'(A9,E13.5)') 'DG = ',DG
        WRITE(LP,'(A9,E13.5)') 'AUX = ',AUX
        WRITE(LP,'(A9)') 'GAMMA ='
        WRITE(LP,'(10E13.5)') (GAMMA(I),I=1,N)
      ENDIF
      IF (DG .GE. AUX) GO TO 410
      THETA=(DBD-AUX)/(DBD-DG)
      THCOMP=1.0D0-THETA
      DO 400 I=1,N
  400 GAMMA(I)=THETA*GAMMA(I)+THCOMP*BDELTA(I)
      DG=AUX
C
C     REVISE THE MATRIX B
C
  410 DO 420 I=1,N
      THETA=BDELTA(I)/DBD
      THCOMP=GAMMA(I)/DG
      DO 420 J=I,N
      B(I,J)=B(I,J)-THETA*BDELTA(J)+THCOMP*GAMMA(J)
  420 B(J,I)=B(I,J)
C
C     APPLY THE WATCHDOG TECHNIQUE
C     THE VALUE ISWDOG=0 GIVES THE RELAXED LINE SEARCH CRITERION
C     THE VALUE ISWDOG=1 GIVES THE STANDARD LINE SEARCH CRITERION
C
      IF (NWDOGL  ==  ITER) TSWDOG=WBEG-SLOPE*(WBEG-WPRED)
      ISWDOG=1
      IF (WVAL .LE. TSWDOG) ISWDOG=0
      IF (ITER .LT. ITERXX+NWDXXX) ISWDOG=1
      ITER=ITER+1
      IF (ISCREEN.NE.0) THEN
        NWDOGL = ITER
        GOTO 50
      ENDIF
      IF (WVAL .LE. WOPT) GO TO 50
      IF (ITER .LE. NWDOGL+NWDOGT) GO TO 90
C
C     THE NEXT INSTRUCTION IS REACHED IF IT IS NECESSARY TO BACK-TRACK.
C       IN THIS CASE X IS REPLACED BY XOPT
C
      ITERXX=ITER
      INFO=-1
      DO 430 I=1,N
  430 X(I)=XOPT(I)
      WVAL=WOPT
      IndexDesign0 = INDEXD1
      GO TO 300
  440 INFO=0
      NWDOGL=ITER
      GO TO 90
C
C     ERROR RETURNS FROM VF13AD
C
  450 INFO=2
      IF (IPRINT.LE.0) GO TO 530
      WRITE(LP,460) NFUN
  460 FORMAT (/5X,'ERROR RETURN FROM VF13AD DUE TO',I5,
     1' CALLS OF VF13AD')
      GO TO 530
  470 INFO=3
      IF (IPRINT.LE.0) GO TO 530
      WRITE(LP,480) MAXLIN
  480 FORMAT (/5X,'WARNING RETURN FROM VF13AD DUE TO',I5,
     1' CALLS OF VF13AD IN A LINE SEARCH')
      GO TO 530
  490 INFO=4
      IF (IPRINT.LE.0) GO TO 530
      WRITE(LP,500)
  500 FORMAT (/5X,'ERROR RETURN FROM VF13AD',
     1 ' DUE TO AN UPHILL SEARCH DIRECTION')
      GO TO 530
  510 IF (IPRINT.LE.0) GO TO 530
      WRITE(LP,520)
  520 FORMAT (/5X,'ERROR RETURN FROM VF13AD',
     1 ' BECAUSE THE CONSTRAINTS SEEM TO BE INCONSISTENT')
      GO TO 530
  523 INFO=6
      IF (IPRINT.LE.0) GO TO 530
      WRITE(LP,526)
  526 FORMAT(/5X,'ERROR RETURN FROM VF13AD DUE TO LOSS OF',/5X,
     1 'POSITIVE DEFINITENESS.  THE SOLUTION MAY BE UNBOUNDED')
  530 IF (WVAL .LE. WOPT) GO TO 550
      DO 540 I=1,N
  540 X(I)=XOPT(I)
      IndexDesign0 = INDEXD1
      GO TO 300
C
C     PRINT THE VALUES OF F, X AND C THAT ARE REQUESTED
C
  550 IF (IPRINT .LE. 0) GO TO 670
      WRITE(LP,560)
  560 FORMAT (/5X,'THE FOLLOWING VALUES OF F, X AND C',
     1 ' ARE RETURNED BY VF13AD')
  570 IF (IPRINT .LE. 0) GO TO 670
      WRITE(LP,580) ITER, NFUN
  580 FORMAT (/5X,'ITERATIONS =',I5,5X,'CALLS OF VF13AD =',I5)
      IF (INFO .GE. 1) GO TO 610
      WRITE(LP,590)
  590 FORMAT (/5X,'THE CURRENT VALUES OF F, X AND C ARE AS FOLLOWS')
      IF (ITER  ==  NWDOGL) WRITE(LP,600)
  600 FORMAT (5X,'THESE VALUES ARE ALSO THE OPTIMAL ONES')
  610 IF (IPRINT .LE. 0) GO TO 670
      WRITE(LP,620) (X(I),I=1,NPR)
  620 FORMAT (/2X,'X =',5D20.10)
      IF (N .GE. 6) WRITE(LP,630) (X(I),I=6,N)
  630 FORMAT (5X,5D20.10)
      WRITE(LP,640) F
  640 FORMAT (2X,'F =',D20.10)
      WRITE(LP,650) (C(K),K=1,MPR)
  650 FORMAT (2X,'C =',5D20.10)
      IF (M .GE. 6) WRITE(LP,630) (C(K),K=6,M)
      WRITE(LP,651) (VMU(K),K=1,MPR)
  651 FORMAT (2X,'VMU =',5D20.10)
      IF (M .GE. 6) WRITE(LP,630) (VMU(K),K=6,M)
      WRITE(LP,652) (DELTA(K),K=1,NPR)
  652 FORMAT (2X,'DELTA =',5D20.10)
      IF (N .GE. 6) WRITE(LP,630) (DELTA(K),K=6,N)
      WRITE(LP,'(A10,I3)') 'ISCREEN = ', ISCREEN
      WRITE(LP,'(A10,I3)') 'NMB     = ', NMB
      IF (INFO .GE. 1) GO TO 670
      ITERP=ITER+IPRINT
      IF (IPRINT .LE. 0) GO TO 240
      WRITE(LP,660)
  660 FORMAT (/5X,'AND THE OPTIMAL VALUES ARE THE FOLLOWING ONES')
      WRITE(LP,620) (XOPT(I),I=1,NPR)
      IF (N .GE. 6) WRITE(LP,630) (XOPT(I),I=6,N)
      WRITE(LP,640) FOPT
      WRITE(LP,650) (COPT(K),K=1,MPR)
      IF (M .GE. 6) WRITE(LP,630) (COPT(K),K=6,M)
      GO TO 240
  670 INF=INFO
C
      IF (INF.LT.3) RETURN
C
      IF (IPRINT.GT.1) THEN
        WRITE(LP,'(A30)') 'SPECIAL BRANCH OF SQP'
        WRITE(LP,'(A9,I5)') 'INF = ',INF
      ENDIF
C
      CTEMP = 0.0D0
      DO I=1,M
        IF (I.LE.MEQ .OR. C(I).LT.0.0D0) THEN
          IF (CTEMP.LT.ABS(C(I))) CTEMP = ABS(C(I))
        ENDIF
      ENDDO
C
      IF (JVMU == 0 .AND. CTEMP.GT.GMAX1) THEN
C       VMU HAS NOT BEEN UPDATED
        JVMU=1
        DO I=1,M
          VMU(I)=1.0D0
        ENDDO
        WOPT=FOPT
        DO I=1,M
          IF (I.LE.MEQ) THEN
            WOPT=WOPT+VMU(I)*ABS(COPT(I))
          ELSE
            WOPT=WOPT+VMU(I)*DMAX1(0.0D0,-COPT(I))
          ENDIF
        ENDDO
        WVAL=WOPT
        F=FOPT
        IndexDesign0 = INDEXD1
        DO I=1,N
          X(I)=XOPT(I)
        ENDDO
        DO I=1,M
          C(I)=COPT(I)
        ENDDO
        DO I=1,NP
          DO J=1,NP
            B(I,J)=0.0D0
          ENDDO
          B(I,I)=1.0D0
        ENDDO
        B(NP,NP)=0.0D0
        INF     = 0
        INFO    = -2
        icall   = 2
        NWDOGL  = ITER
        ITERP   = MIN0(IPRINT-1,0)
        ISWDOG  = 0
        call zerod(GLFIN,N)
        call zerod(GAMMA,N)
        call zerod(VLAM,M+2*N)
        call zerod(W,IWDIM)
        call zeroi(IW,NP)
        GOTO 90
      ENDIF
C
      IF (INF == 3 .AND. ISCREEN == 0 .AND. WVAL.LE.WOPT) THEN
        IF (CTEMP.GT.GMAX1) THEN
          CALL TRYVMU(0,FOPT,FBEG,FPRED,F,WOPT,WBEG,WPRED,WVAL,
     *                COPT,C,VMU,M,MEQ)
        ENDIF
        TEMP = (WBEG-WVAL)/(ABS(WBEG)+1.0D-6)
        IF (TEMP.GT.1.0D-2) THEN
          INF    = 0
          INFO   = 0
          ICALL  = 2
          ISWDOG = 0
          RETURN
        ENDIF
      ENDIF
C
      NMB=NMB+1
      IF (NMB.LE.20) THEN
        K = 0
        IF (NMB == 1) THEN
          K = 1
        ELSE
          IF (CBK.GT.1.0D-3) THEN
            TEMP = (CBK-CTEMP)/(ABS(CBK)+1.0D-6)
            IF (TEMP.GT.1.0D-2) K = 1
          ELSE
            TEMP = (FBK-F)/(ABS(FBK)+1.0D-6)
            IF (TEMP.GT.1.0D-2) K = 1
          ENDIF
        ENDIF
C
        IF ((K == 1).OR.((ISCREEN == 0).AND.(CTEMP.GT.GMAX1))) THEN
          FBK = F
          CBK = CTEMP
          DO I=1,NP
            DO J=1,NP
              B(I,J)=0.0D0
            ENDDO
            B(I,I)=1.0D0
          ENDDO
          B(NP,NP)=0.0D0
          WVAL = WOPT
          F    = FOPT
          IndexDesign0 = INDEXD1
          DO I=1,N
            X(I) = XOPT(I)
          ENDDO
          DO I=1,M
            C(I) = COPT(I)
          ENDDO
          INF     = 0
          INFO    = -2
          icall   = 2
          NWDOGL  = ITER
          ITERP   = MIN0(IPRINT-1,0)
          ISWDOG  = 0
          call zerod(GLFIN,N)
          call zerod(GAMMA,N)
          call zerod(VLAM,M+2*N)
          call zerod(W,IWDIM)
          call zeroi(IW,NP)
          IF (K == 1) THEN
            ISCREEN = 0
            IF (CTEMP.GT.GMAX1) THEN
              CALL TRYVMU(1,FOPT,FBEG,FPRED,F,WOPT,WBEG,WPRED,WVAL,
     *                    COPT,C,VMU,M,MEQ)
            ENDIF
            GOTO 90
          ELSE
C           constrait screening branch for failed cases
            ISCREEN = ISCREEN+1
            GOTO 85
          ENDIF
        ENDIF
      ENDIF
C
      IF (M.GT.0) THEN
        K = 0
        IF (CTEMP.GT.GMAX1) THEN
          IF (CMAXBEST.LT.CTEMP) THEN
            TEMP=(CTEMP-CMAXBEST)/CTEMP
            IF (TEMP.GT.1.0D-3) K=1
          ENDIF
        ELSE
          IF (FBEST.LT.F) THEN
            TEMP=(F-FBEST)/(ABS(F)+1.0D-3)
            IF (TEMP.GT.1.0D-3) K=1
          ENDIF
        ENDIF
        IF (K == 1) THEN
          INF=0
          INFO=-3
          ISCREEN=0
          IndexDesign0 = INDEXD2
          DO I=1,N
            X(I)=XBEST(I)
          ENDDO
          GOTO 300
        ENDIF
      ENDIF
C
      RETURN
      END
C
Chd|====================================================================
Chd|  TRYVMU                        src/mater1/mincf.F            
Chd|-- called by -----------
Chd|        VF13BD                        src/mater1/mincf.F            
Chd|-- calls ---------------
Chd|====================================================================
      SUBROUTINE TRYVMU(K,FOPT,FBEG,FPRED,F,WOPT,WBEG,WPRED,WVAL,
     &                   COPT,C,VMU,M,ME)
      IMPLICIT NONE
      INTEGER K,M,ME
      DOUBLE PRECISION FOPT,FBEG,FPRED,F,WOPT,WBEG,WPRED,WVAL,
     &         COPT(M),C(M),VMU(M)
C
      INTEGER I,J
      DOUBLE PRECISION CB,CP,CV,CT,WB,WP,WV,AUX0,AUX,VMAX,CMAX,CMAXO
C
      VMAX = VMU(1)
      DO I=2,M
        IF (VMU(I).GT.VMAX) VMAX=VMU(I)
      ENDDO
      IF (VMAX.GE.1.0D10) RETURN
C
      IF (K == 0) THEN
C
C       TRY WHETHER AUX CAN BE IMPROVED BY INCREASING VMU
C
        IF (M.GT.0 .AND. WPRED.LT.WBEG) THEN
          AUX0 = (WBEG-WVAL)/(WBEG-WPRED)
          CT = 2.0D0
          CB = WBEG - FBEG
          CP = WPRED - FPRED
          CV = WVAL - F
C         UPDATE
10        WB = FBEG + CT*CB
          WP = FPRED + CT*CP
          WV = F + CT*CV
          IF (WP.LT.WB) THEN
            AUX=(WB-WV)/(WB-WP)
            IF (AUX.GT.AUX0) THEN
              AUX0 = AUX
              CT = 2.0D0*CT
              IF (CT.LT.1.0D2) GOTO 10
            ENDIF
          ENDIF
          CT = 0.5D0*CT
          IF (CT.GT.1.0D0) THEN
            DO I=1,M
              VMU(I) = CT*VMU(I)
            ENDDO
            WBEG  = FBEG + CT*CB
            WPRED = FPRED + CT*CP
            WVAL  = F + CT*CV
            WOPT  = FOPT+CT*(WOPT-FOPT)
          ENDIF
        ENDIF
      ELSE IF (K == 1) THEN
C
C       INCREASE VMU TO ENSURE SMALLER VIOLATION
C
        CT = 1.05D0
        CB = WBEG - FBEG
        CP = WPRED - FPRED
        CV = WVAL - F
        WBEG  = FBEG + CT*CB
        WPRED = FPRED + CT*CP
        WVAL  = F + CT*CV
        WOPT  = FOPT+CT*(WOPT-FOPT)
        DO I=1,M
          VMU(I) = CT*VMU(I)
        ENDDO
      ELSE IF (K == 2) THEN
        WOPT = FOPT
        DO I=1,M
          VMU(I)=VMU(I)*1.0D2
          IF (VMU(I).LT.ABS(FOPT)) VMU(I)=ABS(FOPT)
          IF (I.LE.ME) THEN
            WOPT = WOPT+VMU(I)*ABS(COPT(I))
          ELSE
            WOPT = WOPT+VMU(I)*DMAX1(0.0D0,-COPT(I))
          ENDIF
        ENDDO
        WBEG = WOPT
        WVAL = WOPT
      ELSE IF (K == 3) THEN
C       AFTER ISCREEN (REDUCE VIOLATION), UPDATE VMU, WVAL AND WOPT
        WVAL=F
        WOPT=FOPT
        DO I=1,M
          IF (I.LE.ME) THEN
            WOPT=WOPT+VMU(I)*ABS(COPT(I))
            WVAL=WVAL+VMU(I)*ABS(C(I))
          ELSE
            WOPT=WOPT+VMU(I)*DMAX1(0.0D0,-COPT(I))
            WVAL=WVAL+VMU(I)*DMAX1(0.0D0,-C(I))
          ENDIF
        ENDDO
        IF (WVAL.GE.WOPT) THEN
          CV = WVAL-F
          CB = WOPT-FOPT
          IF (CV.LT.CB) THEN
            CT = 1.05D0*(F-FOPT)/(CB-CV)
            IF (CT.GT.1.0D3) CT=1.0D3
            IF (CT.LT.1.0D0) CT=1.0D0
            DO I=1,M
              VMU(I) = CT*VMU(I)
            ENDDO
            WVAL=F+CT*CV
            WOPT=FOPT+CT*CB
            WBEG=WOPT
          ELSE
            J=0
            CMAXO=0.0D0
            DO I=1,M
              IF ((I.LE.ME).OR.(COPT(I).LT.0.0D0)) THEN
                IF (CMAXO.LT.ABS(COPT(I))) THEN
                  J=I
                  CMAXO=ABS(COPT(I))
                ENDIF
              ENDIF
            ENDDO
            IF (CMAXO.GT.0.0D0) THEN
              IF ((J.LE.ME).OR.(C(J).LT.0.0D0)) THEN
                CMAX = ABS(C(J))
              ELSE
                CMAX = 0.0D0
              ENDIF
              WOPT = WOPT-VMU(J)*CMAXO
              WVAL = WVAL-VMU(J)*CMAX
              IF (VMU(J).LT.1.0D0) VMU(J)=1.0D0
              CB = VMU(J)*CMAXO
              CV = VMU(J)*CMAX
              IF (CV.LT.CB) THEN
                CT = 1.05D0*(WVAL-WOPT)/(CB-CV)
                IF (CT.GT.1.0D3) CT=1.0D3
                IF (CT.LT.1.0D0) CT=1.0D0
                VMU(J)=CT*VMU(J)
                WOPT = WOPT+VMU(J)*CMAXO
                WVAL = WVAL+VMU(J)*CMAX
                WBEG = WOPT
              ENDIF
            ENDIF ! IF (CMAXO.GT.0.0D0)
          ENDIF
        ENDIF
      ENDIF
C
      RETURN
      END




c=======================================================================
Chd|====================================================================
Chd|  VF13AD                        src/mater1/mincf.F            
Chd|-- called by -----------
Chd|        MINCF                         src/mater1/mincf.F            
Chd|-- calls ---------------
Chd|        VF13BD                        src/mater1/mincf.F            
Chd|====================================================================
      SUBROUTINE VF13AD (dvarl,dvaru,N, M, NDC, MEQ, X, F, G, C, CN,
     .              LCN, MAXFUN, MAXLIN, ACC, INF,W, LW, IW, icall,
     .              step_mincf,GMAX1)
c=======================================================================
C*   REVISED F77 VERSION OF VMCWD WITH REVERSE COMMUNICATION *
C*   MORE EFFICIENT QUADRATIC PROGRAMMING SUBROUTINE       *
C*   AND NO RESTRICTIONS ON SIZE OF PROBLEM   (APR 1987)   *
C***********************************************************
C
C     SUBROUTINE VF13AD CALCULATES THE LEAST VALUE OF A FUNCTION OF
C      SEVERAL VARIABLES SUBJECT TO EQUALITY AND INEQUALITY CONSTRAINTS
C      IT MUST BE CALLED REPEATEDLY BY THE USER UNTIL A RETURN
C      WITH INF.GT.0 TAKES PLACE. BEFORE THE INITIAL CALL INF
C      MUST BE SET NEGATIVE AND BEFORE SUBSEQUENT CALLS IT MUST BE
C      LEFT WITH VALUE ZERO.
C     AFTER PARTITIONING THE WORKSPACE SUBROUTINE VF13AD CALLS
C     SUBROUTINE  VF13BD WHICH PERFORMS THE MAJOR WORK OF AN ITERATION
C     SUBROUTINE  VF13BD ALSO CALLS
C                 SUBROUTINE VF13CD
C     A SUITABLE VERSION OF THIS SUBROUTINE IS GIVEN BELOW
C     IT DEPENDS ON THE HARWELL LIBRARY ROUTINE 'VE17AD'
C     THE PURPOSE OF VF13CD IS TO SOLVE THE QUADRATIC PROGRAMMING
C      PROBLEM THAT DEFINES THE SEARCH DIRECTION OF EACH ITERATION OF
C      VF13AD AND THAT PROVIDES ESTIMATES OF THE LAGRANGIAN PARAMETERS
C      OF THE CONSTRAINT FUNCTIONS OF THE MAIN CALCULATION
C
      IMPLICIT NONE
      INTEGER N, M, NDC, MEQ, INF, LW, icall,MAXLIN
      integer NP,MNN,IVMU,IVLAM,IXOPT,IB,ICBEG,IGRDCD,
     *         IGAMMA,IGLFIN,IXBEG,IBDELT,IDELTA,IGRAD,IH,IWDIM,MAXFUN,
     *         LCN,ICOPT,IXBEST

      DOUBLE PRECISION   X(N), F, G(N), C(NDC), CN(LCN,NDC), ACC, W(LW)
      DOUBLE PRECISION   dvarl(N),dvaru(N),step_mincf,GMAX1
      INTEGER IW(N+1)
c       integer lenvec, ivec(lenvec)
c     double precision AK(N),BK(N)
c-----------------------------------------------------------------------
C     N IS THE NUMBER OF VARIABLES
C     M IS THE TOTAL NUMBER OF CONSTRAINTS
C     MEQ IS THE NUMBER OF EQUALITY CONSTRAINTS
C     X IS THE VECTOR OF VARIABLES
C       IT MUST BE SET BY THE USER BEFORE THE INITIAL CALL AND LEFT
C       UNCHANGED THEREAFTER.
C     F IS THE VALUE OF THE OBJECTIVE FUNCTION
C     G IS THE GRADIENT OF THE OBJECTIVE FUNCTION
C     C IS THE VECTOR OF CONSTRAINT FUNCTIONS
C     CN IS THE MATRIX OF CONSTRAINT NORMALS
C     LCN IS THE FIRST DIMENSION OF CN
C     F,G,C,CN MUST ALL BE SET BY THE USER BEFORE EACH CALL.
C     MAXFUN BOUNDS THE NUMBER OF CALLS OF VF13AD
C     ACC CONTROLS THE FINAL ACCURACY - THE CALCULATION ENDS WHEN THE
C      OBJECTIVE FUNCTION PLUS SUITABLY WEIGHTED MULTIPLES OF THE
C      CONSTRAINT FUNCTIONS ARE PREDICTED TO DIFFER FROM THEIR
C      OPTIMAL VALUES BY AT MOST ACC
C     INF CONTROLS THE CALCULATION
C          =NEGATIVE ON INITIAL CALL
C          =0 DURING CALCULATION
C          =1 WHEN REQUIRED ACCURACY IS ACHIEVED
C          =2 WHEN VF13AD IS CALLED MAXFUN TIMES
C          =3 WHEN A LINE SEARCH REQUIRES 20 CALLS OF VF13AD
C          =4 WHEN AN UPHILL SEARCH DIRECTION IS CALCULATED
C          =5 WHEN NO FEASIBLE POINT IS FOUND BY VF13CD
C          =6 WHEN POSITIVE DEFINITENESS IS LOST
C          =7 WHEN LW IS TOO SMALL
C          =8 WHEN N,M OR MEQ HAS A SILLY VALUE
C     INF ALSO ALLOWS OPTIONS THAT CONTROL STEP-LENGTHS AND CONSTRAINT
C       VIOLATIONS.  BOUNDS ON THE CHANGES IN THE COMPONENTS OF X
C       MUST BE SET BY THE USER IF INF=-101 OR INF=-111 INITIALLY,
C       AND INITIAL VALUES OF THE PARAMETERS OF THE LINE SEARCH
C       OBJECTIVE FUNCTION MUST BE SET BY THE USER IF INF=-110 OR
C       INF=-111 INITIALLY - SEE THE COMMENTS ON THE ARRAY W BELOW
C     W IS A WORKSPACE ARRAY OF LENGTH LW
C       IF INF=-101 OR INF=-111 INITIALLY THEN THE USER MUST SET
C       W(M+I) TO AN UPPER BOUND ON THE CHANGE IN X(I) PER ITERATION
C       FOR I=1,2,...,N
C       IF INF=-110 OR INF=-111 INITIALLY THEN THE USER MUST SET THE
C       VALUES OF THE PARAMETERS OF THE LINE SEARCH OBJECTIVE FUNCTION
C       IN W(K) FOR K=1,2,...,M
C     IW IS AN INTEGER WORKSPACE ARRAY OF LENGTH N+1
C
C
C     VF13ED IS USED TO HOLD PRIVATE VARIABLES WHOSE VALUES
C       MUST BE PRESERVED ON SUBSEQUENT CALLS
C
c-----------------------------------------------------------------------
c
      MNN   =M+2*N
      NP    =N+1
      IVMU  =1
      IVLAM =IVMU+MNN
      IXOPT =IVLAM+MNN
      ICOPT =IXOPT+N
      IB    =ICOPT+NDC
      ICBEG =IB+NP*NP
      IGRDCD=ICBEG+MNN
      IGAMMA=IGRDCD+MNN
      IGLFIN=IGAMMA+N
      IXBEG =IGLFIN+N
      IBDELT=IXBEG+N
      IDELTA=IBDELT+N
      IGRAD=IDELTA+NP
      IXBEST=IGRAD+NP
      IH=IXBEST+N
      IWDIM=M+(NP*(3*NP+15))/2
c-----------------------------------------------------------------------
      CALL VF13BD (N, M, NDC, MEQ, X, F, G, C, CN, LCN, MAXFUN,
     1  MAXLIN, ACC, INF, W(IH), IWDIM, IW, NP, W(IVMU),
     2  dvarl, dvaru, W(IVLAM), W(IXOPT), W(ICOPT), W(IB),
     3   W(IDELTA), W(ICBEG), W(IGRDCD), W(IGAMMA), W(IGLFIN),
     4    W(IXBEG), W(IBDELT), W(IGRAD), W(IXBEST), icall,
     5  step_mincf,GMAX1 )
c-----------
      RETURN
      END
C
C-----------------------------------------------------------------------
C
Chd|====================================================================
Chd|  VF13CD                        src/mater1/mincf.F            
Chd|-- called by -----------
Chd|        VF13BD                        src/mater1/mincf.F            
Chd|-- calls ---------------
Chd|        VE17AD                        src/mater1/mincf.F            
Chd|====================================================================
      SUBROUTINE VF13CD (N, M, NDC, MEQ, G, C, CN, LCN, NP,
     1                   VLAM, B, DELTA, dvarl, dvaru, XVAR,
     2                   GRAD, W, IWDIM, IW,INFO)
C
      IMPLICIT NONE
      INTEGER N,M,NDC,MEQ,LCN,NP,IWDIM,INFO
      DOUBLE PRECISION G(N), C(NDC), CN(LCN,NDC)
      DOUBLE PRECISION VLAM(M+2*N), DELTA(NP), B(NP,NP), GRAD(NP),
     . W(IWDIM)
      DOUBLE PRECISION dvarl(N), dvaru(N), XVAR(N)
      INTEGER IW(NP)
C
C     THE PRIVATE ARRAYS THAT DO NOT OCCUR IN VF13AD/BD ARE AS FOLLOWS
C     dvarl AND dvaru GIVE LOWER AND UPPER BOUNDS ON VARIABLES
C     W AND IW ARE USED AS WORKING SPACE BY VE17AD
C
      DOUBLE PRECISION  FEASP, VLARGE, VSMALL, LBD
C
C     SET SOME PARAMETERS THAT ARE USED BY VF13CD AND VE17AD
C     FEASP IS A SCALING FACTOR THAT IS USED TO ACHIEVE FEASIBILITY
C     VLARGE IS ANY LARGE NUMBER
C     VSMALL IS USED TO DECIDE WHETHER THE CONSTRAINTS ARE INCONSISTENT
C
      DOUBLE PRECISION DIAG
      INTEGER I,J,K,IFLAG,INFOA,MTOTAL,MACT
C
C     /* 'DIAG' DECLARED DOUBLE PRESISION FOR VE17AD CALL - pef 920825 */
C
C     MACT WILL BE SET TO THE NUMBER OF ACTIVE CONSTRAINTS
C
      MACT=0
      INFOA=0
      FEASP=0.9D0
      VLARGE=1.0D10
      VSMALL=1.0D-6    
      LBD=1.0D0 !set limit to the expanded variable
      DIAG=0.D0
      IF (INFO .GE. 0) GO TO 50
C
C     EXTEND B BECAUSE OF THE EXTRA VARIABLE THAT IS INTRODUCED TO
C       ALLOW FOR INCONSISTENT LINEAR QP CONSTRAINTS
C
      DO 40 I=1,NP
      B(I,NP)=0.0D0
   40 B(NP,I)=0.0D0
C
C     SET THE ELEMENTS OF GRAD, CRES AND CN(NP,*)
C
   50 DO 60 I=1,N
   60 GRAD(I)=G(I)
      GRAD(NP)=-VLARGE
      DO 80 K=1,M
      IF (K .LE. MEQ) GO TO 70
      IF (C(K) .LT. 0.0D0) GO TO 70
      VLAM(K)=-C(K)
      CN(NP,K)=0.0D0
      GO TO 80
   70 VLAM(K)=0.0D0
      CN(NP,K)=C(K)
   80 CONTINUE
C
C     CALL SUBROUTINE VE17AD
C
      IFLAG=-1
c      BDU(NP)=1.0D0
   90 CALL VE17AD (N, NP,M, NDC, MEQ, dvarl, dvaru, LBD, XVAR, CN, LCN,
     1            VLAM,GRAD, B, NP, DELTA, MACT,IW,INFOA,DIAG,W,IWDIM)
C
C     BRANCH TO LABEL 160 IF THE CONSTRAINTS SEEM INCONSISTENT
C     MAKE ONE EXTRA CALL OF VE17AD IF THE EXTRA VARIABLE IS
C       NOT AT ITS BOUND
C
      IF (INFOA == 111) THEN
        INFO=111
        RETURN
      ENDIF
C
      INFO=6
      IF (DELTA(NP) .LE. VSMALL) GO TO 160
      K=0
      MTOTAL=M+2*NP
      GO TO 110
  109 IF (IW(K)  ==  MTOTAL) GO TO 120
  110 K=K+1
      IF (K .LE. MACT) GO TO 109
      IF (IFLAG  ==  0) GO TO 160
      IFLAG=0
      LBD=FEASP*DELTA(NP)
      GO TO 90
C
C     CALCULATE THE LAGRANGE MULTIPLIERS
C
  120 DO 130 K=1,M+2*N
  130 VLAM(K)=0.0D0
      DO 150 J=1,MACT
      K=IW(J)
      IF ((K == M+NP).OR.(K == M+2*NP)) GOTO 150
      IF (K .GT. M+N) THEN
        K=M+(K-M-NP)*2
      ELSE IF(K .GT. M) THEN
        K=M+(K-M)*2-1
      ENDIF
      VLAM(K)=W(J)
  150 CONTINUE
      GO TO 170
C
  160 INFO=5
  170 CONTINUE
c-----------
      RETURN
      END
C
C
Chd|====================================================================
Chd|  VE17AD                        src/mater1/mincf.F            
Chd|-- called by -----------
Chd|        VF13CD                        src/mater1/mincf.F            
Chd|-- calls ---------------
Chd|====================================================================
      SUBROUTINE VE17AD (N,NP,M,NDC,MEQ,dvarl,dvaru,LBD,XVAR,A,IDIMA,B,
     1                   GRAD,G,IDIMG,X,NACT,IACT,INFO,DIAG,W,IWDIM)
      IMPLICIT NONE
      INTEGER N,NP,M,NDC,MEQ,IDIMA,IDIMG,NACT,INFO,IWDIM
      INTEGER IACT(NP)
      DOUBLE PRECISION dvarl(N),dvaru(N),A(IDIMA,NDC),B(NDC),GRAD(NP),
     1 G(IDIMG,NP),X(NP),DIAG,W(IWDIM),XVAR(N),LBD
C
C     THIS SUBROUTINE CALCULATES THE NP-COMPONENT VECTOR X(.)
C       THAT MINIMIZES THE OBJECTIVE FUNCTION
C            GRAD(.).X(.) + 0.5*X(.).G(.,.).X(.),
C       SUBJECT TO
C            A(.,K).X(.)  ==  B(K)     K=1,2,...,MEQ,
C            A(.,K).X(.) .GE. B(K)     K=MEQ+1,...,M,
C            0 .LE. X(NP) .LE. z       I=1,...,NP.
C     NP IS THE NUMBER OF VARIABLES AND MUST BE POSITIVE.
C     M IS THE TOTAL NUMBER OF CONSTRAINTS.
C     MEQ IS THE NUMBER OF EQUALITY CONSTRAINTS.
C     dvarl(.) MUST CONTAIN THE LOWER BOUNDS ON XVAR(.).
C     dvaru(.) MUST CONTAIN THE UPPER BOUNDS ON XVAR(.).
C     A(.,.) IS A MATRIX WHOSE COLUMNS ARE THE CONSTRAINT NORMALS.
C     IDIMA IS THE FIRST DIMENSION OF THE ARRAY A(.,.) SO ITS VALUE
C       MUST BE AT LEAST NP.
C     B(.) MUST CONTAIN THE RIGHT HAND SIDES OF THE CONSTRAINTS.
C     GRAD(.) MUST CONTAIN THE GRADIENT OF F AT X=0.
C     G(.,.) MUST BE THE SECOND DERIVATIVE MATRIX OF F, AND MUST
C       BE SYMMETRIC.
C     IDIMG IS THE FIRST DIMENSION OF THE ARRAY G(.,.) SO ITS VALUE
C       MUST BE AT LEAST NP.
C     X(.) IS SET TO THE FINAL VECTOR OF VARIABLES.
C     NACT IS THE FINAL NUMBER OF ACTIVE CONSTRAINTS.
C     IACT(K) (K=1,2,...,NACT) ARE THE INDICES OF THE FINAL
C       ACTIVE CONSTRAINTS.  THE INDICES OF THE LOWER AND UPPER
C       BOUND CONSTRAINTS ARE (M+1,...,M+NP) AND (M+NP+1,...,M+2*NP)
C       RESPECTIVELY.
C     INFO INDICATES THE REASON FOR THE RETURN FROM THE SUBROUTINE.
C       THE FOLLOWING CASES CAN OCCUR:
C         INFO>0 INDICATES THAT THE CALCULATION IS SUCCESSFUL, THE
C       VALUE OF INFO BEING THE NUMBER OF ITERATIONS.
C         INFO=0 INDICATES THAT THE ACCURACY IS INSUFFICIENT TO
C       MAINTAIN INCREASING FUNCTION VALUES.
C         INFO<0 INDICATES THAT THE CONSTRAINT WITH INDEX ABS(INFO)
C       AND THE CONSTRAINTS WHOSE INDICES ARE IACT(K) (K=1,2,...,
C       NACT) ARE INCONSISTENT.
C     ON RETURN DIAG IS EQUAL TO THE MULTIPLE OF THE UNIT MATRIX
C        THAT WAS ADDED TO G(.,.) TO GIVE POSITIVE DEFINITENESS.
C     THE ELEMENTS OF W(.) ARE USED FOR WORKING SPACE.  THE LENGTH
C       OF W(.) MUST NOT BE LESS THAN ( 1.5*NP*NP + 7.5*NP + M ).
C     WHEN INFO>0 ON RETURN, THE LAGRANGE MULTIPLIERS OF THE FINAL
C       ACTIVE CONSTRAINTS ARE HELD IN W(K) (K=1,2,...,NACT).
C     THE VALUES OF NP, M, MEQ, IDIMA AND IDIMG AND THE ELEMENTS OF
C       A, B, GRAD AND G ARE NOT ALTERED BY THE CALCULATION.
C
      DOUBLE PRECISION CVMAX,DIAGR,FDIFF,FDIFFA,GA,GB,PARINC,PARNEW,
     *     RATIO,RES,STEP,SUM,SUMA,SUMB,SUMC,SUMX,SUMY,
     *     TEMP,TEMPA,VSMALL,VLARGE,XMAG,XMAGR,ZERO,ONE
C
      INTEGER I,J,K,IWZ,IWR,IWW,IWD,IWX,IWA,IWB,IFINC,KFINC,IA,IB,
     *         ID,II,IR,IRA,IRB,NM,IZ,IZA,ITERC,ITREF,JFINC,IFLAG,IWS,
     *         IW,IS,KK,IL,IU,JU,KFLAG,LFLAG,JFLAG,KDROP,
     *         NU,MFLAG,KM,KNEXT,IX,NFLAG,IWY,IY,JL
C
C     THE FOLLOWING INTEGERS ARE USED TO PARTITION W(.).
C     THE FIRST NP ELEMENTS OF W HOLD LAGRANGE MULTIPLIER ESTIMATES.
C     W(IWZ+I+(NP-1)*J) HOLDS THE MATRIX ELEMENT Z(I,J).
C     W(IWR+I+0.5*J*(J-1)) HOLDS THE UPPER TRIANGULAR MATRIX
C       ELEMENT R(I,J).  THE SUBSEQUENT NP COMPONENTS OF W MAY BE
C       TREATED AS AN EXTRA COLUMN OF R(.,.).
C     W(IWW-NP+I) (I=1,2,...,NP) ARE USED FOR TEMPORARY STORAGE.
C     W(IWW+I) (I=1,2,...,NP) ARE USED FOR TEMPORARY STORAGE.
C     W(IWD+I) (I=1,2,...,NP) HOLDS G(I,I) DURING THE CALCULATION.
C     W(IWX+I) (I=1,2,...,NP) HOLDS VARIABLES THAT WILL BE USED TO
C       TEST THAT THE ITERATIONS INCREASE THE OBJECTIVE FUNCTION.
C     W(IWA+K) (K=1,2,...,M) USUALLY HOLDS THE RECIPROCAL OF THE
C       LENGTH OF THE K-TH CONSTRAINT, BUT ITS SIGN INDICATES
C       WHETHER THE CONSTRAINT IS ACTIVE.
C     W(IWB+K) (K=1,2,...,2*NP) EXTENDS THE PREVIOUS PARTITION OF
C       W(.) FOR THE BOUND CONSTRAINTS.
C
      IWZ=NP
      IWR=IWZ+NP*NP
      IWW=IWR+(NP*(NP+3))/2
      IWD=IWW+NP
      IWX=IWD+NP
      IWA=IWX+NP
      IWB=IWA+M
C
C     SET SOME PARAMETERS.
C     NUMBERS LESS THAN VSMALL ARE ASSUMED TO BE NEGLIGIBLE.
C     VLARGE PREVENTS AN OVERFLOW FOR LARGE X(.).
C     THE MULTIPLE OF I THAT IS ADDED TO G IS AT MOST DIAGR TIMES
C       THE LEAST MULTIPLE OF I THAT GIVES POSITIVE DEFINITENESS.
C     X IS RE-INITIALISED IF ITS MAGNITUDE IS REDUCED BY THE
C       FACTOR XMAGR.
C     A CHECK IS MADE FOR AN INCREASE IN F EVERY IFINC
C       ITERATIONS, AFTER KFINC ITERATIONS ARE COMPLETED.
C
      RATIO=0.0D0
      ZERO=0.0D0
      ONE=1.0D0
      VSMALL=1.D-30
      VLARGE=1.D20
      DIAGR=2.0D0
      XMAGR=0.01D0
      IFINC=3
      KFINC=MAX0(10,NP)
C
      SUMY=0.D0
      STEP=0.D0
      RES=0.D0
      PARNEW=0.D0
      PARINC=0.D0
      TEMP=0.D0
      NU=0
      NFLAG=0
      MFLAG=0
      KNEXT=0
      JFLAG=0
      KDROP=0
C
C     FIND THE RECIPROCALS OF THE LENGTHS OF THE CONSTRAINT NORMALS.
C     RETURN IF A CONSTRAINT IS INFEASIBLE DUE TO A ZERO NORMAL.
C
      NACT=0
      IF (M .LE. 0) GOTO 50
      DO 40 K=1,M
      SUM=ZERO
      DO 10 I=1,NP
   10 SUM=SUM+A(I,K)**2
      IF (SUM .GT. ZERO) GOTO 20
      IF (B(K)  ==  ZERO) GOTO 30
      INFO=-K
      IF (K .LE. MEQ) GOTO 730
      IF (B(K)) 30,30,730
   20 SUM=ONE/DSQRT(DMAX1(VSMALL,SUM))
   30 IA=IWA+K
   40 W(IA)=SUM
C
C     INCLUDE THE NORMALS OF THE BOUNDS.
C
   50 IA=IWB+1
      IB=IWB+NP+NP
      DO 53 I=IA,IB
   53 W(I)=ONE
C
C     IF NECESSARY INCREASE THE DIAGONAL ELEMENTS OF G.
C
      DIAG=ZERO
      DO 60 I=1,NP
      ID=IWD+I
      W(ID)=G(I,I)
      DIAG=DMAX1(DIAG,VSMALL-W(ID))
      IF(I == NP)GO TO 60
      II=I+1
      DO 55 J=II,NP
      GA=-DMIN1(W(ID),G(J,J))
      GB=DABS(W(ID)-G(J,J))+DABS(G(I,J))
      IF(GB.GT.VSMALL)GA=GA+G(I,J)**2/GB
   55 DIAG=DMAX1(DIAG,GA)
   60 CONTINUE
      IF (DIAG .LE. ZERO) GOTO 90
   70 DIAG=DIAGR*DIAG
      DO 80 I=1,NP
      ID=IWD+I
   80 G(I,I)=DIAG+W(ID)
C
C     FORM THE CHOLESKY FACTORISATION OF G.  THE TRANSPOSE
C       OF THE FACTOR WILL BE PLACED IN THE R-PARTITION OF W.
C
   90 IR=IWR
      DO 130 J=1,NP
      IRA=IWR
      IRB=IR+1
      DO 120 I=1,J
      TEMP=G(I,J)
      IF (I  ==  1) GOTO 110
      DO 100 K=IRB,IR
      IRA=IRA+1
  100 TEMP=TEMP-W(K)*W(IRA)
  110 IR=IR+1
      IRA=IRA+1
      IF (I .LT. J) THEN
        IF (DABS(W(IRA)).LT.VSMALL) GOTO 950
        W(IR)=TEMP/W(IRA)
      ENDIF
  120 CONTINUE
      IF (TEMP .LT. VSMALL) GOTO 140
  130 W(IR)=DSQRT(TEMP)
      GOTO 170
C
C     INCREASE FURTHER THE DIAGONAL ELEMENTS OF G.
C
  140 W(J)=ONE
      SUMX=ONE
      K=J
  150 SUM=ZERO
      IRA=IR-1
      DO 160 I=K,J
      SUM=SUM-W(IRA)*W(I)
  160 IRA=IRA+I
      IR=IR-K
      K=K-1
      IF (DABS(W(IR)).LT.VSMALL) GOTO 950
      W(K)=SUM/W(IR)
      SUMX=SUMX+W(K)**2
      IF (K .GE. 2) GOTO 150
      IF (DABS(SUMX).LT.VSMALL) GOTO 950
      DIAG=DIAG+VSMALL-TEMP/SUMX
      GOTO 70
C
C     SET Z TO THE INVERSE OF THE MATRIX IN R.
C
  170 NM=NP-1
      DO 220 I=1,NP
      IZ=IWZ+I
      IF (I  ==  1) GOTO 190
      DO 180 J=2,I
      W(IZ)=ZERO
  180 IZ=IZ+NP
  190 IR=IWR+(I+I*I)/2
      IF (DABS(W(IR)).LT.VSMALL) GOTO 950
      W(IZ)=ONE/W(IR)
      IF (I  ==  NP) GOTO 220
      IZA=IZ
      DO 210 J=I,NM
      IR=IR+I
      SUM=ZERO
      DO 200 K=IZA,IZ,NP
      SUM=SUM+W(K)*W(IR)
  200 IR=IR+1
      IZ=IZ+NP
      IF (DABS(W(IR)).LT.VSMALL) GOTO 950
  210 W(IZ)=-SUM/W(IR)
  220 CONTINUE
C
C     SET THE INITIAL VALUES OF SOME VARIABLES.
C     ITERC COUNTS THE NUMBER OF ITERATIONS.
C     ITREF IS SET TO ONE WHEN ITERATIVE REFINEMENT IS REQUIRED.
C     JFINC INDICATES WHEN TO TEST FOR AN INCREASE IN F.
C
      ITERC=1
      ITREF=0
      JFINC=-KFINC
C
C     SET X TO ZERO AND SET THE CORRESPONDING RESIDUALS OF THE
C       KUHN-TUCKER CONDITIONS.
C
  230 IFLAG=1
      IWS=IWW-NP
      DO 240 I=1,NP
      X(I)=ZERO
      IW=IWW+I
      W(IW)=GRAD(I)
      IF (I .GT. NACT) GOTO 240
      W(I)=ZERO
      IS=IWS+I
      K=IACT(I)
      IA=K-M
      IB=IA-NP
      IF (K .LE. M) W(IS)=B(K)
      IF (IA .GT. 0 .AND. IA .LE. N) W(IS)=dvarl(IA)-XVAR(IA)
      IF (IA .GT. N .AND. IA .LE. NP) W(IS)=0.0D0
      IF (IB .GT. 0 .AND. IB .LE. N) W(IS)=XVAR(IB)-dvaru(IB)
      IF (IB .GT. N .AND. IB .LE. NP) W(IS)=-LBD
  240 CONTINUE
      XMAG=ZERO
      IF (NACT) 340,340,280
C
C     SET THE RESIDUALS OF THE KUHN-TUCKER CONDITIONS FOR GENERAL X.
C
  250 IFLAG=2
      IWS=IWW-NP
      DO 260 I=1,NP
      IW=IWW+I
      W(IW)=GRAD(I)
      DO 260 J=1,NP
  260 W(IW)=W(IW)+G(I,J)*X(J)
      IF (NACT  ==  0) GOTO 340
      DO 277 K=1,NACT
      KK=IACT(K)
      IS=IWS+K
      IF (KK .GT. M) GOTO 273
      W(IS)=B(KK)
      DO 270 I=1,NP
      IW=IWW+I
      W(IW)=W(IW)-W(K)*A(I,KK)
  270 W(IS)=W(IS)-X(I)*A(I,KK)
      GOTO 277
  273 IA=KK-M
      IF (IA .GT. NP) GOTO 275
      IW=IWW+IA
      W(IW)=W(IW)-W(K)
      IF (IA .LE. N) THEN
         W(IS)=dvarl(IA)-XVAR(IA)-X(IA)
      ELSE
         W(IS)=-X(IA)
      ENDIF
      GOTO 277
  275 IB=IA-NP
      IW=IWW+IB
      W(IW)=W(IW)+W(K)
      IF (IB .LE. N) THEN
         W(IS)=X(IB)-(dvaru(IB)-XVAR(IB))
      ELSE
         W(IS)=X(IB)-LBD
      ENDIF
  277 CONTINUE
C
C     PRE-MULTIPLY THE VECTOR IN THE S-PARTITION OF W BY THE
C       INVERSE OF R TRANSPOSE.
C
  280 IR=IWR
      IL=IWS+1
      IU=IWS+NACT
      DO 310 I=IL,IU
      SUM=ZERO
      IF (I  ==  IL) GOTO 300
      JU=I-1
      DO 290 J=IL,JU
      IR=IR+1
  290 SUM=SUM+W(IR)*W(J)
  300 IR=IR+1
      IF (DABS(W(IR)).LT.VSMALL) GOTO 950
  310 W(I)=(W(I)-SUM)/W(IR)
C
C     SHIFT X TO SATISFY THE ACTIVE CONSTRAINTS AND MAKE THE
C       CORRESPONDING CHANGE TO THE GRADIENT RESIDUALS.
C
      DO 330 I=1,NP
      IZ=IWZ+I
      SUM=ZERO
      DO 320 J=IL,IU
      SUM=SUM+W(J)*W(IZ)
  320 IZ=IZ+NP
      X(I)=X(I)+SUM
      DO 330 J=1,NP
      IW=IWW+J
  330 W(IW)=W(IW)+SUM*G(I,J)
C
C     FORM THE SCALAR PRODUCT OF THE CURRENT GRADIENT RESIDUALS
C       WITH EACH COLUMN OF Z.
C
  340 KFLAG=1
      GOTO 930
  350 IF (NACT  ==  NP) GOTO 380
C
C     SHIFT X SO THAT IT SATISFIES THE REMAINING KUHN-TUCKER
C     CONDITIONS.
C     X=X-Z*Z'*GRAD
      IL=IWS+NACT+1
      IZA=IWZ+NACT*NP
      DO 370 I=1,NP
      SUM=ZERO
      IZ=IZA+I
      DO 360 J=IL,IWW
      SUM=SUM+W(IZ)*W(J)
  360 IZ=IZ+NP
  370 X(I)=X(I)-SUM
      INFO=ITERC
      IF (NACT  ==  0) GOTO 410
C
C     UPDATE THE LAGRANGE MULTIPLIERS.
C
  380 LFLAG=3
      GOTO 740
  390 DO 400 K=1,NACT
      IW=IWW+K
  400 W(K)=W(K)+W(IW)
C
C     REVISE THE VALUE OF XMAG.
C     BRANCH IF ITERATIVE REFINEMENT IS REQUIRED.
C
  410 JFLAG=1
      GOTO 910
  420 IF (IFLAG  ==  ITREF) GOTO 250
C
C     DELETE A CONSTRAINT IF A LAGRANGE MULTIPLIER OF AN
C       INEQUALITY CONSTRAINT IS NEGATIVE.
C
      KDROP=0
      GOTO 440
  430 KDROP=KDROP+1
      IF (W(KDROP) .GE. ZERO) GOTO 440
      IF (IACT(KDROP) .LE. MEQ) GOTO 440
      NU=NACT
      MFLAG=1
      GOTO 800
  440 IF (KDROP .LT. NACT) GOTO 430
C
C     SEEK THE GREATEST NORMALISED CONSTRAINT VIOLATION, DISREGARDING
C       ANY THAT MAY BE DUE TO COMPUTER ROUNDING ERRORS.
C     numd*(A'*d-b) reciprocal of constraints normal
  450 CVMAX=ZERO
      KM=M+NP+NP
      DO 480 K=1,KM
      IA=IWA+K
      IF (W(IA) .LE. ZERO) GOTO 480
      IF (K .GT. M) GOTO 473
      SUM=-B(K)
      DO 460 I=1,NP
  460 SUM=SUM+X(I)*A(I,K)
      SUMX=-SUM*W(IA)
      IF (K .LE. MEQ) SUMX=DABS(SUMX)
      IF (SUMX .LE. CVMAX) GOTO 480
      TEMP=DABS(B(K))
      DO 470 I=1,NP
  470 TEMP=TEMP+DABS(X(I)*A(I,K))
      GOTO 477
  473 IA=K-M
      IB=IA-NP
      IF (IA .LE. N) SUM=X(IA)-(dvarl(IA)-XVAR(IA))
      IF (IA .GT. N .AND. IA .LE. NP) SUM=X(IA)
      IF (IB .GT. 0 .AND. IB .LE. N) SUM=(dvaru(IB)-XVAR(IB))-X(IB)
      IF (IB .GT. N .AND. IB .LE. NP) SUM=LBD-X(IB)
      SUMX=-SUM
      IF (SUMX .LE. CVMAX) GOTO 480
      IF (IA .LE. N) TEMP=DABS(dvarl(IA)-XVAR(IA))+DABS(X(IA))
      IF (IA .GT. N .AND. IA .LE. NP) TEMP=DABS(X(IA))
      IF (IB .GT. 0 .AND. IB .LE. N) TEMP=DABS(X(IB))+
     1                                DABS(dvaru(IB)-XVAR(IB))
      IF (IB .GT. N .AND. IB .LE. NP) TEMP=DABS(X(IB))+ABS(LBD)
  477 TEMPA=TEMP+DABS(SUM)
      IF (TEMPA .LE. TEMP) GOTO 480
      TEMP=TEMP+1.5D0*DABS(SUM)
      IF (TEMP .LE. TEMPA) GOTO 480
      CVMAX=SUMX
      RES=SUM
      KNEXT=K
  480 CONTINUE
C
C     TEST FOR CONVERGENCE.
C
      INFO=ITERC
      IF (CVMAX .LE. VSMALL) GOTO 700
C
C     RETURN IF, DUE TO ROUNDING ERRORS, THE ACTUAL CHANGE IN
C       X MAY NOT INCREASE THE OBJECTIVE FUNCTION.
C
      JFINC=JFINC+1
      IF (JFINC  ==  0) GOTO 510
      IF (JFINC .NE. IFINC) GOTO 530
      FDIFF=ZERO
      FDIFFA=ZERO
      DO 500 I=1,NP
      SUM=2.0D0*GRAD(I)
      SUMX=DABS(SUM)
      DO 490 J=1,NP
      IX=IWX+J
      TEMP=G(I,J)*(W(IX)+X(J))
      SUM=SUM+TEMP
  490 SUMX=SUMX+DABS(TEMP)
      IX=IWX+I
      FDIFF=FDIFF+SUM*(X(I)-W(IX))
  500 FDIFFA=FDIFFA+SUMX*DABS(X(I)-W(IX))
      INFO=0
      SUM=FDIFFA+FDIFF
      IF (SUM .LE. FDIFFA) GOTO 700
      TEMP=FDIFFA+1.5D0*FDIFF
      IF (TEMP .LE. SUM) GOTO 700
      JFINC=0
  510 DO 520 I=1,NP
      IX=IWX+I
  520 W(IX)=X(I)
C
C     FORM THE SCALAR PRODUCT OF THE NEW CONSTRAINT NORMAL WITH EACH
C       COLUMN OF Z.  PARNEW WILL BECOME THE LAGRANGE MULTIPLIER OF
C       THE NEW CONSTRAINT.
C
  530 ITERC=ITERC+1
      DO 540 I=1,NP
      IW=IWW+I
      W(IW)=ZERO
      IF (KNEXT .LE. M) W(IW)=A(I,KNEXT)
  540 CONTINUE
      IWS=IWR+(NACT+NACT*NACT)/2
      KFLAG=2
      IF (KNEXT .LE. M) GOTO 930
      IA=KNEXT-M
      IB=IA-NP
      TEMP=ONE
      IF (IA .GT. NP) TEMP=-TEMP
      IF (IA .GT. NP) IA=IB
      IW=IWW+IA
      W(IW)=TEMP
      GOTO 930
  550 PARNEW=ZERO
C
C     APPLY GIVENS ROTATIONS TO MAKE THE LAST (NP-NACT-2) SCALAR
C       PRODUCTS EQUAL TO ZERO.
C
      IF (NACT  ==  NP) GOTO 570
      NU=NP
      NFLAG=1
      GO TO 860
C
C     BRANCH IF THERE IS NO NEED TO DELETE A CONSTRAINT.
C
  560 IS=IWS+NACT
      IF (NACT  ==  0) GOTO 640
      SUMA=ZERO
      SUMB=ZERO
      SUMC=ZERO
      IZ=IWZ+NACT*NP
      DO 563 I=1,NP
      IZ=IZ+1
      IW=IWW+I
      SUMA=SUMA+W(IW)*W(IZ)
      SUMB=SUMB+DABS(W(IW)*W(IZ))
  563 SUMC=SUMC+DABS(W(IZ))
      TEMP=SUMB+0.1D0*DABS(SUMA)
      TEMPA=SUMB+0.2D0*DABS(SUMA)
      IF(TEMP.LE.SUMB)GO TO 570
      IF(TEMPA.LE.TEMP)GO TO 570
      IA=IWA+KNEXT
      IF (DABS(W(IA)).LT.VSMALL) GOTO 950
      SUMC=SUMC/W(IA)
      TEMP=SUMC+0.1D0*DABS(SUMA)
      TEMPA=SUMC+0.2D0*DABS(SUMA)
      IF(TEMP.LE.SUMC)GO TO 567
      IF(TEMPA.LE.TEMP)GO TO 567
      GO TO 640
C
C     CALCULATE THE MULTIPLIERS FOR THE NEW CONSTRAINT NORMAL
C       EXPRESSED IN TERMS OF THE ACTIVE CONSTRAINT NORMALS.
C     THEN WORK OUT WHICH CONSTRAINT TO DROP.
C
  567 LFLAG=4
      GO TO 740
  570 LFLAG=1
      GOTO 740
C
C     COMPLETE THE TEST FOR LINEARLY DEPENDENT CONSTRAINTS.
C
  573 DO 577 I=1,NP
      IA=KNEXT-M
      IB=IA-NP
      SUMA=ZERO
      IF (KNEXT .LE. M) SUMA=A(I,KNEXT)
      IF (I  ==  IA) SUMA=ONE
      IF (I  ==  IB) SUMA=-ONE
      SUMB=DABS(SUMA)
      DO 575 K=1,NACT
      KK=IACT(K)
      IA=KK-M
      IB=IA-NP
      TEMP=ZERO
      IF (KK .LE. M) TEMP=A(I,KK)
      IF (I  ==  IA) TEMP=ONE
      IF (I  ==  IB) TEMP=-ONE
      IW=IWW+K
      TEMP=W(IW)*TEMP
      SUMA=SUMA-TEMP
      SUMB=SUMB+DABS(TEMP)
  575 CONTINUE
      TEMP=SUMB+0.1D0*DABS(SUMA)
      TEMPA=SUMB+0.2D0*DABS(SUMA)
      IF(TEMP.LE.SUMB)GO TO 577
      IF(TEMPA.LE.TEMP)GO TO 577
      GO TO 630
  577 CONTINUE
      LFLAG=1
      GO TO 775
C
C     BRANCH IF THE CONSTRAINTS ARE INCONSISTENT.
C
  580 INFO=-KNEXT
      IF (KDROP  ==  0) GOTO 700
      PARINC=RATIO
      PARNEW=PARINC
C
C     REVISE THE LAGRANGE MULTIPLIERS OF THE ACTIVE CONSTRAINTS.
C
  590 DO 600 K=1,NACT
      IW=IWW+K
      W(K)=W(K)-PARINC*W(IW)
      IF (IACT(K) .GT. MEQ) W(K)=DMAX1(ZERO,W(K))
  600 CONTINUE
      IF (KDROP  ==  0) GOTO 680
C
C     DELETE THE CONSTRAINT TO BE DROPPED.
C     SHIFT THE VECTOR OF SCALAR PRODUCTS.
C     THEN, IF APPROPRIATE, MAKE ONE MORE SCALAR PRODUCT ZERO.
C
      NU=NACT+1
      MFLAG=2
      GOTO 800
  610 IWS=IWS-NACT-1
      NU=MIN0(NP,NU)
      DO 620 I=1,NU
      IS=IWS+I
      J=IS+NACT
  620 W(IS)=W(J+1)
      NFLAG=2
      GOTO 860
C
C     CALCULATE THE STEP TO THE VIOLATED CONSTRAINT.
C
  630 IS=IWS+NACT
  640 SUMY=W(IS+1)
      IF (DABS(SUMY).LT.VSMALL) GOTO 950
      STEP=-RES/SUMY
      PARINC=STEP/SUMY
      IF (NACT  ==  0) GOTO 660
C
C     CALCULATE THE CHANGES TO THE LAGRANGE MULTIPLIERS, AND REDUCE
C     THE STEP ALONG THE NEW SEARCH DIRECTION IF NECESSARY.
C
      LFLAG=2
      GOTO 740
  650 IF (KDROP  ==  0) GOTO 660
      IF (DABS(PARINC).LT.VSMALL) GOTO 950
      TEMP=ONE-RATIO/PARINC
      IF (TEMP .LE. ZERO) KDROP=0
      IF (KDROP  ==  0) GOTO 660
      STEP=RATIO*SUMY
      PARINC=RATIO
      RES=TEMP*RES
C
C     UPDATE X AND THE LAGRANGE MULTIPLIERS.
C     DROP A CONSTRAINT IF THE FULL STEP IS NOT TAKEN.
C
  660 IWY=IWZ+NACT*NP
      DO 670 I=1,NP
      IY=IWY+I
      X(I)=X(I)+STEP*W(IY)
      if (.not.(X(I).lt.1.0d30 .and. X(I).gt.-1.0d30)) then
        INFO = 111                                 ! check NaN data
        return
      endif
  670 CONTINUE
      PARNEW=PARNEW+PARINC
      IF (NACT .GE. 1) GOTO 590
C
C     ADD THE NEW CONSTRAINT TO THE ACTIVE SET.
C
  680 NACT=NACT+1
      W(NACT)=PARNEW
      IACT(NACT)=KNEXT
      IA=IWA+KNEXT
      W(IA)=-W(IA)
C
C     ESTIMATE THE MAGNITUDE OF X.  THEN BEGIN A NEW ITERATION,
C       RE-INITIALISING X IF THIS MAGNITUDE IS SMALL.
C
      JFLAG=2
      GOTO 910
  690 IF (SUM .LT. XMAGR*XMAG) GOTO 230
      IF (ITREF) 450,450,250
C
C     INITIATE ITERATIVE REFINEMENT IF IT HAS NOT YET BEEN USED,
C       OR RETURN AFTER RESTORING THE DIAGONAL ELEMENTS OF G.
C     ALSO SATISFY ANY BOUND CONSTRAINTS EXACTLY.
C
  700 ITREF=ITREF+1
      JFINC=-1
      IF (ITREF  ==  1) GOTO 250
      DO 720 I=1,NP
      ID=IWD+I
      G(I,I)=W(ID)
      IF (I .GT. NACT) GOTO 720
      K=IACT(I)
      IF (K .LE. M) GOTO 720
      IA=K-M
      IB=IA-NP
      IF (IA .LE. N) X(IA)=dvarl(IA)-XVAR(IA)
      IF (IA .GT. N .AND. IA .LE. NP) X(IA)=0.0D0
      IF (IB .GT. 0 .AND. IB .LE. N) X(IB)=dvaru(IB)-XVAR(IB)
      IF (IB .GT. N) X(IB)=LBD
  720 CONTINUE
  730 GO TO 950
C
C     THE REMAINING INSTRUCTIONS ARE USED AS SUBROUTINES.
C
C*****************************************************************
C
C     CALCULATE THE LAGRANGE MULTIPLIERS BY PRE-MULTIPLYING THE
C       VECTOR IN THE S-PARTITION OF W BY THE INVERSE OF R.
C
  740 IR=IWR+(NACT+NACT*NACT)/2
      I=NACT
      SUM=ZERO
      GOTO 770
  750 IRA=IR-1
      SUM=ZERO
      DO 760 J=I,NACT
      IW=IWW+J
      SUM=SUM+W(IRA)*W(IW)
  760 IRA=IRA+J
      IR=IR-I
      I=I-1
  770 IW=IWW+I
      IS=IWS+I
      IF (DABS(W(IR)).LT.VSMALL) GOTO 950
      W(IW)=(W(IS)-SUM)/W(IR)
      IF (I .GT. 1) GOTO 750
      IF (LFLAG  ==  3) GOTO 390
      IF (LFLAG  ==  4) GOTO 573
C
C     CALCULATE THE NEXT CONSTRAINT TO DROP.
C
  775 KDROP=0
      DO 790 K=1,NACT
      IF (IACT(K) .LE. MEQ) GOTO 790
      IW=IWW+K
      IF (RES*W(IW) .GE. ZERO) GOTO 790
      IF (DABS(W(IW)).LT.VSMALL) GOTO 790
      TEMP=W(K)/W(IW)
      IF (KDROP  ==  0) GOTO 780
      IF (DABS(TEMP) .GE. DABS(RATIO)) GOTO 790
  780 KDROP=K
      RATIO=TEMP
  790 CONTINUE
      GOTO (580,650), LFLAG
C
C*****************************************************************
C
C     DROP THE CONSTRAINT IN POSITION KDROP IN THE ACTIVE SET.
C
  800 IA=IWA+IACT(KDROP)
      W(IA)=-W(IA)
      IF (KDROP  ==  NACT) GOTO 850
C
C     SET SOME INDICES AND CALCULATE THE ELEMENTS OF THE NEXT
C       GIVENS ROTATION.
C
      IZ=IWZ+KDROP*NP
      IR=IWR+(KDROP+KDROP*KDROP)/2
  810 IRA=IR
      IR=IR+KDROP+1
      TEMP=DMAX1(DABS(W(IR-1)),DABS(W(IR)))
      IF (DABS(TEMP).LT.VSMALL) GOTO 950
      SUM=TEMP*DSQRT((W(IR-1)/TEMP)**2+(W(IR)/TEMP)**2)
      IF (DABS(SUM).LT.VSMALL) GOTO 950
      GA=W(IR-1)/SUM
      GB=W(IR)/SUM
C
C     EXCHANGE THE COLUMNS OF R.
C
      DO 820 I=1,KDROP
      IRA=IRA+1
      J=IRA-KDROP
      TEMP=W(IRA)
      W(IRA)=W(J)
  820 W(J)=TEMP
      W(IR)=ZERO
C
C     APPLY THE ROTATION TO THE ROWS OF R.
C
      W(J)=SUM
      KDROP=KDROP+1
      DO 830 I=KDROP,NU
      TEMP=GA*W(IRA)+GB*W(IRA+1)
      W(IRA+1)=GA*W(IRA+1)-GB*W(IRA)
      W(IRA)=TEMP
  830 IRA=IRA+I
C
C     APPLY THE ROTATION TO THE COLUMNS OF Z.
C
      DO 840 I=1,NP
      IZ=IZ+1
      J=IZ-NP
      TEMP=GA*W(J)+GB*W(IZ)
      W(IZ)=GA*W(IZ)-GB*W(J)
  840 W(J)=TEMP
C
C     REVISE IACT AND THE LAGRANGE MULTIPLIERS.
C
      IACT(KDROP-1)=IACT(KDROP)
      W(KDROP-1)=W(KDROP)
      IF (KDROP .LT. NACT) GOTO 810
  850 NACT=NACT-1
      GOTO (250,610), MFLAG
C
C*****************************************************************
C
C     APPLY GIVENS ROTATIONS TO REDUCE SOME OF THE SCALAR
C       PRODUCTS IN THE S-PARTITION OF W TO ZERO.
C     [Z,Z'*A(*,KNEXT)] ROTATION
  860 IZ=IWZ+NU*NP
  870 IZ=IZ-NP
  880 IS=IWS+NU
      NU=NU-1
      IF (NU  ==  NACT) GOTO 900
      IF (W(IS)  ==  ZERO) GOTO 870
      TEMP=DMAX1(DABS(W(IS-1)),DABS(W(IS)))
      IF (DABS(TEMP).LT.VSMALL) GOTO 950
      SUM=TEMP*DSQRT((W(IS-1)/TEMP)**2+(W(IS)/TEMP)**2)
      IF (DABS(SUM).LT.VSMALL) GOTO 950
      GA=W(IS-1)/SUM
      GB=W(IS)/SUM
      W(IS-1)=SUM
      DO 890 I=1,NP
      K=IZ+NP
      TEMP=GA*W(IZ)+GB*W(K)
      W(K)=GA*W(K)-GB*W(IZ)
      W(IZ)=TEMP
  890 IZ=IZ-1
      GOTO 880
  900 GOTO (560,630), NFLAG
C
C*****************************************************************
C
C     CALCULATE THE MAGNITUDE OF X AND REVISE XMAG.
C
  910 SUM=ZERO
      DO 920 I=1,NP
      TEMP=DABS(X(I))
      IF (TEMP .LE. VLARGE) TEMP=TEMP*(DABS(GRAD(I))+
     1 DABS(G(I,I)*X(I)))
  920 SUM=SUM+TEMP
      XMAG=DMAX1(XMAG,SUM)
      GOTO (420,690), JFLAG
C
C*****************************************************************
C
C     PRE-MULTIPLY THE VECTOR IN THE W-PARTITION OF W BY Z TRANSPOSE.
C
  930 JL=IWW+1
      IZ=IWZ
      DO 940 I=1,NP
      IS=IWS+I
      W(IS)=ZERO
      DO 940 J=JL,IWD
      IZ=IZ+1
  940 W(IS)=W(IS)+W(IZ)*W(J)
      GOTO (350,550), KFLAG
  950 CONTINUE
!     ONE MORE CHECK THE DV BOUND
      TEMPA = 1.0D0
      DO I=1,N
        TEMP = X(I)+XVAR(I)
        IF (TEMP.LT.dvarl(I)) THEN
          TEMP = ABS((dvarl(I)-XVAR(I))/X(I))
          IF (TEMP.LT.1.0D-10) THEN
c           screen out numerical calculation error
            X(I) = dvarl(I)-XVAR(I)
          ELSE
            IF (TEMP.LT.TEMPA) TEMPA=TEMP
          ENDIF
        ELSE IF (TEMP.GT.dvaru(I)) THEN
          TEMP = ABS((dvaru(I)-XVAR(I))/X(I))
          IF (TEMP.LT.1.0D-10) THEN
c           screen out numerical calculation error
            X(I) = dvaru(I)-XVAR(I)
          ELSE
            IF (TEMP.LT.TEMPA) TEMPA=TEMP
          ENDIF
        ENDIF
      ENDDO
      IF (TEMPA.LT.1.0D0) THEN
        DO I=1,N
          X(I) = TEMPA*X(I)
        ENDDO
      ENDIF
      RETURN
      END
C
Chd|====================================================================
Chd|  F01CRF                        src/mater1/mincf.F            
Chd|-- called by -----------
Chd|-- calls ---------------
Chd|====================================================================
      SUBROUTINE F01CRF(A,M,N,MN,MOVE,IWRK,IFAIL)
      IMPLICIT NONE
      INTEGER M,N,MN,IWRK,IFAIL
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     *****
C     ALGORITHM 380 - REVISED
C     *****
C     A IS A ONE-DIMENSIONAL ARRAY OF LENGTH MN = M*N, WHICH
C     CONTAINS THE M X N MATRIX TO BE F01CRFPOSED
C     (STORED COLUMNWISE).
C     MOVE IS A ONE-DIMENSIONAL ARRAY OF LENGTH IWRK USED
C     TO STORE INFORMATION TO SPEED UP THE PROCESS. THE
C     VALUE IWRK = (M+N)/2 IS RECOMMENDED.
C     IFAIL INDICATES THE SUCCESS OR FAILURE OF THE ROUTINE.
C     NORMAL RETURN  IFAIL = 0
C     ERRORS         IFAIL = 1, MN NOT EQUAL TO M*N
C                    IFAIL = 2, IWRK NEGATIVE OR ZERO
C                    IFAIL.LT.ZERO, (SHOULD NEVER OCCUR),
C     IN THIS CASE WE SET IFAIL EQUAL TO THE FINAL VALUE OF I
C     WHEN THE SEARCH IS COMPLETED BUT SOME LOOPS HAVE NOT BEEN
C     MOVED NOTE * MOVE(I) WILL STAY ZERO FOR FIXED POINTS.
C     CHECK ARGUMENTS AND INITIALIZE.
C     .. Scalar Arguments ..
C     .. Array Arguments ..
      DOUBLE PRECISION  A(MN)
      INTEGER           MOVE(IWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, D
      INTEGER           I, I1, I1C, I2, I2C, IM, IR0, IR1, IR2, J, J1,
     *                  K, KMI, MAX, N1, NCOUNT
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      I1 = 0
      IF (M.LT.2 .OR. N.LT.2) GO TO 240
      IF (MN.NE.M*N) GO TO 360
      IF (IWRK.LT.1) GO TO 380
      IF (M == N) GO TO 260
      NCOUNT = 2
      K = MN - 1
      DO 20 I = 1, IWRK
         MOVE(I) = 0
   20 CONTINUE
      IF (M.LT.3 .OR. N.LT.3) GO TO 60
C     CALCULATE THE NUMBER OF FIXED POINTS, EUCLIDS ALGORITHM
C     FOR GCD(M-1,N-1).
      IR2 = M - 1
      IR1 = N - 1
   40 IR0 = MOD(IR2,IR1)
      IR2 = IR1
      IR1 = IR0
      IF (IR0.NE.0) GO TO 40
      NCOUNT = NCOUNT + IR2 - 1
C     SET INITIAL VALUES FOR SEARCH
   60 I = 1
      IM = M
C     AT LEAST ONE LOOP MUST BE RE-ARRANGED
      GO TO 160
C     SEARCH FOR LOOPS TO REARRANGE
   80 MAX = K - I
      I = I + 1
      IF (I.GT.MAX) GO TO 320
      IM = IM + M
      IF (IM.GT.K) IM = IM - K
      I2 = IM
      IF (I == I2) GO TO 80
      IF (I.GT.IWRK) GO TO 120
      IF (MOVE(I) == 0) GO TO 160
      GO TO 80
  100 I2 = M*I1 - K*(I1/N)
  120 IF (I2.LE.I .OR. I2.GE.MAX) GO TO 140
      I1 = I2
      GO TO 100
  140 IF (I2.NE.I) GO TO 80
C     REARRANGE THE ELEMENTS OF A LOOP AND ITS COMPANION LOOP
  160 I1 = I
      KMI = K - I
      B = A(I1+1)
      I1C = KMI
      C = A(I1C+1)
  180 I2 = M*I1 - K*(I1/N)
      I2C = K - I2
      IF (I1.LE.IWRK) MOVE(I1) = 2
      IF (I1C.LE.IWRK) MOVE(I1C) = 2
      NCOUNT = NCOUNT + 2
      IF (I2 == I) GO TO 220
      IF (I2 == KMI) GO TO 200
      A(I1+1) = A(I2+1)
      A(I1C+1) = A(I2C+1)
      I1 = I2
      I1C = I2C
      GO TO 180
C     FINAL STORE AND TEST FOR FINISHED
  200 D = B
      B = C
      C = D
  220 A(I1+1) = B
      A(I1C+1) = C
      IF (NCOUNT.LT.MN) GO TO 80
C     NORMAL RETURN
  240 IFAIL = 0
      RETURN
C     IF MATRIX IS SQUARE, EXCHANGE ELEMENTS A(I,J) AND A(J,I).
  260 N1 = N - 1
      DO 300 I = 1, N1
         J1 = I + 1
         DO 280 J = J1, N
            I1 = I + (J-1)*N
            I2 = J + (I-1)*M
            B = A(I1)
            A(I1) = A(I2)
            A(I2) = B
  280    CONTINUE
  300 CONTINUE
      GO TO 240
C     ERROR RETURNS.
  320 IFAIL = -I
  340 RETURN
  360 IFAIL = 1
      GO TO 340
  380 IFAIL = 2
      GO TO 340
      END

      
Chd|====================================================================
Chd|  ZEROD                         src/mater1/mincf.F            
Chd|-- called by -----------
Chd|        VF13BD                        src/mater1/mincf.F            
Chd|-- calls ---------------
Chd|====================================================================
      SUBROUTINE ZEROD(TAB,LEN)
        INTEGER LEN
        DOUBLE PRECISION TAB(LEN)
        TAB(1:LEN) = 0.D0
c-----------
      RETURN
      END 
     
      
Chd|====================================================================
Chd|  ZEROI                         src/mater1/mincf.F            
Chd|-- called by -----------
Chd|-- calls ---------------
Chd|====================================================================
      SUBROUTINE ZEROI(TAB,LEN)
        INTEGER LEN,TAB(LEN)
        TAB(1:LEN) = 0
c-----------
      RETURN
      END 

