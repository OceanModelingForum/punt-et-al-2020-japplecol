      IMPLICIT NONE
      INCLUDE "PROJECT4.INC"
C
C     Local variables
      INTEGER ISEED1,Isim,Iyear,II,RunType,Itype
      REAL*8 ProbRec100,ProbRec20,CatchMin,CatchMax,TheCatch,Product
      REAL*8 PercMin,PercMax,PercAve,Target1,Target2,Target
      REAL*8 ProbFinA,ProbFinB,Lower5A,Lower5B,Pred,ProbFin,Lower5
C
C     Read general specs
      OPEN(UNIT=13,FILE="Input.dat")
      READ(13,*) RunType
      READ(13,*) Target1
      READ(13,*) Target2
      CLOSE(13)
C
C     Read in the inputs
      CALL ReadIN()
C
C     Create an Output File
      OPEN(UNIT=15,FILE="DUMP.OUT")
      OPEN(UNIT=16,FILE="FINAL.OUT",ACCESS="APPEND")
C
C     Solve historical catch
      CatchMin = 0
      CatchMax = Carry*0.1
      N(-50) = Carry
      DO 90000 II = 1,30
       TheCatch = (CatchMin+CatchMax)/2.0
       DO 91000 Iyear = -50,-1
        Cact(Iyear) = TheCatch
        Product = Rpar*N(Iyear)*(1.0-(N(Iyear)/Carry)**Zpar)
        N(Iyear+1) = N(Iyear) + Product - Cact(Iyear)
        IF (N(Iyear+1).LT.0) N(Iyear+1) = 0
91000  CONTINUE        
       IF (N(0).LT.InitDepl*Carry) THEN
        CatchMax = TheCatch
       ELSE
        CatchMin = TheCatch
       ENDIF
90000 CONTINUE       

      IF (RunType.EQ.1) THEN
       DO 10000 Itype = 1,2
        IF (Itype.EQ.1) Target = Target1
        IF (Itype.EQ.2) Target = Target2
  
        PercMin = 0
        if (ControlRule.EQ.1) PercMax = 1
        if (ControlRule.EQ.2) PercMax = 1
        if (ControlRule.EQ.3) PercMax = 2
        if (ControlRule.EQ.4) PercMax = 1
C
C       Do bisection
        DO 11000 II = 1,20
C       
         PercAve = (PercMin+PercMax)/2.0
         if (ControlRule.EQ.1) PercDCAC = PercAve
         if (ControlRule.EQ.2) PercRY = PercAve
         if (ControlRule.EQ.3) MultSlope = PercAve
         if (ControlRule.EQ.4) FRPerfect = PercAve
C        
         SurveyCV = 0.2
         SurveyCVEst = SurveyCV
         ProbRec100 = 0  
         ISEED1 = -89191
         DO 11100 Isim = 1,Nsim
          CALL Project1(ISEED1,Isim,0)
          IF (FinDepl(Isim).GT.MNPL) ProbRec100 = ProbRec100 + 1
11100    CONTINUE       
         CALL SORT(FinDepl,Nsim)
         ProbFinA = ProbRec100/FLOAT(NSIM)
         Lower5A = FinDepl(Nsim*0.05)
       
         SurveyCV = 0.8
         SurveyCVEst = SurveyCV
         ProbRec100 = 0  
         ISEED1 = -89191
         DO 11200 Isim = 1,Nsim
          CALL Project1(ISEED1,Isim,0)
          IF (FinDepl(Isim).GT.MNPL) ProbRec100 = ProbRec100 + 1
11200    CONTINUE       
         CALL SORT(FinDepl,Nsim)
         ProbFinB = ProbRec100/FLOAT(NSIM)
         Lower5B = FinDepl(Nsim*0.05)

         ProbFin = min(ProbFinA,ProbFinB)
         Lower5 = min(Lower5A,Lower5B)
         WRITE(*,700) PercAve,ProbFinA,ProbFinB,Lower5A,Lower5B
700      FORMAT(1x,5(F5.3,1x))
       
         IF (Itype.EQ.1) Pred = ProbFin
         IF (Itype.EQ.2) Pred = Lower5
         WRITE(*,701) PercAve,Target,Pred
701      FORMAT(1x,5(F7.5,1x))
        
C         IF (ABS(Pred-Target).LT.0.001) GOTO 10999
         IF (ControlRule.EQ.1.OR.ControlRule.EQ.2.OR.
     +       ControlRule.EQ.4) THEN
          IF (Pred.LT.Target) THEN
           PercMax = PercAve
          ELSE 
           PercMin= PercAve
          ENDIF 
         ENDIF
         IF (ControlRule.EQ.3) THEN
          IF (Pred.LT.Target) THEN
           PercMin = PercAve
          ELSE 
           PercMax= PercAve
          ENDIF 
         ENDIF
11000   CONTINUE       
10999   CONTINUE
        WRITE(*,702) Itype,PercAve,ProbFin,Lower5
702     FORMAT(1x,I2,1x,5(F7.5,1x))       
10000  CONTINUE      
C      
      ENDIF
C
      IF (RunType.EQ.0) THEN

       ProbRec100 = 0  
       ProbRec20 = 0  
       ISEED1 = -89191
       DO 20000 Isim = 1,Nsim
C       WRITE(*,*) Isim
        CALL Project1(ISEED1,Isim,1)
        IF (FinDepl(Isim).GT.MNPL) ProbRec100 = ProbRec100 + 1
        IF (Depl20(Isim).GT.MNPL) ProbRec20 = ProbRec20 + 1
20000  CONTINUE       
       CALL SORT(FinDepl,Nsim)
       CALL SORT(Depl20,Nsim)
       CALL SORT(AAV,Nsim)
       CALL SORT(TEPBR,Nsim)
       WRITE(16,600) CASEA,CASEB,CASEC,ControlRule, 
     +  InitDepl,
     +  AbundBias,SurveyFreq,SurveyCV,SurveyCVEst,
     +  LambdaMax,MNPL,
     +  CatchBias,CatchEstBias,CatchCV,CatchObsCV,
     +  NMINMult,PercDCAC,PercRY,MultSlope,FRPerfect,
     +  ProbRec100/FLOAT(NSIM),ProbRec20/FLOAT(NSIM),
     +  FinDepl(Nsim*0.05),Depl20(Nsim*0.05),
     +  MIN(10000.0d0,TEPBR(Nsim*0.5)),AAV(Nsim*0.5)
C      
       CLOSE(15)
       CLOSE(16)
C       
      ENDIF
C
      STOP
600   FORMAT(1x,3(A2,1x),I2,1x,2(F6.3,1x),I4,1x,8(F6.3,1x),
     +       5(1x,F6.4),4(1x,F5.3),1x,F10.3,1x,F5.3)
      END
C
C =================================================================================================================
C
      SUBROUTINE ReadIN()
C
C This subroutine reads in the basic data
C
      IMPLICIT NONE
      INCLUDE "PROJECT4.INC"
C
C     Local variables
      INTEGER Isim,ISEED
      REAL*8 Error,LowM,HiM
      REAL*8 RAN2,XNORM
      EXTERNAL RAN2,XNORM
C
C     Set parameter
      Carry = 1000.0
      
      OPEN(UNIT=13,FILE="OPMODEL.DAT")
      READ(13,'(40X,I4)') Nsim
      READ(13,'(40X,I4)') Nyear
      READ(13,'(40X,A2)') CASEA
      READ(13,'(40X,A2)') CASEB
      READ(13,'(40X,A2)') CASEC
      READ(13,'(40X,I4)') FirstYrCat
      READ(13,'(40X,I4)') FirstYrAbund
           
      UpdateFreq = 4
     
      InitDepl = 0.3
      SurveyCV = 0.2
      SurveyCVEst = SurveyCV
      MNPL = 0.5
      Zpar = 1.0
      SurveyFreq = 4
      CatchCV = 0.3
      CAtchEstBias = 1
      
      CatchObsCV = 0.3
      CatchBias = 1.0
      AbundBias = 1.0
      IF (CASEA.EQ.'Ce') LambdaMax = 1.04
      IF (CASEA.EQ.'Pi') LambdaMax = 1.12
      IF (CASEB.EQ.'08') THEN
       SurveyCV = 0.8
       SurveyCVEst = 0.8
      ENDIF 
      IF (CASEC.EQ.'01') CatchBias = 2.0
      IF (CASEC.EQ.'02') AbundBias = 2.0
      IF (CASEC.EQ.'03') LambdaMax = 1.0+(LambdaMax-1.0) / 2.0
      IF (CASEC.EQ.'04'.AND.SurveyCV.EQ.0.8) SurveyCV = 1.6
      IF (CASEC.EQ.'04'.AND.SurveyCV.EQ.0.2) SurveyCV = 0.8
      IF (CASEC.EQ.'05') CatchCV =  1.20
      IF (CASEC.EQ.'06') SurveyFreq = 8
      IF (CASEC.EQ.'07') THEN
       MNPL = 0.45
       Zpar = 0.53
      ENDIF
      IF (CASEC.EQ.'08') THEN
       MNPL = 0.70
       Zpar = 5.04
       CatchBias = 2.0
      ENDIF
      IF (CASEC.EQ.'09') InitDepl = 0.5
      
      IF (CASEC.EQ.'11') AbundBias = 0.1
      IF (CASEC.EQ.'12') AbundBias = 0.5
      IF (CASEC.EQ.'13') AbundBias = 1.0
      IF (CASEC.EQ.'14') AbundBias = 2.0
      IF (CASEC.EQ.'15') AbundBias = 4.0
      IF (CASEC.EQ.'16') AbundBias = 6.0
      IF (CASEC.EQ.'17') AbundBias = 8.0
      IF (CASEC.EQ.'18') AbundBias = 10.0
      
      IF (CASEC.EQ.'20') SurveyFreq = 1
      IF (CASEC.EQ.'21') SurveyFreq = 2
      IF (CASEC.EQ.'22') SurveyFreq = 3
      IF (CASEC.EQ.'23') SurveyFreq = 5
      IF (CASEC.EQ.'24') SurveyFreq = 6
      IF (CASEC.EQ.'25') SurveyFreq = 7
      IF (CASEC.EQ.'26') SurveyFreq = 9
      
      IF (CASEC.EQ.'31') CatchEstBias = 0.1
      IF (CASEC.EQ.'32') CatchEstBias = 0.5
      IF (CASEC.EQ.'33') CatchEstBias = 1.0
      IF (CASEC.EQ.'34') CatchEstBias = 2.0
      IF (CASEC.EQ.'35') CatchEstBias = 4.0
      IF (CASEC.EQ.'36') CatchEstBias = 6.0
      IF (CASEC.EQ.'37') CatchEstBias = 8.0
      IF (CASEC.EQ.'38') CatchEstBias = 10.0

      IF (CASEC.EQ.'41') CatchObsCV = 0.1
      IF (CASEC.EQ.'42') CatchObsCV = 0.3
      IF (CASEC.EQ.'43') CatchObsCV = 0.5
      IF (CASEC.EQ.'44') CatchObsCV = 0.7
      IF (CASEC.EQ.'45') CatchObsCV = 0.9
      IF (CASEC.EQ.'46') CatchObsCV = 1.1
      IF (CASEC.EQ.'47') CatchObsCV = 1.3
      IF (CASEC.EQ.'48') CatchObsCV = 1.5
      IF (CASEC.EQ.'49') CatchObsCV = 1.7
      
      Rpar = LambdaMax - 1.0
      CLOSE(13)
      
      OPEN(UNIT=13,FILE="CONTROLRULE.DAT")
      READ(13,'(40X,I4)') ControlRule
C
C     PBR specifications
      READ(13,*)
      READ(13,'(40X,F10.0)') NMINMult
      READ(13,'(40X,F10.0)') FrPBR
      READ(13,'(40X,F10.0)') RmaxPBR
C
C     DCAC specifications
      READ(13,*)
      READ(13,'(40X,F10.0)') PercDCAC
      READ(13,'(40X,I5)') NsimDCAC
      READ(13,'(40X,F10.0)') LowM
      READ(13,'(40X,F10.0)') HiM
C
C     RY specifications
      READ(13,*)
      READ(13,'(40X,F10.0)') PercRY
      READ(13,'(40X,I5)') NsimRY
C
C     Slope specifications
      READ(13,*)
      READ(13,'(40X,F10.0)') MultSlope
C
C     Perfect specifications
      READ(13,*)
      READ(13,'(40X,F10.0)') FRPerfect
C      
      WRITE(*,*)
      
      CLOSE(13)
C
      ISEED = -10101
      DO 10000 Isim = 1,NsimDCAC
       MDCAC(Isim) = LowM+RAN2(ISEED)*(HiM-LowM)
       Error = XNORM(2,0.0d0,0.2d0,ISEED)
       cDCAC(Isim) = EXP(Error - 0.2*0.2/2.0)
10000 CONTINUE
C
      RETURN
      END
C
C =================================================================================================================
C
      Subroutine Project1(ISEED1,Isim,DoDump)
C
C This conducts a single projection
C
      IMPLICIT NONE
      INCLUDE "PROJECT4.INC"
C
C     Global variables
      INTEGER ISEED1,Isim,DoDump
C
C     Local variables
      REAL*8 Product,Error,EPBR,Sigma,TOTALC,TOTALEPBR
      REAL*8 TheLimit(-100:MAXYEAR)
      INTEGER Iyear
      REAL*8 XNORM,CalcPBR,CalcDCAC,CalcRY,CalcSlope
      EXTERNAL XNORM,CalcPBR,CalcDCAC,CalcRY,CalcSlope
C
C     Reset observed data
      Nhat = -1
      ObsCV = -1
      Cest = -1
      DO 2000 Iyear = -50,-1
        Error = XNORM(1,0.0d0,1.0d0,ISEED1)
        Sigma = SQRT(log(1.0+CatchObsCV**2.0))
        Cest(Iyear) = Cact(Iyear)*CATCHESTBIAS*
     +     EXP(Error*Sigma-Sigma**2.0/2.0)
        IF (MOD(Iyear,SurveyFreq).EQ.0) THEN
         Error = XNORM(1,0.0d0,1.0d0,ISEED1)
         Sigma = SQRT(log(1.0+SurveyCV**2.0))
         Nhat(Iyear) = AbundBias*N(Iyear)*
     +    EXP(Error*Sigma-Sigma**2.0/2.0)
         ObsCV(Iyear) = SurveyCVEst
        ENDIF
2000  CONTINUE      
C
C     FInd the average catch over year -1 to -20
      AveHistCat = 0
      DO 2100 Iyear = -20,-1
       AveHistCat = AveHistCat + Cest(Iyear)
2100  CONTINUE
      AveHistCat = AveHistCat / 20.0
C
C     Initial depletion
      TOTALEPBR = 0
      DO 1000 Iyear = 0,Nyear
C
C      Generate the data
       IF (MOD(Iyear,SurveyFreq).EQ.0) THEN
        Error = XNORM(1,0.0d0,1.0d0,ISEED1)
        Sigma = SQRT(log(1.0+SurveyCV**2.0))
        Nhat(Iyear) = AbundBias*N(Iyear)*
     +   EXP(Error*Sigma-Sigma**2.0/2.0)
        ObsCV(Iyear) = SurveyCVEst
       ENDIF
C
C      Apply the control rule
       IF (ControlRule.EQ.0.AND.MOD(Iyear,UpdateFreq).EQ.0) 
     +     EPBR = CalcPBR(Iyear)
       IF (ControlRule.EQ.1.AND.MOD(Iyear,UpdateFreq).EQ.0) 
     +     EPBR = CalcDCAC(Iyear)
       IF (ControlRule.EQ.2.AND.MOD(Iyear,UpdateFreq).EQ.0) 
     +     EPBR = CalcRY(Iyear)
       IF (ControlRule.EQ.3.AND.MOD(Iyear,UpdateFreq).EQ.0) 
     +     EPBR = CalcSlope(Iyear)
       IF (ControlRule.EQ.4.AND.MOD(Iyear,UpdateFreq).EQ.0) 
     +     EPBR = FRPerfect*0.5*Rpar*N(Iyear)
       IF (EPBR.LT.0) EPBR = 0
       IF (Iyear.LT.100) TOTALEPBR = TOTALEPBR + EPBR
       TheLimit(Iyear)= EPBR
       Error = XNORM(1,0.0d0,1.0d0,ISEED1)
       Cact(Iyear) = CatchBias*EPBR*(1.0 + Error*CatchCV)
       IF (Cact(Iyear).LT.0) Cact(Iyear) = 0
C
C      Generate the catches
       Error = XNORM(1,0.0d0,1.0d0,ISEED1)
       Sigma = SQRT(log(1.0+CatchObsCV**2.0))
       Cest(Iyear) = Cact(Iyear)*CATCHESTBIAS*
     +  EXP(Error*Sigma-Sigma**2.0/2.0)
      
       Product = Rpar*N(Iyear)*(1.0-(N(Iyear)/Carry)**Zpar)
       N(Iyear+1) = N(Iyear) + Product - Cact(Iyear)
       IF (N(Iyear+1).LT.0) N(Iyear+1) = 0
       IF (DoDump.EQ.1.AND.Isim.LE.10000) WRITE(15,600) Isim,Iyear,
     +  N(Iyear),Nhat(Iyear),EPBR,Cact(Iyear),Cest(Iyear)
1000  CONTINUE      
      FinDepl(Isim) = N(Nyear)/Carry
      Depl20(Isim) = N(20)/CARRY
      AAV(Isim) = 0
      TOTALC = 0
      DO 3000 Iyear = 1,Nyear-1
       AAV(Isim) = AAV(Isim) + ABS(TheLimit(Iyear)-TheLimit(Iyear-1))
       TOTALC = TOTALC + TheLimit(Iyear-1)
3000  CONTINUE      
      AAV(Isim) = AAV(Isim)/(TOTALC+0.00001)
      TEPBR(Isim) = TOTALEPBR
C      
C      STOP
C      
      RETURN
600   FORMAT(1x,I4,1x,I4,2(1x,F12.4),3(1x,F9.4))
      END
C
C ==========================================================================================================================
C
      REAL*8 FUNCTION CalcPBR(TheYear)
C
C This subroutine computes a PBR based on recent abundance data
C
      IMPLICIT NONE
      INCLUDE "PROJECT4.INC"
C
C     Global variables
      INTEGER TheYear
C
C     Local variables
      INTEGER Year,LastYear
      REAL*8 Nmin,Term1,AbundVal,AbundCV,NetCV
C
      LastYear = -1000
      DO 51000 Year = 0,TheYear
        IF (ObsCV(Year).GT.0) LastYear = Year
51000 CONTINUE       
      AbundVal = NHat(LastYear)
      AbundCV = ObsCV(LastYear)
C
C     Compute Nmin
      NetCV = AbundCV
      Term1 = NminMult*SQRT(Log(1.0+NetCV**2.0))
      Nmin = AbundVal/EXP(Term1)
C
C     Multiply by Fmsy
      CalcPBR = FrPBR*Nmin*0.5*RmaxPBR
C      
      RETURN
      END
C
C ==========================================================================================================================
C
      REAL*8 FUNCTION CalcDCAC(TheYear)
C
C This subroutine computes a PBR based on recent abundance data
C
      IMPLICIT NONE
      INCLUDE "PROJECT4.INC"
C
C     Global variables
      INTEGER TheYear
C
C     Local variables
      INTEGER ISEED,ISIM,Iyear,Nfound
      REAL*8 TotCat,SumX,SumY,SumXX,SumXY,Slope,Change,EChange
      REAL*8 Error,Ncatch,DCAC(21000)
      REAL*8 XBAR,Intercept,SS,Denom,Pred,SigmaSlope,VarChange,CVSlope
      REAL*8 RAN2,XNORM
      EXTERNAL RAN2,XNORM
C
C     Find the catch data to use
      TotCat = 0
      Ncatch = 0
      DO 10000 Iyear = TheYear-1,MAX(TheYear-20,FirstYrCat),-1
       TotCat = TotCat + Cest(Iyear)
       Ncatch = Ncatch + 1.0
10000 CONTINUE
C
C     Find the abundance 
      Nfound = 0
      SumX = 0
      SumXX = 0
      SumY = 0
      SumXY = 0
      DO 20000 Iyear = TheYear,MAX(-50,FirstYrAbund),-1
       IF (ObsCV(Iyear).GT.0.AND.
     +   (Nfound.LT.5.OR.TheYear-Iyear.LT.20)) THEN
        Nfound = Nfound + 1
        SumX = SumX + Iyear
        SumXY = SumXY + log(0.1d0+NHat(Iyear))*FLOAT(Iyear)
        SumY = SumY + log(0.1d0+NHat(Iyear))
        SumXX = SumXX + Iyear*Iyear
       ENDIF
20000 CONTINUE
      Slope = (SumXY-SumX*SumY/FLOAT(Nfound))/
     +        (SumXX-SumX*SumX/FLOAT(Nfound))
      XBAR = SumX/FLOAT(Nfound)
      Intercept = SumY/FLOAT(Nfound) - Slope*XBAR
      EChange = 1.0-EXP(Slope*16)
C
C     Sigma of beta
      SS = 0
      Denom = 0
      Nfound = 0
      DO 21000 Iyear = TheYear,MAX(-50,FirstYrAbund),-1
       IF (ObsCV(Iyear).GT.0.AND.
     +   (Nfound.LT.5.OR.TheYear-Iyear.LT.20)) THEN
        Nfound = Nfound + 1
        Pred = Intercept + Slope*Iyear
        SS = SS + (Pred-log(0.1d0+NHat(Iyear)))**2.0
        Denom = Denom + (Iyear-XBAR)**2.0
       ENDIF
21000 CONTINUE
      IF (Nfound.LE.2) THEN
       SigmaSlope = 0.2
      ELSE
       SigmaSlope = SQRT(SS/(FLOAT(Nfound-2))/Denom)
      ENDIF
      VarChange = (16*EXP(Slope*16))**2.0*SigmaSlope**2.0
C      WRITE(*,*) "E",EChange
      CVSlope = ABS(SQRT(VarChange)/ABS(EChange+0.0001))
C      WRITE(*,*) "C",CVSLOPE
C
C     Now simulate
      ISEED = -190101
      DO 30000 Isim = 1,NsimDCAC
       Error = XNORM(2,0.0d0,CVSLOPE,ISEED)
       Change = EChange*EXP(Error - CVSLOPE*CVSLOPE/2.0)
       IF (Change.GT.0) Change = 1
       DCAC(Isim) = TotCat/
     +             (Ncatch + Change/(0.5*cDCAC(Isim)*MDCAC(Isim)))
30000 CONTINUE      
      CALL SORT(DCAC,NsimDCAC)  
      CalcDCAC = DCAC(PercDCAC*NsimDCAC)
C      WRITE(*,*) CalcDCAC
C      
      RETURN
      END
C
C ==========================================================================================================================
C
      REAL*8 FUNCTION CalcRY(TheYear)
C
C This subroutine computes a PBR based on the replacement yield method
      IMPLICIT NONE
      INCLUDE "PROJECT4.INC"
C
C     Global variables
      INTEGER TheYear
C
C     Local variables
      REAL*8 RY(21000),LOGL(21000),B(-100:MAXYEAR)
      REAL*8 q,logq,Npnts,Sigma,Residual,RYV,LogLV,TotalL,SumL
      INTEGER Isim,ISEED,Iyear,Yr1
      REAL*8 RAN2
      EXTERNAL RAN2
C
C     Now simulate
      ISEED = -190101
      TOTALL = 0
      Yr1 = MAX(FirstYrAbund,FirstYrCat,TheYear-20)
      DO 30000 Isim = 1,NsimRY
       B(Yr1) = AveHistCat+(2000.0-AveHistCat)*RAN2(ISEED)
       RYV = AveHistCat*2*RAN2(ISEED)
       Logq = 0
       Npnts = 0
       DO 31000 Iyear = Yr1,TheYear
        B(Iyear+1) = B(Iyear) + RYV - Cest(Iyear)
        IF (B(Iyear+1).LT.1) B(Iyear+1) = 1
        IF (ObsCV(Iyear).GT.0) THEN
         Npnts = Npnts + 1
         logq = logq + log(0.1d0+NHat(Iyear)) - log(B(Iyear))
        ENDIF
31000  CONTINUE
       q = exp(logq/Npnts)
       Sigma = 0
       DO 32000 Iyear = MAX(FirstYrAbund,FirstYrCat,TheYear-20),TheYear
        IF (ObsCV(Iyear).GT.0) THEN
         Residual = log(0.1d0+NHat(Iyear)) - log(q*B(Iyear))
         Sigma = Sigma + Residual*Residual
        ENDIF
32000  CONTINUE
       Sigma = SQRT(Sigma/Npnts)
       LogLV = Npnts*log(Sigma)
       RY(Isim) = RYV
       LogL(Isim) = EXP(-LogLV)
C       WRITE(*,*) RY(Isim),LogL(Isim),Sigma
       TotalL = TotalL + logL(Isim)
30000 CONTINUE      
C
      CALL SORT2(RY,LogL,NsimRY)
      SumL = 0
      DO 40000 Isim = 1,NsimRY
       SumL = SumL + LogL(Isim)/TotalL
C       WRITE(*,*) RY(Isim),LogL(Isim),SumL
       IF (SumL.GT.PercRY) THEN
        CalcRY = RY(Isim)
        GOTO 40001
       ENDIF 
40000 CONTINUE
40001 CONTINUE
C      WRITE(*,*) CalcRY
C
      RETURN
      END
C
C ==========================================================================================================================
C
      REAL*8 FUNCTION CalcSlope(TheYear)
C
C This subroutine computes a PBR based on the replacement yield method
      IMPLICIT NONE
      INCLUDE "PROJECT4.INC"
C
C     Global variables
      INTEGER TheYear
C
C     Local variables
      INTEGER Nfound,Iyear
      REAL*8 SUMX,SUMY,SUMXX,SUMXY,Slope,XBAR,Intercept
      REAL*8 SS,Denom,SigmaSlope,TotCat,Ncatch,AveCat,Pred
C
C     Find the catch data to use
      TotCat = 0
      Ncatch = 0
      DO 5000 Iyear = TheYear-1,MAX(TheYear-20,FirstYrCat),-1
       TotCat = TotCat + Cest(Iyear)
       Ncatch = Ncatch + 1.0
5000  CONTINUE
      Avecat = Totcat / Ncatch
C
C     Find the abundance 
      Nfound = 0
      SumX = 0
      SumXX = 0
      SumY = 0
      SumXY = 0
      DO 10000 Iyear = TheYear,MAX(-50,FirstYrAbund),-1
       IF (ObsCV(Iyear).GT.0.AND.
     +   (Nfound.LT.5.OR.TheYear-Iyear.LT.20)) THEN
        Nfound = Nfound + 1
        SumX = SumX + Iyear
        SumXY = SumXY + log(0.1d0+NHat(Iyear))*FLOAT(Iyear)
        SumY = SumY + log(0.1d0+NHat(Iyear))
        SumXX = SumXX + Iyear*Iyear
       ENDIF
10000 CONTINUE
      Slope = (SumXY-SumX*SumY/FLOAT(Nfound))/
     +        (SumXX-SumX*SumX/FLOAT(Nfound))
      XBAR = SumX/FLOAT(Nfound)
      Intercept = SumY/FLOAT(Nfound) - Slope*XBAR
C
C     Sigma of beta
      SS = 0
      Denom = 0
      Nfound = 0
      DO 20000 Iyear = TheYear,-50,-1
       IF (ObsCV(Iyear).GT.0.AND.
     +  (Nfound.LT.5.OR.TheYear-Iyear.LT.20)) THEN
        Nfound = Nfound + 1
        Pred = Intercept + Slope*Iyear
        SS = SS + (Pred-log(0.1d0+NHat(Iyear)))**2.0
        Denom = Denom + (Iyear-XBAR)**2.0
       ENDIF
20000 CONTINUE
      IF (Nfound.LE.2) THEN
       SigmaSlope = 0.2
      ELSE
       SigmaSlope = SQRT(SS/(FLOAT(Nfound-2))/Denom)
      ENDIF
C
C     Update the catch limit
      CalcSlope = Avecat*(1.0+20*Slope-20*SigmaSlope*MultSlope)
C
      RETURN
      END
C
C =================================================================================================================
C
      INCLUDE "COMMON.FOR"
      