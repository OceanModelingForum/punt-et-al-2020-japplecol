      INTEGER MAXYEAR,MAXSIM
      PARAMETER(MAXYEAR=1000,MAXSIM=10000)
C
C     Control parameters
      INTEGER Nsim,Nyear,ControlRule
      REAL*8 FinDepl(MAXSIM),Depl20(MAXSIM),AAV(MAXSIM),TEPBR(MAXSIM)
      CHARACTER*2 CASEA,CASEB,CASEC
      COMMON / Cntl / Nsim,Nyear,FinDepl,Depl20,ControlRule,AAV,TEPBR,
     +               CASEA,CASEB,CASEC 
C      Nsim                                        - Number of simulations
C      Nyear                                       - Number of years
C      FinDepl(s)                                  - Final depletion (year Nyear)
C      Depl20(s)                                   - Depletion after 20 years
C      ControlRule                                 - Which control rule to apply
C      AAV(s)                                      - AAV 
C      TEPBR(s)                                    - total limits
C
C     General biological parameters
      REAL*8 InitDepl,LambdaMax,Rpar,Carry,Zpar,MNPL
      COMMON / Pars / InitDepl,LambdaMax,Rpar,Carry,Zpar,MNPL
C      InitDepl                                   - Initial Depletion
C      LambdaMax                                  - Lamda(max)
C      rpar                                       - Rpar
C      Carry                                      - Carrying capacity
C      zpar                                       - Degree of compensation
C      MNPL                                       - MNPL
C
C     True population size
      REAL*8 N(-100:MaxYear),Cact(-100:MAXYEAR),CatchBias,CatchEstBias
      COMMON /TrueP / N,Cact,CatchBias,CatchEstBias
C      N(y)                                       - True population size in year y
C      Cact(y)                                    - True catch in year y
C      CatchBias                                  - Difference between PBR and expected catch
C      CatchEstBias                               - Bias in estimated catches  

C     Monitoring
      REAL*8 SurveyCV,Nhat(-100:MAXYEAR),ObsCV(-100:MAXYEAR)
      REAL*8 Nhat2(-100:MAXYEAR),ObsCV2(-100:MAXYEAR)
      REAL*8 SurveyCVEst,AbundBias
      REAL*8 CatchCV,CatchObsCV,Cest(-100:MAXYEAR)
      REAL*8 AveHistCat
      INTEGER SurveyFreq,UpdateFreq,FirstYrAbund,FirstYrCat
      COMMON / SurvStuff / SurveyCV,SurveyFreq,Nhat,ObsCV,CatchCV,
     +    CatchObsCV,Cest,AveHistCat,SurveyCVEst,AbundBias,
     +    UpdateFreq,FirstYrAbund,FirstYrCat,Nhat2,ObsCV2
C      SurveyCV                                   - Survey CV
C      SurveyFreq                                 - Survey frequency
C      Nhat(y)                                    - Estimated N
C      ObsCV(y)                                   - Observed CV
C      CatchCV                                    - Catch CV
C      CatchObsCV                                 - CV of the obseved catch about the true catch
C      Cest(y)                                    - Estimate of the catch
C      AveHistCat                                 - Average catch
C      SurveyCVEst                                - Esimate of the survey CV 
C      AbundBias                                  - Bias for the survey estimates
C      UpdateFreq                                 - Frequency of updates to PBR
C      FirstYrAbund                               - First year with abundance data
C      FirstYrCat                                 - First year with catches

C     Apply PBR
      REAL*8 NminMult,FrPBR,RmaxPBR
      COMMON / IsPBR / NminMult,FrPBR,RmaxPBR
C      MminMult                                   - Minimum multiplier
C      FrPBR                                      - Fr
C      RmaxPBR                                    - Assumed Rmax

C     Apply DCAC
      REAL*8 PercDCAC,MDCAC(50000),CDCAC(50000)
      INTEGER NsimDCAC
      COMMON / IsDCAC / PercDCAC,NsimDCAC,MDCAC,CDCAC
C      PercDCAC                                   - Percentile to use
C      NsimDCAC                                   - Number of sims
C      MDCAC(i)                                   - M values
C      CDCAC9i)                                   - c values
      
C     Apply RY
      REAL*8 PercRY
      INTEGER NsimRY
      COMMON / IsRY / PercRY,NsimRY
C      PercRY                                     - Percentile to use
C      NsimRY                                     - Number of sims
      
C     Apply Slope
      REAL*8 MultSlope
      COMMON / IsSlope / MultSlope
C      MultSlope                                  - Multiplier on the slope

C     Apply Perftec
      REAL*8 FRPERFECT
      COMMON / Isperfect / FRPERFECT
C      FRPERFECT                                  - FR  
