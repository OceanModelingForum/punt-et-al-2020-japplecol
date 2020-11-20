      IMPLICIT NONE      
      REAL*8 POP(2000,0:100),PBR(2000,0:100),R1
      REAL*8 RES(20),ERROR(20000)
      INTEGER Isim,I1,I2,Iyear,II
      
      OPEN(UNIT=13,FILE="DUMP.OUT")
      DO 10000 Isim = 1,2000
       DO 10000 Iyear = 0,100
        READ(13,*) I1,I2,POP(Isim,Iyear),R1,PBR(Isim,Iyear) 
10000 CONTINUE
      CLOSE(13)

      OPEN(UNIT=14,FILE="DUMP2.OUT")
      DO 20000 Iyear = 0,100
       DO 21000 II = 1,2000
        Error(II) = POP(II,Iyear)
21000  CONTINUE
       CALL SORT(Error,2000)
       RES(1) = Error(2000*0.05)
       RES(2) = Error(2000*0.25)
       RES(3) = Error(2000*0.5)
       RES(4) = Error(2000*0.75)
       RES(5) = Error(2000*0.95)
       RES(6) = POP(1,Iyear)
       RES(7) = POP(2,Iyear)
       RES(8) = POP(3,Iyear)
       RES(9) = POP(4,Iyear)
       RES(10) = POP(5,Iyear)
       DO 22000 II = 1,2000
        Error(II) = PBR(II,Iyear)
22000  CONTINUE
       CALL SORT(Error,2000)
       RES(11) = Error(2000*0.05)
       RES(12) = Error(2000*0.25)
       RES(13) = Error(2000*0.5)
       RES(14) = Error(2000*0.75)
       RES(15) = Error(2000*0.95)
       RES(16) = PBR(1,Iyear)
       RES(17) = PBR(2,Iyear)
       RES(18) = PBR(3,Iyear)
       RES(19) = PBR(4,Iyear)
       RES(20) = PBR(5,Iyear)
       WRITE(*,600) Iyear,(Res(II),II=1,20)
20000 CONTINUE      
      CLOSE(14)
      
      STOP
600   FORMAT(1x,I4,1x,20(F10.4,1x))      
      END
C
C ===========================================================================================================
C
      INCLUDE "COMMON.FOR"
      
      