      SUBROUTINE MULT (TT,B,RR,MAXA,NN,NWM)                             MUT00001
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . MUT00002
C .                                                                   . MUT00003
C .   P R O G R A M                                                   . MUT00004
C .        TO EVALUATE PRODUCT OF B TIMES RR AND STORE RESULT IN TT   . MUT00005
C .                                                                   . MUT00006
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . MUT00007
C                                                                       MUT00008
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               MUT00009
      DIMENSION TT(1),B(1),RR(1),MAXA(1)                                MUT00010
C                                                                       MUT00011
      IF (NWM.GT.NN) GO TO 20                                           MUT00012
      DO 10 I=1,NN                                                      MUT00013
   10 TT(I)=B(I)*RR(I)                                                  MUT00014
      GO TO 900                                                         MUT00015
C                                                                       MUT00016
   20 DO 40 I=1,NN                                                      MUT00017
   40 TT(I)=0.                                                          MUT00018
      DO 100 I=1,NN                                                     MUT00019
      KL=MAXA(I)                                                        MUT00020
      KU=MAXA(I+1) - 1                                                  MUT00021
      II=I + 1                                                          MUT00022
      CC=RR(I)                                                          MUT00023
      DO 100 KK=KL,KU                                                   MUT00024
      II=II - 1                                                         MUT00025
  100 TT(II)=TT(II) + B(KK)*CC                                          MUT00026
      IF (NN.EQ.1) GO TO 900                                            MUT00027
      DO 200 I=2,NN                                                     MUT00028
      KL=MAXA(I) + 1                                                    MUT00029
      KU=MAXA(I+1) - 1                                                  MUT00030
      IF (KU-KL) 200,210,210                                            MUT00031
  210 II=I                                                              MUT00032
      AA=0.                                                             MUT00033
      DO 220 KK=KL,KU                                                   MUT00034
      II=II - 1                                                         MUT00035
  220 AA=AA + B(KK)*RR(II)                                              MUT00036
      TT(I)=TT(I) + AA                                                  MUT00037
  200 CONTINUE                                                          MUT00038
C                                                                       MUT00039
  900 RETURN                                                            MUT00040
      END                                                               MUT00041
