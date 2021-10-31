      SUBROUTINE DECOMP (A,MAXA,NN,ISH,IOUT)                            DEC00001
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . DEC00002
C .                                                                   . DEC00003
C .   P R O G R A M                                                   . DEC00004
C .        TO CALCULATE (L)*(D)*(L)(T) FACTORIZATION OF               . DEC00005
C .        STIFFNESS MATRIX                                           . DEC00006
C .                                                                   . DEC00007
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . DEC00008
C                                                                       DEC00009
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               DEC00010
      DIMENSION A(1),MAXA(1)                                            DEC00011
      IF (NN.EQ.1) GO TO 900                                            DEC00012
C                                                                       DEC00013
      DO 200 N=1,NN                                                     DEC00014
      KN=MAXA(N)                                                        DEC00015
      KL=KN + 1                                                         DEC00016
      KU=MAXA(N+1) - 1                                                  DEC00017
      KH=KU - KL                                                        DEC00018
      IF (KH) 304,240,210                                               DEC00019
  210 K=N - KH                                                          DEC00020
      IC=0                                                              DEC00021
      KLT=KU                                                            DEC00022
      DO 260 J=1,KH                                                     DEC00023
      IC=IC + 1                                                         DEC00024
      KLT=KLT - 1                                                       DEC00025
      KI=MAXA(K)                                                        DEC00026
      ND=MAXA(K+1) - KI - 1                                             DEC00027
      IF (ND) 260,260,270                                               DEC00028
  270 KK=MIN0(IC,ND)                                                    DEC00029
      C=0.                                                              DEC00030
      DO 280 L=1,KK                                                     DEC00031
  280 C=C + A(KI+L)*A(KLT+L)                                            DEC00032
      A(KLT)=A(KLT) - C                                                 DEC00033
  260 K=K + 1                                                           DEC00034
  240 K=N                                                               DEC00035
      B=0.                                                              DEC00036
      DO 300 KK=KL,KU                                                   DEC00037
      K=K - 1                                                           DEC00038
      KI=MAXA(K)                                                        DEC00039
      C=A(KK)/A(KI)                                                     DEC00040
      IF (ABS(C).LT.1.E07) GO TO 290                                    DEC00041
      WRITE (IOUT,2010) N,C                                             DEC00042
      GO TO 800                                                         DEC00043
  290 B=B + C*A(KK)                                                     DEC00044
  300 A(KK)=C                                                           DEC00045
      A(KN)=A(KN) - B                                                   DEC00046
  304 IF (A(KN)) 310,310,200                                            DEC00047
  310 IF (ISH.EQ.0) GO TO 320                                           DEC00048
      IF (A(KN).EQ.0.) A(KN)=-1.E-16                                    DEC00049
      GO TO 200                                                         DEC00050
  320 WRITE (IOUT,2000) N,A(KN)                                         DEC00051
      GO TO 800                                                         DEC00052
  200 CONTINUE                                                          DEC00053
      GO TO 900                                                         DEC00054
C                                                                       DEC00055
  800 STOP                                                              DEC00056
  900 RETURN                                                            DEC00057
C                                                                       DEC00058
 2000 FORMAT (//' STOP - STIFFNESS MATRIX NOT POSITIVE DEFINITE',//,    DEC00059
     1          ' NONPOSITIVE PIVOT FOR EQUATION ',I8,//,               DEC00060
     2          ' PIVOT = ',E20.12)                                     DEC00061
 2010 FORMAT (//' STOP - STURM SEQUENCE CHECK FAILED BECAUSE OF',       DEC00062
     1          ' MULTIPLIER GROWTH FOR COLUMN NUMBER ',I8,//,          DEC00063
     2          ' MULTIPLIER = ',E20.8)                                 DEC00064
      END                                                               DEC00065
