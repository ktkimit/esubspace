      SUBROUTINE REDBAK (A,V,MAXA,NN)                                   RED00001
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . RED00002
C .                                                                   . RED00003
C .   P R O G R A M                                                   . RED00004
C .        TO REDUCE AND BACK-SUBSTITUTE ITERATION VECTORS            . RED00005
C .                                                                   . RED00006
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . RED00007
C                                                                       RED00008
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               RED00009
      DIMENSION A(1),V(1),MAXA(1)                                       RED00010
C                                                                       RED00011
      DO 400 N=1,NN                                                     RED00012
      KL=MAXA(N) + 1                                                    RED00013
      KU=MAXA(N+1) - 1                                                  RED00014
      IF (KU-KL) 400,410,410                                            RED00015
  410 K=N                                                               RED00016
      C=0.                                                              RED00017
      DO 420 KK=KL,KU                                                   RED00018
      K=K - 1                                                           RED00019
  420 C=C + A(KK)*V(K)                                                  RED00020
      V(N)=V(N) - C                                                     RED00021
  400 CONTINUE                                                          RED00022
C                                                                       RED00023
      DO 480 N=1,NN                                                     RED00024
      K=MAXA(N)                                                         RED00025
  480 V(N)=V(N)/A(K)                                                    RED00026
      IF (NN.EQ.1) GO TO 900                                            RED00027
      N=NN                                                              RED00028
      DO 500 L=2,NN                                                     RED00029
      KL=MAXA(N) + 1                                                    RED00030
      KU=MAXA(N+1) - 1                                                  RED00031
      IF (KU-KL) 500,510,510                                            RED00032
  510 K=N                                                               RED00033
      DO 520 KK=KL,KU                                                   RED00034
      K=K - 1                                                           RED00035
  520 V(K)=V(K) - A(KK)*V(N)                                            RED00036
  500 N=N - 1                                                           RED00037
C                                                                       RED00038
  900 RETURN                                                            RED00039
      END                                                               RED00040
