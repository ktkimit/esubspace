      SUBROUTINE JACOBI (A,B,X,EIGV,D,N,NWA,RTOL,NSMAX,IFPR,IOUT)       JAC00001
C ..................................................................... JAC00002
C .                                                                   . JAC00003
C .   P R O G R A M                                                   . JAC00004
C .        TO SOLVE THE GENERALIZED EIGENPROBLEM USING THE            . JAC00005
C .        GENERALIZED JACOBI ITERATION                               . JAC00006
C ..................................................................... JAC00007
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               JAC00008
      DIMENSION A(NWA),B(NWA),X(N,N),EIGV(N),D(N)                       JAC00009
C                                                                       JAC00010
C     INITIALIZE EIGENVALUE AND EIGENVECTOR MATRICES                    JAC00011
C                                                                       JAC00012
      N1=N + 1                                                          JAC00013
      II=1                                                              JAC00014
      DO 10 I=1,N                                                       JAC00015
      IF (A(II).GT.0. .AND. B(II).GT.0.) GO TO 4                        JAC00016
      WRITE (IOUT,2020) II,A(II),B(II)                                  JAC00017
      GO TO 800                                                         JAC00018
    4 D(I)=A(II)/B(II)                                                  JAC00019
      EIGV(I)=D(I)                                                      JAC00020
   10 II=II + N1 - I                                                    JAC00021
      DO 30 I=1,N                                                       JAC00022
      DO 20 J=1,N                                                       JAC00023
   20 X(I,J)=0.                                                         JAC00024
   30 X(I,I)=1.                                                         JAC00025
      IF (N.EQ.1) GO TO 900                                             JAC00026
C                                                                       JAC00027
C     INITIALIZE SWEEP COUNTER AND BEGIN ITERATION                      JAC00028
C                                                                       JAC00029
      NSWEEP=0                                                          JAC00030
      NR=N - 1                                                          JAC00031
   40 NSWEEP=NSWEEP + 1                                                 JAC00032
      IF (IFPR.EQ.1) WRITE (IOUT,2000) NSWEEP                           JAC00033
C                                                                       JAC00034
C     CHECK IF PRESENT OFF-DIAGONAL ELEMENT IS LARGE ENOUGH TO REQUIRE  JAC00035
C     ZEROING                                                           JAC00036
C                                                                       JAC00037
      EPS=(.01)**(NSWEEP*2)                                             JAC00038
      DO 210 J=1,NR                                                     JAC00039
      JP1=J + 1                                                         JAC00040
      JM1=J - 1                                                         JAC00041
      LJK=JM1*N - JM1*J/2                                               JAC00042
      JJ=LJK + J                                                        JAC00043
      DO 210 K=JP1,N                                                    JAC00044
      KP1=K + 1                                                         JAC00045
      KM1=K - 1                                                         JAC00046
      JK=LJK + K                                                        JAC00047
      KK=KM1*N - KM1*K/2 + K                                            JAC00048
      EPTOLA=(A(JK)/A(JJ))*(A(JK)/A(KK))                                JAC00049
      EPTOLB=(B(JK)/B(JJ))*(B(JK)/B(KK))                                JAC00050
      IF (EPTOLA.LT.EPS .AND. EPTOLB.LT.EPS) GO TO 210                  JAC00051
C                                                                       JAC00052
C     IF ZEROING IS REQUIRED, CALCULATE THE ROTATION MATRIX ELEMENTS CA JAC00053
C     AND CG                                                            JAC00054
C                                                                       JAC00055
      AKK=A(KK)*B(JK) - B(KK)*A(JK)                                     JAC00056
      AJJ=A(JJ)*B(JK) - B(JJ)*A(JK)                                     JAC00057
      AB=A(JJ)*B(KK) - A(KK)*B(JJ)                                      JAC00058
      SCALE=A(KK)*B(KK)                                                 JAC00059
      ABCH=AB/SCALE                                                     JAC00060
      AKKCH=AKK/SCALE                                                   JAC00061
      AJJCH=AJJ/SCALE                                                   JAC00062
      CHECK=(ABCH*ABCH+4.0*AKKCH*AJJCH)/4.0                             JAC00063
      IF (CHECK) 50,60,60                                               JAC00064
   50 WRITE (IOUT,2020) JJ,A(JJ),B(JJ)                                  JAC00065
      GO TO 800                                                         JAC00066
   60 SQCH=SCALE*SQRT(CHECK)                                            JAC00067
      D1=AB/2. + SQCH                                                   JAC00068
      D2=AB/2. - SQCH                                                   JAC00069
      DEN=D1                                                            JAC00070
      IF (ABS(D2).GT.ABS(D1)) DEN=D2                                    JAC00071
      IF (DEN) 80,70,80                                                 JAC00072
   70 CA=0.                                                             JAC00073
      CG=-A(JK)/A(KK)                                                   JAC00074
      GO TO 90                                                          JAC00075
   80 CA=AKK/DEN                                                        JAC00076
      CG=-AJJ/DEN                                                       JAC00077
C                                                                       JAC00078
C     PERFORM THE GENERALIZED ROTATION TO ZERO THE PRESENT OFF-DIAGONAL JAC00079
C     ELEMENT                                                           JAC00080
C                                                                       JAC00081
   90 IF (N-2) 100,190,100                                              JAC00082
  100 IF (JM1-1) 130,110,110                                            JAC00083
  110 DO 120 I=1,JM1                                                    JAC00084
      IM1=I - 1                                                         JAC00085
      IJ=IM1*N - IM1*I/2 + J                                            JAC00086
      IK=IM1*N - IM1*I/2 + K                                            JAC00087
      AJ=A(IJ)                                                          JAC00088
      BJ=B(IJ)                                                          JAC00089
      AK=A(IK)                                                          JAC00090
      BK=B(IK)                                                          JAC00091
      A(IJ)=AJ + CG*AK                                                  JAC00092
      B(IJ)=BJ + CG*BK                                                  JAC00093
      A(IK)=AK + CA*AJ                                                  JAC00094
  120 B(IK)=BK + CA*BJ                                                  JAC00095
  130 IF (KP1-N) 140,140,160                                            JAC00096
  140 LJI=JM1*N - JM1*J/2                                               JAC00097
      LKI=KM1*N - KM1*K/2                                               JAC00098
      DO 150 I=KP1,N                                                    JAC00099
      JI=LJI + I                                                        JAC00100
      KI=LKI + I                                                        JAC00101
      AJ=A(JI)                                                          JAC00102
      BJ=B(JI)                                                          JAC00103
      AK=A(KI)                                                          JAC00104
      BK=B(KI)                                                          JAC00105
      A(JI)=AJ + CG*AK                                                  JAC00106
      B(JI)=BJ + CG*BK                                                  JAC00107
      A(KI)=AK + CA*AJ                                                  JAC00108
  150 B(KI)=BK + CA*BJ                                                  JAC00109
  160 IF (JP1-KM1) 170,170,190                                          JAC00110
  170 LJI=JM1*N - JM1*J/2                                               JAC00111
      DO 180 I=JP1,KM1                                                  JAC00112
      JI=LJI + I                                                        JAC00113
      IM1=I - 1                                                         JAC00114
      IK=IM1*N - IM1*I/2 + K                                            JAC00115
      AJ=A(JI)                                                          JAC00116
      BJ=B(JI)                                                          JAC00117
      AK=A(IK)                                                          JAC00118
      BK=B(IK)                                                          JAC00119
      A(JI)=AJ + CG*AK                                                  JAC00120
      B(JI)=BJ + CG*BK                                                  JAC00121
      A(IK)=AK + CA*AJ                                                  JAC00122
  180 B(IK)=BK + CA*BJ                                                  JAC00123
  190 AK=A(KK)                                                          JAC00124
      BK=B(KK)                                                          JAC00125
      A(KK)=AK + 2.*CA*A(JK) + CA*CA*A(JJ)                              JAC00126
      B(KK)=BK + 2.*CA*B(JK) + CA*CA*B(JJ)                              JAC00127
      A(JJ)=A(JJ) + 2.*CG*A(JK) + CG*CG*AK                              JAC00128
      B(JJ)=B(JJ) + 2.*CG*B(JK) + CG*CG*BK                              JAC00129
      A(JK)=0.                                                          JAC00130
      B(JK)=0.                                                          JAC00131
C                                                                       JAC00132
C     UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION                 JAC00133
C                                                                       JAC00134
      DO 200 I=1,N                                                      JAC00135
      XJ=X(I,J)                                                         JAC00136
      XK=X(I,K)                                                         JAC00137
      X(I,J)=XJ + CG*XK                                                 JAC00138
  200 X(I,K)=XK + CA*XJ                                                 JAC00139
  210 CONTINUE                                                          JAC00140
C                                                                       JAC00141
C     UPDATE THE EIGENVALUES AFTER EACH SWEEP                           JAC00142
C                                                                       JAC00143
      II=1                                                              JAC00144
      DO 220 I=1,N                                                      JAC00145
      IF (A(II).GT.0. .AND. B(II).GT.0.) GO TO 215                      JAC00146
      WRITE (IOUT,2020) II,A(II),B(II)                                  JAC00147
      GO TO 800                                                         JAC00148
  215 EIGV(I)=A(II)/B(II)                                               JAC00149
  220 II=II + N1 - I                                                    JAC00150
      IF (IFPR.EQ.0) GO TO 230                                          JAC00151
      WRITE (IOUT,2030)                                                 JAC00152
      WRITE (IOUT,2010) (EIGV(I),I=1,N)                                 JAC00153
C                                                                       JAC00154
C     CHECK FOR CONVERGENCE                                             JAC00155
C                                                                       JAC00156
  230 DO 240 I=1,N                                                      JAC00157
      TOL=RTOL*D(I)                                                     JAC00158
      DIF=ABS(EIGV(I)-D(I))                                             JAC00159
      IF (DIF.GT.TOL) GO TO 280                                         JAC00160
  240 CONTINUE                                                          JAC00161
C                                                                       JAC00162
C     CHECK ALL OFF-DIAGONAL ELEMENTS TO SEE IF ANOTHER SWEEP IS        JAC00163
C     REQUIRED                                                          JAC00164
C                                                                       JAC00165
      EPS=RTOL**2                                                       JAC00166
      DO 250 J=1,NR                                                     JAC00167
      JM1=J - 1                                                         JAC00168
      JP1=J + 1                                                         JAC00169
      LJK=JM1*N - JM1*J/2                                               JAC00170
      JJ=LJK + J                                                        JAC00171
      DO 250 K=JP1,N                                                    JAC00172
      KM1=K - 1                                                         JAC00173
      JK=LJK + K                                                        JAC00174
      KK=KM1*N - KM1*K/2 + K                                            JAC00175
      EPSA=(A(JK)/A(JJ))*(A(JK)/A(KK))                                  JAC00176
      EPSB=(B(JK)/B(JJ))*(B(JK)/B(KK))                                  JAC00177
      IF (EPSA.LT.EPS .AND. EPSB.LT.EPS) GO TO 250                      JAC00178
      GO TO 280                                                         JAC00179
  250 CONTINUE                                                          JAC00180
C                                                                       JAC00181
C     SCALE EIGENVECTORS                                                JAC00182
C                                                                       JAC00183
  255 II=1                                                              JAC00184
      DO 275 I=1,N                                                      JAC00185
      BB=SQRT(B(II))                                                    JAC00186
      DO 270 K=1,N                                                      JAC00187
  270 X(K,I)=X(K,I)/BB                                                  JAC00188
  275 II=II + N1 - I                                                    JAC00189
      GO TO 900                                                         JAC00190
C                                                                       JAC00191
C     UPDATE  D  MATRIX AND START NEW SWEEP, IF ALLOWED                 JAC00192
C                                                                       JAC00193
  280 DO 290 I=1,N                                                      JAC00194
  290 D(I)=EIGV(I)                                                      JAC00195
      IF (NSWEEP.LT.NSMAX) GO TO 40                                     JAC00196
      GO TO 255                                                         JAC00197
C                                                                       JAC00198
  800 STOP                                                              JAC00199
  900 RETURN                                                            JAC00200
C                                                                       JAC00201
 2000 FORMAT (//,' SWEEP NUMBER IN *JACOBI* = ',I8)                     JAC00202
 2010 FORMAT (' ',6E20.12)                                              JAC00203
 2020 FORMAT (' *** ERROR *** SOLUTION STOP',/,                         JAC00204
     1        ' MATRICES NOT POSITIVE DEFINITE',/,                      JAC00205
     2        ' II = ',I8,' A(II) = ',E20.12,' B(II) = ',E20.12)        JAC00206
 2030 FORMAT (/,' CURRENT EIGENVALUES IN *JACOBI* ARE',/)               JAC00207
      END                                                               JAC00208
