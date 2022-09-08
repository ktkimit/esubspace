      SUBROUTINE STARTV (A,B,MAXA,R,NN,NWM,NC,IOUT)                     STA00001
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . STA00002
C                                                                       STA00003
C     ESTABLISH STARTING ITERATION VECTORS                              STA00004
C                                                                       STA00005
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . STA00006
C                                                                       STA00007
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               STA00008
      DIMENSION A(1),B(1),MAXA(1),R(NN,NC),TT(NN),W(NN)                 STA00009
C                                                                       STA00010
      ND=NN/NC                                                          STA00011
      IF (NWM.GT.NN) GO TO 4                                            STA00012
      J=0                                                               STA00013
      DO 6 I=1,NN                                                       STA00014
      II=MAXA(I)                                                        STA00015
      R(I,1)=B(I)                                                       STA00016
      IF (B(I).GT.0) J=J + 1                                            STA00017
    6 W(I)=B(I)/A(II)                                                   STA00018
      IF (NC.LE.J) GO TO 16                                             STA00019
      WRITE (IOUT,1007)                                                 STA00020
      GO TO 800                                                         STA00021
    4 DO 10 I=1,NN                                                      STA00022
      II=MAXA(I)                                                        STA00023
      R(I,1)=B(II)                                                      STA00024
   10 W(I)=B(II)/A(II)                                                  STA00025
   16 DO 20 J=2,NC                                                      STA00026
      DO 20 I=1,NN                                                      STA00027
   20 R(I,J)=0.                                                         STA00028
C                                                                       STA00029
      L=NN - ND                                                         STA00030
      DO 30 J=2,NC                                                      STA00031
      RT=0.                                                             STA00032
      DO 40 I=1,L                                                       STA00033
      IF (W(I).LT.RT) GO TO 40                                          STA00034
      RT=W(I)                                                           STA00035
      IJ=I                                                              STA00036
   40 CONTINUE                                                          STA00037
      DO 50 I=L,NN                                                      STA00038
      IF (W(I).LE.RT) GO TO 50                                          STA00039
      RT=W(I)                                                           STA00040
      IJ=I                                                              STA00041
   50 CONTINUE                                                          STA00042
      TT(J)=FLOAT(IJ)                                                   STA00043
      W(IJ)=0.                                                          STA00044
      L=L - ND                                                          STA00045
   30 R(IJ,J)=1.                                                        STA00046
C                                                                       STA00047
      WRITE (IOUT,1008)                                                 STA00048
      WRITE (IOUT,1002) (TT(J),J=2,NC)                                  STA00049
C                                                                       STA00050
C     A RANDOM VECTOR IS ADDED TO THE LAST VECTOR                       STA00051
C                                                                       STA00052
      PI=3.141592654D0                                                  STA00053
      XX=0.5D0                                                          STA00054
      DO 60 K=1,NN                                                      STA00055
      XX=(PI + XX)**5                                                   STA00056
      IX=INT(XX)                                                        STA00057
      XX=XX - FLOAT(IX)                                                 STA00058
   60 R(K,NC)=R(K,NC) + XX                                              STA00059
      GO TO 900                                                         STA00060
C                                                                       STA00061
  800 STOP                                                              STA00062
  900 RETURN                                                            STA00063
C                                                                       STA00064
 1002 FORMAT (' ',10F10.0)                                              STA00065
 1007 FORMAT (///,' STOP, NC IS LARGER THAN THE NUMBER OF MASS ',       STA00066
     1        'DEGREES OF FREEDOM')                                     STA00067
 1008 FORMAT (///,' DEGREES OF FREEDOM EXCITED BY UNIT STARTING ',      STA00068
     1        'ITERATION VECTORS')                                      STA00069
      END                                                               STA00070
