      SUBROUTINE SCHECK (EIGV,RTOLV,BUP,BLO,BUPC,NEIV,NC,NEI,RTOL,      SCH00001
     1                   SHIFT,IOUT)                                    SCH00002
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SCH00003
C .                                                                   . SCH00004
C .   P R O G R A M                                                   . SCH00005
C .        TO EVALUATE SHIFT FOR STURM SEQUENCE CHECK                 . SCH00006
C .                                                                   . SCH00007
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . SCH00008
C                                                                       SCH00009
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               SCH00010
      DIMENSION EIGV(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC),NEIV(NC)    SCH00011
C                                                                       SCH00012
      FTOL=0.01                                                         SCH00013
C                                                                       SCH00014
      DO 100 I=1,NC                                                     SCH00015
      BUP(I)=EIGV(I)*(1.+FTOL)                                          SCH00016
  100 BLO(I)=EIGV(I)*(1.-FTOL)                                          SCH00017
      NROOT=0                                                           SCH00018
      DO 120 I=1,NC                                                     SCH00019
  120 IF (RTOLV(I).LT.RTOL) NROOT=NROOT + 1                             SCH00020
      IF (NROOT.GE.1) GO TO 200                                         SCH00021
      WRITE (IOUT,1010)                                                 SCH00022
      GO TO 800                                                         SCH00023
C                                                                       SCH00024
C      FIND UPPER BOUNDS ON EIGENVALUE CLUSTERS                         SCH00025
C                                                                       SCH00026
  200 DO 240 I=1,NROOT                                                  SCH00027
  240 NEIV(I)=1                                                         SCH00028
      IF (NROOT.NE.1) GO TO 260                                         SCH00029
      BUPC(1)=BUP(1)                                                    SCH00030
      LM=1                                                              SCH00031
      L=1                                                               SCH00032
      I=2                                                               SCH00033
      GO TO 295                                                         SCH00034
  260 L=1                                                               SCH00035
      I=2                                                               SCH00036
  270 IF (BUP(I-1).LE.BLO(I)) GO TO 280                                 SCH00037
      NEIV(L)=NEIV(L) + 1                                               SCH00038
      I=I + 1                                                           SCH00039
      IF (I.LE.NROOT) GO TO 270                                         SCH00040
  280 BUPC(L)=BUP(I-1)                                                  SCH00041
      IF (I.GT.NROOT) GO TO 290                                         SCH00042
      L=L + 1                                                           SCH00043
      I=I + 1                                                           SCH00044
      IF (I.LE.NROOT) GO TO 270                                         SCH00045
      BUPC(L)=BUP(I-1)                                                  SCH00046
  290 LM=L                                                              SCH00047
      IF (NROOT.EQ.NC) GO TO 300                                        SCH00048
  295 IF (BUP(I-1).LE.BLO(I)) GO TO 300                                 SCH00049
      IF (RTOLV(I).GT.RTOL) GO TO 300                                   SCH00050
      BUPC(L)=BUP(I)                                                    SCH00051
      NEIV(L)=NEIV(L) + 1                                               SCH00052
      NROOT=NROOT + 1                                                   SCH00053
      IF (NROOT.EQ.NC) GO TO 300                                        SCH00054
      I=I + 1                                                           SCH00055
      GO TO 295                                                         SCH00056
C                                                                       SCH00057
C      FIND SHIFT                                                       SCH00058
C                                                                       SCH00059
  300 WRITE (IOUT,1020)                                                 SCH00060
      WRITE (IOUT,1005) (BUPC(I),I=1,LM)                                SCH00061
      WRITE (IOUT,1030)                                                 SCH00062
      WRITE (IOUT,1006) (NEIV(I),I=1,LM)                                SCH00063
      LL=LM - 1                                                         SCH00064
      IF (LM.EQ.1) GO TO 310                                            SCH00065
  330 DO 320 I=1,LL                                                     SCH00066
  320 NEIV(L)=NEIV(L) + NEIV(I)                                         SCH00067
      L=L - 1                                                           SCH00068
      LL=LL - 1                                                         SCH00069
      IF (L.NE.1) GO TO 330                                             SCH00070
  310 WRITE (IOUT,1040)                                                 SCH00071
      WRITE (IOUT,1006) (NEIV(I),I=1,LM)                                SCH00072
      L=0                                                               SCH00073
      DO 340 I=1,LM                                                     SCH00074
      L=L + 1                                                           SCH00075
      IF (NEIV(I).GE.NROOT) GO TO 350                                   SCH00076
  340 CONTINUE                                                          SCH00077
  350 SHIFT=BUPC(L)                                                     SCH00078
      NEI=NEIV(L)                                                       SCH00079
      GO TO 900                                                         SCH00080
C                                                                       SCH00081
  800 STOP                                                              SCH00082
  900 RETURN                                                            SCH00083
C                                                                       SCH00084
 1005 FORMAT (' ',6E22.14)                                              SCH00085
 1006 FORMAT (' ',6I22)                                                 SCH00086
 1010 FORMAT (' *** ERROR ***  SOLUTION STOP IN *SCHECK*',/,            SCH00087
     1        ' NO EIGENVALUES FOUND',/)                                SCH00088
 1020 FORMAT (///,' UPPER BOUNDS ON EIGENVALUE CLUSTERS')               SCH00089
 1030 FORMAT (//,' NO. OF EIGENVALUES IN EACH CLUSTER')                 SCH00090
 1040 FORMAT (' NO. OF EIGENVALUES LESS THAN UPPER BOUNDS')             SCH00091
      END                                                               SCH00092
