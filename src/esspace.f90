!===============================================================================
! M O D U L E
!   Enriched subspace iteration method to solve for the smallest eigenvalues 
!   and corresponding eigenvectors of the generalized symmetric positive 
!   definite eigenproblem
!   
!   References:
!   [1] K.T. Kim and K.J. Bathe, The Bathe subspace iteration method enriched 
!       by turning vectors, Computers & Structures, 2017
!   [2] K.J. Bathe, Finite element procedures, 2nd ed., 2014
!
! NEED:
!   STARTV
!   DECOMP
!   REDBAK
!   MULT
!   JACOBI
!   SCHECK
!
! AUTHOR     : Ki-Tae Kim
!              Klaus-Jürgen Bathe
! EMAIL      : ktkim@mit.edu or qlsn55@gmail.com
! AFFILIATION: Massachusetts Institute of Technology
! LANGUAGE   : FORTRAN 90
! CREATED    : Thu 20 Apr 2017 02:06:10 PM EDT
!===============================================================================
MODULE ESSPACE
    IMPLICIT NONE

    PUBLIC ESUBSPACE_ITER
    PRIVATE ESUBSPACE_PROJ, CHECKNULL, DETERMINE_CONV, BSUBSPACE_ITER1,&
        INSERTION_SORT, PRINT_HEADING, PRINT_SETTING

CONTAINS
    !===========================================================================
    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Driving subroutine of the enriched subspace iteration method
    !
    ! NEED:
    !   PRINT_HEADING
    !   PRINT_SETTING
    !   STARTV
    !   DECOMP
    !   BSUBSPACE_ITER1
    !   ESUBSPACE_PROJ
    !   JACOBI
    !   DETERMINE_CONV
    !   INSERTION_SORT
    !   SCHECK
    !
    ! INPUT:
    !   KK(nwk)  = factorized stiffness matrix (L*D*L^T by COLSOL) in skyline
    !              format
    !   MM(nwm)  = mass matrix in skyline format
    !   DK(nn+1) = address of diagonal elements of factorized stiffness matrix
    !   DM(nn+1) = address of diagonal elements of mass matrix
    !   nn       = order of system
    !   nwk      = number of elements below skyline of stiffness matrix
    !   nwm      = number of elements below skyline of mass matrix;
    !              EQ. nwk for consistent mass matrix
    !              EQ. nn for lumped mass matrix
    !   nroot    = number of required eigenpairs
    !   nitm     = maximum number of iteration permitted
    !   nq       = number of iteration vectors
    !   nnq      = nq*(nq + 1)/2
    !   tolc     = tolerance for convergence;
    !              usually set to 1.0E-6 or smaller
    !   tolt     = tolerance for selecting the turning vectors
    !              usually set to 1.0E-8
    !   ifss     = flag for Sturm sequence check;
    !              EQ. 0 no check
    !              EQ. 1 check
    !   ilog     = unit used for formatted file that prints results during
    !              iteration;
    !              EQ. 0 no printing
    !   nstiff   = unit used for scratch file of stiffness matrix
    !
    ! OUTPUT:
    !   EIGV(nq)        = calculated eigenvalues;
    !                     first nroot values are the required eigenvalues
    !                     in ascending order
    !   RV(nn,nq+nroot) = working vectors;
    !                     first nroot columns are the corresponding required 
    !                     eigenvectors
    !---------------------------------------------------------------------------
    SUBROUTINE ESUBSPACE_ITER(KK,MM,DK,DM,nn,nwk,nwm,nroot,nitm,nq,nnq,tolc,&
            tolt,ifss,ilog,nstiff,EIGV,RV)
        integer, intent(in) :: nn, nwk, nwm, nroot, nitm, nq, nnq, ifss,&
            ilog, nstiff
        integer, intent(in) :: DK(nn+1), DM(nn+1)
        double precision, intent(in)    :: tolc, tolt
        double precision, intent(inout) :: KK(nwk), MM(nwm)
        double precision, intent(out)   :: EIGV(nq), RV(nn,nq+nroot)

        integer :: i, ii, j, k, nn1, nq2, np, np0, np1, nx, nr1, nnr1, nr2,&
            nnx, nt, tp, ifpr, nsmax, iconv, nei, nsch, nmis
        integer :: NEIV(nq), TID(nq)
        double precision :: tolj, tol24, val, vdot, shift
        double precision :: RTOLV(nq), KR(nnq), MR(nnq), WN(nn), WV(nn),&
            QQ(nq,nq), WQ(nq), A(nq,nq), B(nq,nq), C(nq,nq), D(nq,nq),&
            ALPHA(nnq), BUP(nq), BLO(nq), BUPC(nq), TVEC(nq), TEIGV(nq),&
            TRTOLV(nq)

        !---------------!
        ! Print heading !
        !---------------!
        call PRINT_HEADING(ilog)

        !--------------------!
        ! Print input values !
        !--------------------!
        call PRINT_SETTING(KK,MM,DK,DM,RV,nn,nroot,nq,tolc,tolt,ilog)

        !----------------!
        ! Initialization !
        !----------------!
        nn1  = nn + 1
        nq2  = nq*2
        np   = 0
        np1  = 1
        nx   = nq
        nr1  = nx/2
        nnr1 = nr1*(nr1 + 1)/2
        nr2  = nx - nr1
        nnx  = nx*(nx + 1)/2

        iconv = 0
        nsmax = 12
        tolj  = 1.0D-12
        tol24 = 1.0D-24

        if (ilog == 0) then
            ifpr = 0
        else
            ifpr = 1
        end if

        !----------------------------------------------------!
        ! Construct M-orthonormal starting iteration vectors !
        !----------------------------------------------------!
        write(ilog,100)

        ! Generate linearly independent starting vectors
        call STARTV(KK,MM,DK,RV(:,1:nq),nn,nwm,nq,ilog)

        ! Save KK
        write(nstiff,*) KK

        ! LDL^T decomposition of KK
        call DECOMP(KK,DK,nn,0,ilog)

        ! Run a single basic subspace iteration
        call BSUBSPACE_ITER1(KK,MM,DK,DM,nn,nn1,nwk,nwm,nq,nnq,RV(:,1:nq),&
            EIGV,RTOLV,0,KR,MR,WN,QQ,WQ)

        !----------------------!
        ! Start iteration loop !
        !----------------------!
        write(ilog,101)

        do i = 1, nitm
            write(ilog,102) i

            !----------------------------------------------------------------!
            ! Enrich the subspace and project the matrices onto the subspace !
            !----------------------------------------------------------------!
            call ESUBSPACE_PROJ(KK,MM,DK,DM,nn,nn1,nwk,nwm,EIGV(1:np),&
                RV(:,1:np),RV(:,np1:nq),np,nx,nnx,nr1,nnr1,nr2,nq,nnq,&
                nt,tolt,RV(:,nq+1:nq+nr1),KR,MR,TID(1:nr1),A(1:nr1,1:nr1),&
                B(1:nr2,1:nr1),C(1:nr1,1:nr1),D(1:np,1:nr1),WN,WV,WQ,&
                ALPHA(1:nnr1))

            !--------------------------!
            ! Print projected matrices !
            !--------------------------!
            if (ilog .ne. 0) then
                write(ilog,103)
                do j = 1, nq
                    tp = (nq2 - j)*(j-1)/2
                    do k = j, nq
                        write(ilog,200,advance='no') KR(k+tp)
                    end do
                    write(ilog,300)
                end do

                write(ilog,104)
                do j = 1, nq
                    tp = (nq2 - j)*(j-1)/2
                    do k = j, nq
                        write(ilog,200,advance='no') MR(k+tp)
                    end do
                    write(ilog,300)
                end do
            end if

            !-----------------------------------------------------!
            ! Solve for the eigensystem of the projected matrices !
            !-----------------------------------------------------!
            call JACOBI(KR,MR,QQ,EIGV,WQ,nq,nnq,tolj,nsmax,ifpr,ilog)

            !--------------------------!
            ! Print projected matrices !
            !--------------------------!
            if (ilog .ne. 0) then
                write(ilog,105)
                do j = 1, nq
                    tp = (nq2 - j)*(j-1)/2
                    do k = j, nq
                        write(ilog,200,advance='no') KR(k+tp)
                    end do
                    write(ilog,300)
                end do

                write(ilog,106)
                do j = 1, nq
                    tp = (nq2 - j)*(j-1)/2
                    do k = j, nq
                        write(ilog,200,advance='no') MR(k+tp)
                    end do
                    write(ilog,300)
                end do
            end if

            !------------------------!
            ! Calculate error bounds !
            !------------------------!
            do j = np1, nq
                vdot = 0.0D0
                do k = np1, nq
                    vdot = vdot + QQ(k,j)*QQ(k,j)
                end do

                if (vdot == 0.0D0) then
                    stop 'VDOT IN ESUBSPACE_ITER SHOULD NOT BE .EQ. ZERO'
                end if

                val = vdot - EIGV(j)*EIGV(j)
                val = val / vdot
                val = MAX(val,tol24)
                RTOLV(j) = SQRT(val)
            end do

            !-----------------------------------------------!
            ! Print calculated eigenvalues and error bounds !
            !-----------------------------------------------!
            if (ilog .ne. 0) then
                write(ilog,107)
                write(ilog,203) EIGV
            end if

            if (ilog .ne. 0) then
                write(ilog,108)
                write(ilog,202) RTOLV
            end if

            !-------------------------------------------------------------!
            ! Since the projected matrices have been calculated with the  !
            ! factor of 1/EIGV(1:np) for the converged iteration vectors, !
            ! scaling is now corrected for                                !
            !-------------------------------------------------------------!
            do j = 1, np
                if (EIGV(j) == 0.0D0) then
                    stop 'EIGV(:) IN ESUBSPACE_ITER SHOULD NOT BE .EQ. ZERO'
                end if

                QQ(j,:) = QQ(j,:) / EIGV(j)
            end do

            !---------------------------------------------------!
            ! Determine converged eigenpairs and rearrange them !
            !---------------------------------------------------!
            np0 = np
            call DETERMINE_CONV(nq,nroot,EIGV,QQ,RTOLV,np,tolc,iconv,TEIGV,&
                TRTOLV,TVEC)

            !-----------------!
            ! Print np and nt !
            !-----------------!
            write(ilog,109) np0
            write(ilog,110) np - np0
            write(ilog,111) nt

            if (iconv == 0) then
                !----------------!
                ! Update indices !
                !----------------!
                np1  = np + 1
                nx   = nq - np
                nnx  = nx*(nx +  1)/2
                nr1  = nx/2
                nnr1 = nr1*(nr1 + 1)/2
                nr2  = nx - nr1

                !--------------------------------------------------!
                ! Arrange calculated eigenpairs in ascending order !
                !--------------------------------------------------!
                call INSERTION_SORT(EIGV(1:np),QQ(:,1:np),np,nq,TVEC)
                call INSERTION_SORT(EIGV(np1:nq),QQ(:,np1:nq),nx,nq,TVEC)
            end if

            !-----------------------------------------------!
            ! Improve the approximation to the eigenvectors !
            !-----------------------------------------------!
            do ii = 1, nn
                WQ = RV(ii,1:nq)
                do j = 1, nq
                    RV(ii,j) = 0.0D0
                    do k = 1, nq
                        RV(ii,j) = RV(ii,j) + WQ(k)*QQ(k,j)
                    end do
                end do
            end do

            !----------------------------------------------------------------!
            ! If the calculated smallest eigenpairs converged, exit the loop !
            !----------------------------------------------------------------!
            if (iconv == 1) then
                write(ilog,112) tolc
                exit
            end if
        end do

        !------------------------------------------------------!
        ! Print comments when convergence has not been reached !
        !------------------------------------------------------!
        if (iconv == 0) then
            write(ilog,113)
        end if

        !--------------------------------------------------------------------!
        ! Perform a single basic subspace iteration to obtain approximate    !
        ! eigenvectors from iteration vectors (iteration vectors in the loop !
        ! are of the form M*X)                                               !
        !--------------------------------------------------------------------!
        write(ilog,114)
        call BSUBSPACE_ITER1(KK,MM,DK,DM,nn,nn1,nwk,nwm,nq,nnq,RV(:,1:nq),&
            EIGV,RTOLV,1,KR,MR,WN,QQ,WQ)

        !---------------!
        ! Print results !
        !---------------!
        write(ilog,115)
        write(ilog,203) EIGV(1:nroot)
        write(ilog,108)
        write(ilog,202) RTOLV(1:nroot)

        !----------------------------!
        ! Apply Sturm sequence check !
        !----------------------------!
        write(ilog,116)
        if (ifss == 1 .and. iconv == 1) then
            call SCHECK(EIGV,RTOLV,BUP,BLO,BUPC,NEIV,nq,nei,tolc,shift,ilog)
            write(ilog,117) shift

            ! Read stiffness matrix
            rewind(unit=nstiff)
            read(nstiff,*) KK

            ! Shift stiffness matrix
            if (nwm > nn) then
                do i = 1, nwk
                    KK(i) = KK(i) - MM(i)*shift
                end do
            else
                do i = 1, nn
                    ii = DK(i)
                    KK(ii) = KK(ii) - MM(i)*shift
                end do
            end if

            ! Factorize shifted stiffness matrix
            call DECOMP(KK,DK,nn,1,ilog)

            ! Count number of negative diagonal elements
            nsch = 0
            do i = 1, nn
                ii = DK(i)
                if (KK(ii) < 0.0D0) then
                    nsch = nsch + 1
                end if
            end do

            if (nsch == nei) then
                write(ilog,119) nsch
            else
                nmis = nsch - nei
                write(ilog,118) nmis
            end if

        end if

        !------------------!
        ! Formats to write !
        !------------------!
        100 format(//, X, "----- ESTABLISHMENT OF M-ORTHONORMAL ",& 
            "STARTING ITERATION VECTORS -----")
        101 format(//, X, "----- START OF ITERATION ------")
        102 format(//, X, "-- I T E R A T I O N  N U M B E R ",&
            I2)
        103 format(/, X, "PROJECTION OF K (MATRIX KR)")
        104 format(/, X, "PROJECTION OF M (MATRIX MR)")
        105 format(/, X, "KR AFTER JACOBI DIAGONALIZATION")
        106 format(/, X, "MR AFTER JACOBI DIAGONALIZATION")
        107 format(/, X, "EIGENVALUES OF PROJECTED MATRICES")
        108 format(/, X, "ERROR BOUNDS REACHED ON EIGENVALUES")
        109 format(//, X, "NUMBER OF PREVIOUSLY CONVERGED VECTORS (NP)", 14X,&
            "=", 2X, I8)
        110 format(X, "NUMBER OF ADDITIONALLY CONVERGED VECTORS AFTER ",& 
            "ITERATION =", 2X, I8)
        111 format(X, "NUMBER OF TURNING VECTORS USED (NT)", 22X, "=", 2X, I8)
        112 format(//, X, "*** CONVERGENCE REACHED FOR TOL. = ", 2X, ES11.4)
        113 format(//, X, "*** NO CONVERGENCE IN MAXIMUM NUMBER OF ITERATIONS",&
            " PERMITTED", /, 5X, "WE ACCEPT CURRENT ITERATION VALUES", /,&
            5X, "THE STURM SEQUENCE CHECK IS NOT PERFORMED")
        114 format(//, X, "----- SINGLE BASIC SUBSPACE ITERATION ------")
        115 format(//, X, "THE CALCULATED EIGENVALUES ARE")
        116 format(//, X, "----- STURM SEQUENCE CHECK -----")
        117 format(//, X, "CHECK APPLIED AT SHIFT = ", ES22.15)
        118 format(//, X, "THERE ARE ", I8, " EIGENVALUES MISSING")
        119 format(//, X, "WE FOUND THE LOWEST ", I8, " EIGENVALUES")
        200 format(X, ES11.4)
        201 format(X, ES22.15)
        202 format(5(X,ES11.4))
        203 format(5(X,ES22.15))
        300 format()
    END SUBROUTINE ESUBSPACE_ITER

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Construct the enriched subspace and project stiffness and mass matrices 
    !   onto this subspace
    !
    ! NEED:
    !   REDBAK
    !   MULT
    !   CHECKNULL
    !
    ! INPUT:
    !   KK(nwk)    = factorized stiffness matrix (L*D*L^T by COLSOL) in skyline
    !                format
    !   MM(nwm)    = mass matrix in skyline format
    !   DK(nn1)    = address of diagonal elements of factorized stiffness matrix
    !   DM(nn1)    = address of diagonal elements of mass matrix
    !   nn         = order of system
    !   nn1        = nn + 1
    !   nwk        = number of elements below skyline of stiffness matrix
    !   nwm        = number of elements below skyline of mass matrix;
    !                EQ. nwk for consistent mass matrix
    !                EQ. nn for lumped mass matrix
    !   LAMBDA(np) = approximate eigenvalues converged to the required
    !                tolerance
    !   PHI(nn,np) = iteration vectors converged to the required tolerance;
    !                passed as MM*PHI
    !   X(nn,nx)   = iteration vectors to be updated;
    !                passed as MM*X
    !   np         = number of converged iteration vectors
    !   nx         = number of iteration vectors to be updated;
    !                EQ. nr1 + nr2
    !   nnx        = nx*(nx + 1)/2
    !   nr1        = number of iteration vectors used to compute the turning 
    !                vectors
    !   nnr1       = nr1*(nr1 + 1)/2
    !   nr2        = number of the rest of iteration vectors;
    !                GE. nr1
    !   nq         = np + nx
    !   nnq        = nq*(nq + 1)/2
    !   tolt       = tolerance for selecting the turning vectors
    !
    ! OUTPUT:
    !   nt          = number of turning vectors actually used
    !   X(nn,nx)    = updated iteration vectors;
    !                 passed as MM*X
    !   MX(nn,nr1)  = working vectors
    !   KR(nnq)     = projection of stiffness matrix
    !   MR(nnq)     = projection of mass matrix
    !                 Both KR and MR are stored as packed format, e.g.,
    !                 KR(i+(2*nq-j)*(j-1)/2) for j <= i is the (i,j) element of
    !                 the corresponding matrix
    !   TID(nr1)    = working array
    !   A(nr1,nr1)  = working array
    !   B(nr2,nr1)  = working array
    !   C(nr1,nr1)  = working array
    !   D(np,nr1)   = working array
    !   ALPHA(nnr1) = working array
    !   WN(nn)      = working array
    !   WV(nn)      = working array
    !   WQ(nq)      = working array
    !---------------------------------------------------------------------------
    SUBROUTINE ESUBSPACE_PROJ(KK,MM,DK,DM,nn,nn1,nwk,nwm,LAMBDA,PHI,X,np,nx,&
            nnx,nr1,nnr1,nr2,nq,nnq,nt,tolt,MX,KR,MR,TID,A,B,C,D,WN,WV,WQ,ALPHA)
        integer, intent(in)  :: nn, nn1, nwk, nwm, np, nx, nnx, nr1, nnr1,&
            nr2, nq, nnq
        integer, intent(in)  :: DK(nn1), DM(nn1)
        integer, intent(out) :: nt, TID(nr1)
        double precision, intent(in)    :: tolt
        double precision, intent(in)    :: KK(nwk), MM(nwm), LAMBDA(np),&
            PHI(nn,np)
        double precision, intent(inout) :: X(nn,nx)
        double precision, intent(out)   :: MX(nn,nr1), KR(nnq), MR(nnq),&
            A(nr1,nr1), B(nr2,nr1), C(nr1,nr1), D(np,nr1), WN(nn), WV(nn),&
            WQ(nq), ALPHA(nnr1)

        integer :: i, j, k, ad, nxp1, iftv, nxmnt, nx2, ntp1, cc, ad2,&
            ct, nq2, np1, ctt, sct
        double precision :: normt, val

        !----------------!
        ! Initialization !
        !----------------!
        nt    = 0
        TID   = 0
        ALPHA = 0.0D0
        A     = 0.0D0
        B     = 0.0D0
        C     = 0.0D0
        D     = 0.0D0

        !----------------------------------------------------------------------!
        ! Perform the inverse iteration on X(:,1:nr1) and calculate A, B and D !
        !                                                                      !
        !   A = X(:,1:nr1)^T * BX,                                             !
        !   B = X(:,nr1+1:nx)^T * BX,                                          !
        !   D = PHI(:,1:np)^T * BX,                                            !
        ! where                                                                !
        !   BX = KK^(-1) * X(:,1:nr1)                                          !
        !----------------------------------------------------------------------!
        do j = 1, nr1
            !----------------!
            ! WN = KK^(-1)*X !
            !----------------!
            WN = X(:,j)
            call REDBAK(KK,WN,DK,nn)

            !------------!
            ! Compute Aj !
            !------------!
            do i = j, nr1
                do k = 1, nn
                    A(i,j) = A(i,j) + X(k,i)*WN(k)
                end do

                if (i > j) then
                    A(j,i) = A(i,j)
                end if
            end do

            !------------!
            ! Compute Dj !
            !------------!
            do i = 1, np
                do k = 1, nn
                    D(i,j) = D(i,j) + PHI(k,i)*WN(k)
                end do
            end do

            !------------!
            ! Compute Bj !
            !------------!
            do i = 1, nr2
                ad = nr1 + i
                do k = 1, nn
                    B(i,j) = B(i,j) + X(k,ad)*WN(k)
                end do
            end do

            !--------------------------!
            ! Update iteration vectors !
            !--------------------------!
            MX(:,j) = X(:,j)
            X(:,j)  = WN
        end do

        !-----------------------------------------------------------!
        ! Calculate MM*X and C, and check the amount of turning and !
        ! choose numerically stable turning vectors                 !
        !                                                           !
        !   C = X^T * MM * X                                        !
        !-----------------------------------------------------------!
        nxp1 = nx + 1
        do j = nr1, 1, -1
            !-----------!
            ! WN = MM*X !
            !-----------!
            call MULT(WN,MM,X(:,j),DM,nn,nwm)

            !------------!
            ! Compute Cj !
            !------------!
            do i = 1, j
                do k = 1, nn
                    C(i,j) = C(i,j) + X(k,i)*WN(k)
                end do

                if (i < j) then
                    C(j,i) = C(i,j)
                end if
            end do

            !---------------------------!
            ! Step (a) in Reference [1] !
            !---------------------------!
            call CHECKNULL(A,B,C,D,j,nt,TID,ALPHA,tolt,iftv,WQ)

            if (iftv == 1) then
                X(:,nxp1-nt) = X(:,j)
            end if

            !--------------------------!
            ! Update iteration vectors !
            !--------------------------!
            X(:,j) = WN
        end do

        !------------------------------------!
        ! Calculate selected turning vectors !
        !------------------------------------!
        nxmnt = nx - nt
        nx2   = nr2 - nt
        ntp1  = nt + 1
        do j = nt, 1, -1
            cc = nxmnt + j
            ad = TID(ntp1-j)
            WV = X(:,ad)

            !----------------------!
            ! Orthogonalize to PHI !
            !----------------------!
            do i = 1, np
                WV = WV - PHI(:,i)*D(i,ad)
            end do

            !---------------------!
            ! Orthogonalize to Xb !
            !---------------------!
            do i = 1, nx2
                ad2 = nr1 + i
                WV = WV - X(:,ad2)*B(i,ad)
            end do

            !---------------------!
            ! Orthogonalize to Xa !
            !---------------------!
            do i = 1, nr1
                WV = WV - MX(:,i)*A(i,ad)
            end do

            !----------------------------------!
            ! Orthogonalize to turning vectors !
            !----------------------------------!
            WN = X(:,cc)

            do i = nt, j+1, -1
                ad2 = nxmnt + i
                val = 0.0D0
                do k = 1, nn
                    val = val + X(k,ad2)*WN(k)
                end do

                WV = WV - X(:,ad2)*val
            end do

            normt = 0.0D0
            do k = 1, nn
                normt = normt + WV(k)*WN(k)
            end do

            if (normt <= 0.0D0) then
                stop 'NORMT IN ESUBSPACE_PROJ SHOULD BE .GT. ZERO'
            end if

            X(:,cc) = WV / SQRT(normt)
        end do

        !--------------------------------------------!
        ! Project stiffness matrix onto the subspace !
        !--------------------------------------------!
        nq2 = nq*2
        np1 = np + 1

        KR = 0.0D0

        !--------------------!
        ! First block column !
        !--------------------!
        ct = 0
        do j = 1, np
            ad = j + (nq2-j)*(j-1)/2
            ct = ct + np1 - j
            KR(ad) = 1.0D0 / LAMBDA(j)

            ct = ct + nx
        end do

        !--------------------!
        ! Third block column !
        !--------------------!
        ad = np + nr1
        ctt = nq + (nq2-ad)*(ad-1)/2
        do j = 1, nr2
            !----------------!
            ! WN = KK^(-1)*X !
            !----------------!
            cc = nr1+j
            WN = X(:,cc)
            call REDBAK(KK,WN,DK,nn)

            !----------------!
            ! Compute WN^T*X !
            !----------------!
            do i = j, nr2
                ctt = ctt + 1
                ad  = nr1 + i
                do k = 1, nn
                    KR(ctt) = KR(ctt) + X(k,ad)*WN(k)
                end do
            end do

            !--------------------------!
            ! Update iteration vectors !
            !--------------------------!
            X(:,cc) = WN
        end do

        !---------------------!
        ! Second block column !
        !---------------------!
        do j = 1, nr1
            !-------------!
            ! Copy from A
            !-------------!
            sct = ct + 1
            ct  = sct + nr1 - j
            KR(sct:ct) = A(j:nr1,j)

            do i = 1, nr2
                ct = ct + 1
                if (i <= nx2) then
                    !-------------!
                    ! Copy from B !
                    !-------------!
                    KR(ct) = B(i,j)
                else
                    !----------------!
                    ! Compute X^T*MX !
                    !----------------!
                    ad = nr1 + i
                    do k = 1, nn
                        KR(ct) = KR(ct) + X(k,ad)*MX(k,j)
                    end do
                end if
            end do
        end do

        !---------------------------------------!
        ! Project mass matrix onto the subspace !
        !---------------------------------------!
        MR = 0.0D0
        ct = 0

        !--------------------!
        ! First block column !
        !--------------------!
        do j = 1, np
            ad = j + (nq2-j)*(j-1)/2
            ct = ct + np1 - j
            MR(ad) = 1.0D0 / (LAMBDA(j)*LAMBDA(j))

            !-------------!
            ! Copy from D !
            !-------------!
            do i = 1, nr1
                ct = ct + 1
                MR(ct) = D(j,i) / LAMBDA(j)
            end do

            !-----------------!
            ! Compute X^T*PHI !
            !-----------------!
            do i = 1, nr2
                ct = ct + 1
                ad = nr1 + i
                do k = 1, nn
                    MR(ct) = MR(ct) + X(k,ad)*PHI(k,j)
                end do

                MR(ct) = MR(ct) / LAMBDA(j)
            end do
        end do

        !---------------------!
        ! Second block column !
        !---------------------!
        do j = 1, nr1
            !-------------!
            ! Copy from C !
            !-------------!
            sct = ct + 1
            ct  = sct + nr1 - j
            MR(sct:ct) = C(j:nr1,j)

            !---------------!
            ! Compute X^T*X !
            !---------------!
            do i = 1, nr2
                ct = ct + 1
                ad = nr1 + i
                do k = 1, nn
                    MR(ct) = MR(ct) + X(k,ad)*X(k,j)
                end do
            end do
        end do

        !--------------------!
        ! Third block column !
        !--------------------!
        do j = 1, nr2
            cc = nr1 + j

            !-----------!
            ! WN = MM*X !
            !-----------!
            call MULT(WN,MM,X(:,cc),DM,nn,nwm)

            !----------------!
            ! Compute X^T*WN !
            !----------------!
            do i = j, nr2
                ct = ct + 1
                ad = nr1 + i
                do k = 1, nn
                    MR(ct) = MR(ct) + X(k,ad)*WN(k)
                end do
            end do

            !--------------------------!
            ! Update iteration vectors !
            !--------------------------!
            X(:,cc) = WN
        end do
    END SUBROUTINE ESUBSPACE_PROJ

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Check whether the turning vector is null
    !
    ! INPUT:
    !   A(:,:)   = Xa^T * M * BX
    !   B(:,:)   = Xb^T * M * BX
    !   C(:,:)   = BX^T * M * BX
    !   D(:,:)   = PHI^T * M * BX
    !   rr       = index of the turing vector for the check
    !   nt       = number of previously calculated turning vectors
    !   tol      = tolerance for selecting the turning vector
    !   TID(:)   = indices of selected turning vectors
    !   ALPHA(:) = auxiliary array
    !
    ! OUTPUT:
    !   nt
    !   TID(:)
    !   ALPHA(:)
    !   iftv   = flag for the selection;
    !            EQ. 0 for discard
    !            EQ. 1 for use
    !   RIJ(:) = working array
    !---------------------------------------------------------------------------
    SUBROUTINE CHECKNULL(A,B,C,D,rr,nt,TID,ALPHA,tol,iftv,RIJ)
        integer, intent(in)    :: rr
        integer, intent(inout) :: nt, TID(:)
        integer, intent(out)   :: iftv
        double precision, intent(in)    :: A(:,:), B(:,:), C(:,:), D(:,:), tol
        double precision, intent(inout) :: ALPHA(:)
        double precision, intent(out)   :: RIJ(:)

        integer :: j, k, ii, cc
        integer :: tp
        double precision :: val
        double precision :: rik, rjk, rkk

        iftv = 0
        ii   = nt + 1
        do j = 1, ii
            if (j == ii) then
                cc = rr
            else
                cc = TID(j)
            end if

            RIJ(j) = C(rr,cc) - DOT_PRODUCT(A(:,rr),A(:,cc))&
                - DOT_PRODUCT(B(:,rr),B(:,cc))&
                - DOT_PRODUCT(D(:,rr),D(:,cc))

            tp = (j-1)*j/2
            do k = 1, j-1
                rik = RIJ(k)
                if (j == ii) then
                    rjk = RIJ(k)
                else
                    rjk = ALPHA(k + tp)
                end if

                rkk = ALPHA( k*(k+1)/2 )

                if (rkk == 0.0D0) then
                    stop 'RKK IN CHECKNULL SHOULD NOT BE .EQ. ZERO'
                end if

                RIJ(j) = RIJ(j) - rik*rjk/rkk
            end do
        end do

        if (C(rr,rr) == 0.0D0) then
            stop 'C(RR,RR) IN CHECKNULL SHOULD NOT BE .EQ. ZERO'
        end if

        val = RIJ(ii) / C(rr,rr)
        if (val > tol) then
            nt      = ii
            TID(ii) = rr

            tp = (ii-1)*ii/2
            do j = 1, ii
                ALPHA(j + tp) = RIJ(j)
            end do

            iftv = 1
        end if
    END SUBROUTINE CHECKNULL

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Rearrange the calculated eigenpairs and check whether the smallest nroot 
    !   eigenpairs converged to the tolerance
    !
    ! NEED:
    !   INSERTION_SORT
    !
    ! INPUT:
    !   nq        = number of iteration vectors
    !   nroot     = number of required eigenpairs
    !   EIGV(nq)  = calculated eigenvalues
    !   QQ(nq,nq) = Ritz coordinates
    !   RTOLV(nq) = calculated error bounds
    !   np        = number of previously converged eigenpairs
    !   tolc      = tolerance for convergence
    !
    ! OUTPUT:
    !   iconv      = flag for convergence;
    !                EQ. 0 no convergence
    !                EQ. 1 indicates that the smallest nroot eigenpairs 
    !                      converged to the tolerance
    !   EIGV(nq)   = rearranged eigenvalues
    !   QQ(nq,nq)  = rearranged Ritz coordinates
    !   RTOLV(nq)  = rearranged error bounds
    !   np         = number of converged eigenpairs
    !   TEIGV(nq)  = eigenvalues in ascending order
    !   TRTOLV(nq) = error bounds corresponding to TEIGV
    !   TVEC(nq)   = working array
    !---------------------------------------------------------------------------
    SUBROUTINE DETERMINE_CONV(nq,nroot,EIGV,QQ,RTOLV,np,tolc,iconv,TEIGV,TRTOLV,&
            TVEC)
        integer, intent(in)    :: nq, nroot
        integer, intent(inout) :: np
        integer, intent(out)   :: iconv
        double precision, intent(in)    :: tolc
        double precision, intent(inout) :: EIGV(nq), QQ(nq,nq), RTOLV(nq)
        double precision, intent(out)   :: TEIGV(nq), TRTOLV(nq), TVEC(nq)

        integer :: j, add
        double precision :: tval1, tval2

        !----------------------------------------------------------------------!
        ! Check for convergence of iteration vectors and rearrange the results !
        !----------------------------------------------------------------------!
        add =  np
        do j = np+1, nq
            if (RTOLV(j) <= tolc) then
                add = add + 1
                np  = np + 1
                if (j > add) then
                    tval1 = EIGV(j)
                    tval2 = RTOLV(j)
                    TVEC  = QQ(:,j)

                    EIGV(j)  = EIGV(add)
                    RTOLV(j) = RTOLV(add)
                    QQ(:,j)  = QQ(:,add)

                    EIGV(add)  = tval1
                    RTOLV(add) = tval2
                    QQ(:,add)  = TVEC
                end if
            end if
        end do

        !-------------------------------------------------------!
        ! Check whether the smallest nroot eigenpairs converged !
        !-------------------------------------------------------!
        iconv  = 0
        TEIGV  = EIGV
        TRTOLV = RTOLV

        call INSERTION_SORT(TEIGV,TRTOLV,nq,1,TVEC(1))

        do j = 1, nroot
            if (TRTOLV(j) > tolc) then
                exit
            end if
        end do

        j = j - 1
        if (j == nroot) then
            iconv = 1
            EIGV = TEIGV
            TRTOLV = RTOLV
        end if
    END SUBROUTINE DETERMINE_CONV

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Perform single basic subspace iteration
    !
    ! NEED:
    !   REDBAK
    !   MULT
    !   JACOBI
    !   INSERTION_SORT
    !
    ! INPUT:
    !   KK(nwk)   = factorized stiffness matrix (L*D*L^T by COLSOL) in skyline
    !               format
    !   MM(nwm)   = mass matrix in skyline format
    !   DK(nn1)   = address of diagonal elements of factorized stiffness matrix
    !   DM(nn1)   = address of diagonal elements of mass matrix
    !   nn        = order of system
    !   nn1       = nn + 1
    !   nwk       = number of elements below skyline of stiffness matrix
    !   nwm       = number of elements below skyline of mass matrix;
    !               EQ. nwk for consistent mass matrix
    !               EQ. nn for lumped mass matrix
    !   nq        = number of iteration vectors
    !   nnq       = nq*(nq + 1)/2
    !   RV(nn,nq) = iteration vectors
    !   ifrv      = flag for the form of RV;
    !               EQ. 0 for RV in the form of M*X
    !               EQ. 1 for RV in the form of X
    !
    ! OUTPUT:
    !   EIGV(nq)  = calculated eigenvalues
    !   RV(nn,nq) = M-orthonormal iteration vectors
    !   RTOLV(nq) = calculated error bounds
    !   KR(nnq)   = working array
    !   MR(nnq)   = working array
    !   WN(nn)    = working array
    !   QQ(nq,nq) = working array
    !   WQ(nq)    = working array
    !---------------------------------------------------------------------------
    SUBROUTINE BSUBSPACE_ITER1(KK,MM,DK,DM,nn,nn1,nwk,nwm,nq,nnq,RV,EIGV,RTOLV,&
            ifrv,KR,MR,WN,QQ,WQ)
        integer, intent(in) :: nn, nn1, nwk, nwm, nq, nnq, ifrv
        integer, intent(in) :: DK(nn1), DM(nn1)
        double precision, intent(in)    :: KK(nwk), MM(nwm)
        double precision, intent(inout) :: RV(nn,nq)
        double precision, intent(out)   :: EIGV(nq), RTOLV(nq), KR(nnq),&
            MR(nnq), WN(nn), QQ(nq,nq), WQ(nq)

        integer, parameter :: ifpr = 0, nsmax = 12
        double precision, parameter :: tolj = 1.0D-12, tol24 = 1.0D-24

        integer :: i, j, k, ct
        double precision :: vdot, val

        !------------------------------!
        ! Project KK onto the subspace !
        !------------------------------!
        ct = 0
        KR = 0.0D0
        do j = 1, nq
            !-----------------!
            ! WN = KK^(-1)*RV !
            !-----------------!
            WN = RV(:,j)
            call REDBAK(KK,WN,DK,nn) 
            do i = j, nq
                ct = ct + 1
                do k = 1, nn
                    KR(ct) = KR(ct) + RV(k,i)*WN(k)
                end do
            end do

            !--------------------------!
            ! Update iteration vectors !
            !--------------------------!
            RV(:,j) = WN
        end do

        !------------------------------!
        ! Project MM onto the subspace !
        !------------------------------!
        ct = 0
        MR = 0.0D0
        do j = 1, nq
            !------------!
            ! WN = MM*RV !
            !------------!
            call MULT(WN,MM,RV(:,j),DM,nn,nwm)
            do i = j, nq
                ct = ct + 1
                do k = 1, nn
                    MR(ct) = MR(ct) + RV(k,i)*WN(k)
                end do
            end do

            !--------------------------!
            ! Update iteration vectors !
            !--------------------------!
        if (ifrv .ne. 1) then
                RV(:,j) = WN
            end if
        end do

        !-----------------------------------------------------!
        ! Solve for the eigensystem of the projected matrices !
        !-----------------------------------------------------!
        call JACOBI(KR,MR,QQ,EIGV,WQ,nq,nnq,tolj,nsmax,ifpr,0)

        !-----------------------------------------------!
        ! Arrange calculated results in ascending order !
        !-----------------------------------------------!
        call INSERTION_SORT(EIGV,QQ,nq,nq,WQ)

        !---------------------------!
        ! Improve iteration vectors !
        !---------------------------!
        do i = 1, nn
            WQ = RV(i,:)
            do j = 1, nq
                RV(i,j) = 0.0D0
                do k = 1, nq
                    RV(i,j) = RV(i,j) + WQ(k)*QQ(k,j)
                end do
            end do
        end do

        !------------------------!
        ! Calculate error bounds !
        !------------------------!
        do j = 1, nq
            vdot = 0.0D0
            do k = 1, nq
                vdot = vdot + QQ(k,j)*QQ(k,j)
            end do

            if (vdot == 0.0D0) then
                stop 'VDOT IN BSUBSPACE_ITER1 SHOULD NOT BE .EQ. ZERO'
            end if

            val = vdot - EIGV(j)*EIGV(j)
            val = val / vdot
            val = MAX(val,tol24)
            RTOLV(j) = SQRT(val)
        end do
    END SUBROUTINE BSUBSPACE_ITER1

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Arrange the calculated eigenvalues and corresponding Ritz coordinates 
    !   in ascending order
    !
    ! INPUT:
    !   X(n)   = calculated eigenvalues
    !   A(m,n) = Ritz coordinates
    !   n      = number of calculated eigenvalues
    !   m      = first dimension of A
    !
    ! OUTPUT:
    !   X(n)   = arranged eigenvalues
    !   A(m,n) = arranged Ritz coordinates
    !   TT(m)  = working array
    !---------------------------------------------------------------------------
    SUBROUTINE INSERTION_SORT(X,A,n,m,TT)
        integer, intent(in) :: n, m
        double precision, intent(inout) :: X(n), A(m,n)
        double precision, intent(out)   :: TT(m)

        integer :: i, j, j1
        integer :: ti
        double precision :: tv

        do i = 2, n
            tv = X(i)
            TT = A(:,i)

            j = i - 1
            do
                if (j < 1) then
                    exit
                else if (X(j) <= tv) then
                    exit
                end if

                j1      = j + 1
                X(j1)   = X(j)
                A(:,j1) = A(:,j)
                j       = j - 1
            end do

            j1      = j + 1
            X(j1)   = tv
            A(:,j1) = TT
        end do
    END SUBROUTINE INSERTION_SORT

    ! ************************
    ! * Subroutines to print *
    ! ************************
    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Print heading
    !
    ! INPUT:
    !   n = unit specifier in write statement
    !---------------------------------------------------------------------------
    SUBROUTINE PRINT_HEADING(n)
        integer, intent(in) :: n

        character(10) :: date(3)
        integer :: dateval(8)

        call DATE_AND_TIME( date(1), date(2), date(3), dateval )

        write(n, '( 8X, 66("*") )')
        write(n, '( 8X, "*", 64X, "*" )')
        write(n, '( 8X, "* SOLVE GENERALIZED EIGENPROBLEMS BY ENRICHED", X,&
            "SUBSPACE ITERATION *")')
        write(n, '( 8X, "*", 64X, "*" )')
        write(n, '( 8X, "*", 64X, "*" )')
        write(n, '( 8X, "* AUTHORS    : KI-TAE KIM,", 39X, "*" )')
        write(n, '( 8X, "*              KLAUS-JÜRGEN BATHE", 32X, "*" )')
        write(n, '( 8X, "* AFFILIATION:", X,&
            "DEPARTMENT OF MECHANICAL ENGINEERING,", 13X, "*" )')
        write(n, '( 8X, "*", 14X, "MASSACHUSETTS INSTITUTE OF TECHNOLOGY", 13X,&
            "*" )')
        write(n, '( 8X, "* E-MAIL     : ktkim@mit.edu or qlsn55@gmail.com",&
            17X, "*" )')
        write(n, '( 8X, "*", 64X, "*" )')
        write(n, '( 8X, "* TIME OF ANALYSIS : ", I2, ":", I2, ", ", I2,& 
            "/", I2, "/", I4, 26X, " *" )' ) dateval(5), dateval(6),&
            dateval(2), dateval(3), dateval(1)
        write(n, '( 8X, "*", 64X, "*" )')
        write(n, '( 8X, 66("*") )')
    END SUBROUTINE PRINT_HEADING

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Print input values
    !
    ! INPUT:
    !   n       = unit specifier in write statement
    !   nn      = order of system
    !   nroot   = number of required eigenpairs
    !   nq      = number of iteration vectors
    !   tolc    = tolerance for convergence;
    !   tolt    = tolerance for selecting the turning vectors
    !   KK(:)   = factorized stiffness matrix (L*D*L^T by COLSOL) in skyline
    !             format
    !   MM(:)   = mass matrix in skyline format
    !   DK(:)   = address of diagonal elements of factorized stiffness matrix
    !   DM(:)   = address of diagonal elements of mass matrix
    !   RV(:,:) = iteration vectors
    !---------------------------------------------------------------------------
    SUBROUTINE PRINT_SETTING(KK,MM,DK,DM,RV,nn,nroot,nq,tolc,tolt,n)
        integer, intent(in) :: nn, nroot, nq, n
        integer, intent(in) :: DK(:), DM(:)
        double precision, intent(in) :: tolc, tolt
        double precision, intent(in) :: KK(:), MM(:), RV(:,:)

        integer :: msize, wsize

        msize = (SIZE(KK) + SIZE(MM))*8 + (SIZE(DK) + SIZE(DM))*4
        wsize = (SIZE(RV,1)*SIZE(RV,2))*8
        msize = msize/1000
        wsize = wsize/1000

        write(n, '(/," TOTAL DEGREES OF FREEDOM                =", 2X, I8)') nn
        write(n, '(" NUMBER OF REQUIRED EIGENPAIRS           =", 2X, I8)')&
            nroot
        write(n, '(" NUMBER OF ITERATION VECTORS USED        =", 2X, I8)') nq
        write(n, '(" TOLERANCE FOR SELECTING TURNING VECTORS =", 2X, ES11.4)')&
           tolt
        write(n, '(" TOLERANCE FOR CONVERGENCE               =", 2X, ES11.4)')&
            tolc
        write(n, '(/," ALLOCATED MEMORY FOR")')
        write(n, '(4X, "MASS & STIFFNESS ARRAYS =", 2X, I8, X, "KB")') msize
        write(n, '(4X, "WORKING ARRAY ", 10X, "=" 2X, I8, X, "KB")') wsize
    END SUBROUTINE PRINT_SETTING
    !===========================================================================
END MODULE ESSPACE
