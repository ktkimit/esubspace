!===============================================================================
! P R O G R A M
!   Example program for the use of ESUBSPACE_ITER
!
!   References:
!   [1] K.J. Bathe, Finite element procedures, 2nd ed., 2014
!
! NEED:
!   ESUBSPACE_ITER
!
! AUTHORS    : Ki-Tae Kim
!              Klaus-Jürgen Bathe
! EMAIL      : ktkim@mit.edu or qlsn55@gmail.com
! AFFILIATION: Massachusetts Institute of Technology
! LANGUAGE   : FORTRAN 90
! CREATED    : Fri 21 Apr 2017 12:36:28 AM EDT
!===============================================================================
PROGRAM ESSPACE_MAIN
    USE ESSPACE, ONLY : ESUBSPACE_ITER

    IMPLICIT NONE

    character(20) :: FNAME

    integer              :: nn, nwk, nwm, nroot, nitm, nq, nnq, ifss, ilog,&
        nstiff
    integer, allocatable :: DK(:), DM(:)
    double precision              :: tolc, tolt
    double precision, allocatable :: KK(:), MM(:), EIGV(:), RV(:,:)

    call MAIN_TEST

CONTAINS
    !===========================================================================
    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Example subroutine to show how to run ESUBSPACE_ITER
    !---------------------------------------------------------------------------
    SUBROUTINE MAIN_TEST()
        integer :: dum

        call PRINT_HEADING()

        ! Read stiffness matrix and mass matrix data and input values from files
        call READ_INPUTDATA()

        ! Check the exception
        call CHECK_EXCEPTION()

        ! Set number of iteration vectors used
        nq  = nroot*2
        nnq = nq*(nq + 1)/2

        allocate( EIGV(nq), RV(nn,nq+nroot) )

        ! Generate an output formatted log file
        open(unit=ilog, file="frame2d_ex.log", status='replace',&
            form='formatted',action='write', position='rewind')

        ! Open a scratch file for stiffness matrix
        nstiff = ilog + 10
        open(unit=nstiff, status='scratch')

        ! Run the subroutine of the enriched subspace iteration method
        call ESUBSPACE_ITER(KK,MM,DK,DM,nn,nwk,nwm,nroot,nitm,nq,nnq,tolc,&
            tolt,ifss,ilog,nstiff,EIGV,RV)

        ! Deallocate arrays
        deallocate( KK, MM, DK, DM, EIGV, RV )

        ! Close file
        close(unit=ilog)
        close(unit=nstiff)
    END SUBROUTINE MAIN_TEST

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Check the exception
    !---------------------------------------------------------------------------
    SUBROUTINE CHECK_EXCEPTION()
        if (nn < nroot*2) then
            write(*,101)
            write(*,102)
            stop
        else if (nn <= 1) then
            write(*,101)
            write(*,103)
            stop
        end if

        101 format(//, X, "*** ERROR *** PROGRAM STOP")
        102 format(/, 4X, "'NUMBER OF EIGENPAIRS SOUGHT*2' IS GREATER THAN", & 
            " 'THE SYSTEM SIZE'")
        103 format(/, 4X, "'THE SYSTEM SIZE' IS LESS THAN OR EQUAL TO 1")
    END SUBROUTINE CHECK_EXCEPTION

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Read input data
    !
    ! NEED:
    !   *.stiff
    !   *.mass
    !   *.in
    !
    ! ALLOCATED VARIABLES:
    !   DK, KK
    !   DM, MM
    !---------------------------------------------------------------------------
    SUBROUTINE READ_INPUTDATA()
        integer :: dumi
        character(50) :: dumc

        !-----------------!
        ! Input file name !
        !-----------------!
        write(*, '(4X, "Input the file name (without extension)")')
        write(*, '(8X, ">")', advance='no')
        read(*,*) FNAME

        !------------!
        ! Open files !
        !------------!
        open(unit=1, file=TRIM(FNAME)//".stiff", status='old',&
            action='read', form='unformatted', access='stream')
        open(unit=2, file=TRIM(FNAME)//".mass", status='old',&
            action='read', form='unformatted', access='stream')
        open(unit=3, file=TRIM(FNAME)//".in", status='old',&
            action='read', form='formatted')

        ! Read stiffness matrix data
        read(1) nn
        read(1) nwk
        allocate( DK(nn+1), KK(nwk) )
        read(1) DK
        read(1) KK

        ! Read mass matrix data
        read(2) dumi
        read(2) nwm
        allocate( DM(nn+1), MM(nwm) )
        read(2) DM
        read(2) MM

        ! Read input values from *.in
        do
            read(3,'(A31)') dumc(1:31)
            if (dumc(1:31) == "* NUMBER OF EIGENPAIRS SOUGHT *") then
                exit
            end if
        end do
        read(3,*) nroot

        do
            read(3,'(A42)') dumc(1:42)
            if (dumc(1:42) == "* MAXIMUM NUMBER OF ITERATIONS PERMITTED *") then
                exit
            end if
        end do
        read(3,*) nitm

        do
            read(3,'(A40)') dumc(1:40)
            if (dumc(1:40) == "* CONVERGENCE TOLERANCE ON EIGENVALUES *") then
                exit
            end if
        end do
        read(3,*) tolc

        do
            read(3,'(A43)') dumc(1:43)
            if (dumc(1:43) == "* TOLERANCE FOR SELECTING TURNING VECTORS *")&
                then
                exit
            end if
        end do
        read(3,*) tolt

        do
            read(3,'(A24)') dumc(1:24)
            if (dumc(1:24) == "* STURM SEQUENCE CHECK *") then
                exit
            end if
        end do
        read(3,*)
        read(3,*)
        read(3,*) ifss

        do
            read(3,'(A46)') dumc(1:46)
            if (dumc(1:46) == "* PRINTING LOG IN FORMATTED FILE             *")&
                then
                exit
            end if
        end do
        read(3,*)
        read(3,*)
        read(3,*) ilog

        !-------------!
        ! Close files !
        !-------------!
        close(unit=1)
        close(unit=2)
        close(unit=3)
    END SUBROUTINE READ_INPUTDATA

    !---------------------------------------------------------------------------
    ! S U B R O U T I N E
    !   Print heading
    !---------------------------------------------------------------------------
    SUBROUTINE PRINT_HEADING()
        character(10) :: date(3)
        integer :: dateval(8)

        call DATE_AND_TIME( date(1), date(2), date(3), dateval )

        write(*, '( 8X, 66("*") )')
        write(*, '( 8X, "*", 64X, "*" )')
        write(*, '( 8X, "* SOLVE GENERALIZED EIGENPROBLEMS BY ENRICHED", X,&
            "SUBSPACE ITERATION *")')
        write(*, '( 8X, "*", 64X, "*" )')
        write(*, '( 8X, "*", 64X, "*" )')
        write(*, '( 8X, "* AUTHORS    : KI-TAE KIM,", 39X, "*" )')
        write(*, '( 8X, "*              KLAUS-JÜRGEN BATHE", 32X, "*" )')
        write(*, '( 8X, "* AFFILIATION:", X,&
            "DEPARTMENT OF MECHANICAL ENGINEERING,", 13X, "*" )')
        write(*, '( 8X, "*", 14X, "MASSACHUSETTS INSTITUTE OF TECHNOLOGY", 13X,&
            "*" )')
        write(*, '( 8X, "* E-MAIL     : ktkim@mit.edu or qlsn55@gmail.com",&
            17X, "*" )')
        write(*, '( 8X, "*", 64X, "*" )')
        write(*, '( 8X, "* TIME OF ANALYSIS : ", I2, ":", I2, ", ", I2,& 
            "/", I2, "/", I4, 26X, " *" )' ) dateval(5), dateval(6),&
            dateval(2), dateval(3), dateval(1)
        write(*, '( 8X, "*", 64X, "*" )')
        write(*, '( 8X, 66("*") )')
    END SUBROUTINE PRINT_HEADING
    !===========================================================================
END PROGRAM ESSPACE_MAIN
