!
!  L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License"
!  or "3-clause license")
!  Please read attached file License.txt
!
!                             DRIVER 3  in Fortran 90
!     --------------------------------------------------------------
!            TIME-CONTROLLED DRIVER FOR L-BFGS-B
!     --------------------------------------------------------------
!
!        L-BFGS-B is a code for solving large nonlinear optimization
!             problems with simple bounds on the variables.
!
!        The code can also be used for unconstrained problems and is
!        as efficient for these problems as the earlier limited memory
!                          code L-BFGS.
!
!        This driver shows how to terminate a run after some prescribed
!        CPU time has elapsed, and how to print the desired information
!        before exiting.
!
!     References:
!
!        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!        memory algorithm for bound constrained optimization'',
!        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!        Subroutines for Large Scale Bound Constrained Optimization''
!        Tech. Report, NAM-11, EECS Department, Northwestern University,
!        1994.
!
!
!          (Postscript files of these papers are available via anonymous
!           ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                              *  *  *
!
!         February 2011   (latest revision)
!         Optimization Center at Northwestern University
!         Instituto Tecnologico Autonomo de Mexico
!
!         Jorge Nocedal and Jose Luis Morales

      program driver3

      use lbfgsb_module, wp => lbfgsp_wp
      use iso_fortran_env, only: output_unit

      ! This time-controlled driver shows that it is possible to terminate
      ! a run by elapsed CPU time, and yet be able to print all desired
      ! information. This driver also illustrates the use of two
      ! stopping criteria that may be used in conjunction with a limit
      ! on execution time. The sample problem used here is the same as in
      ! driver1 and driver2 (the extended Rosenbrock function with bounds
      ! on the variables).

      implicit none

      ! We specify a limit on the CPU time (tlimit = 10 seconds)
      !
      ! We suppress the default output (iprint = -1). The user could
      !   also elect to use the default output by choosing iprint >= 0.)
      ! We suppress the code-supplied stopping tests because we will
      !   provide our own termination conditions
      ! We specify the dimension n of the sample problem and the number
      !    m of limited memory corrections stored.

      integer,  parameter    :: n = 1000, m = 10, iprint = -1
      real(wp), parameter    :: factr  = 0.0_wp, pgtol  = 0.0_wp, &
                                tlimit = 10.0_wp
!
      character(len=60)      :: task, csave
      logical                :: lsave(4)
      integer                :: isave(44)
      real(wp)               :: f
      real(wp)               :: dsave(29)
      integer,  allocatable  :: nbd(:), iwa(:)
      real(wp), allocatable  :: x(:), l(:), u(:), g(:), wa(:)
!
      real(wp) :: t1, t2, time1, time2
      integer  :: i, j

      allocate ( nbd(n), x(n), l(n), u(n), g(n) )
      allocate ( iwa(3*n) )
      allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )

      ! This time-controlled driver shows that it is possible to terminate
      ! a run by elapsed CPU time, and yet be able to print all desired
      ! information. This driver also illustrates the use of two
      ! stopping criteria that may be used in conjunction with a limit
      ! on execution time. The sample problem used here is the same as in
      ! driver1 and driver2 (the extended Rosenbrock function with bounds
      ! on the variables).

      ! We now specify nbd which defines the bounds on the variables:
      !                l   specifies the lower bounds,
      !                u   specifies the upper bounds.

      ! First set bounds on the odd-numbered variables.

      do i=1, n,2
         nbd(i)=2
         l(i)=1.0_wp
         u(i)=100.0_wp
      end do

      ! Next set bounds on the even-numbered variables.

      do i=2, n,2
         nbd(i)=2
         l(i)=-100.0_wp
         u(i)=100.0_wp
      end do

      ! We now define the starting point.

      do i=1, n
         x(i)=3.0_wp
      end do

      ! We now write the heading of the output.

      write (output_unit,'(/,5x, a,/,5x, a,/)') &
               'Solving sample problem.',&
               ' (f = 0.0 at the optimal solution.)'

      ! We start the iteration by initializing task.

      task = 'START'

      ! ------- the beginning of the loop ----------

      ! We begin counting the CPU time.

      call cpu_time(time1)

      do while( task(1:2)=='FG'.or.task=='NEW_X'.or.task=='START')

         ! This is the call to the L-BFGS-B code.

         call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa, &
                     task,iprint, csave,lsave,isave,dsave)

         if (task(1:2) == 'FG') then

         ! the minimization routine has returned to request the
         ! function f and gradient g values at the current x.
         ! Before evaluating f and g we check the CPU time spent.

         call cpu_time(time2)
         if (time2-time1 > tlimit) then

            task='STOP: CPU EXCEEDING THE TIME LIMIT.'

            ! Note: Assigning task(1:4)='STOP' will terminate the run;
            ! setting task(7:9)='CPU' will restore the information at
            ! the latest iterate generated by the code so that it can
            ! be correctly printed by the driver.

            ! In this driver we have chosen to disable the
            ! printing options of the code (we set iprint=-1);
            ! instead we are using customized output: we print the
            ! latest value of x, the corresponding function value f and
            ! the norm of the projected gradient |proj g|.

            ! We print out the information contained in task.

              write (output_unit,*) task

            ! We print the latest iterate contained in wa(j+1:j+n), where

              j = 3*n+2*m*n+11*m**2
              write (output_unit,*) 'Latest iterate X ='
              write (output_unit,'((1x,1p, 6(1x,d11.4)))') (wa(i),i = j+1,j+n)

            ! We print the function value f and the norm of the projected
            ! gradient |proj g| at the last iterate; they are stored in
            ! dsave(2) and dsave(13) respectively.

              write (output_unit,'(a,1p,d12.5,4x,a,1p,d12.5)') &
              'At latest iterate   f =',dsave(2),'|proj g| =',dsave(13)

         else

              ! The time limit has not been reached and we compute
              ! the function value f for the sample problem.

              f=0.25_wp*(x(1)-1.0_wp)**2
              do i=2, n
                 f=f+(x(i)-x(i-1)**2)**2
              end do
              f=4.0_wp*f

            ! Compute gradient g for the sample problem.

               t1 = x(2) - x(1)**2
               g(1) = 2.0_wp*(x(1)-1.0_wp)-16.0_wp*x(1)*t1
               do i=2,n-1
                  t2=t1
                  t1=x(i+1)-x(i)**2
                  g(i)=8.0_wp*t2-16.0_wp*x(i)*t1
               end do
               g(n)=8.0_wp*t1
            endif

            ! go back to the minimization routine.
      else

           if (task(1:5) == 'NEW_X') then

            ! the minimization routine has returned with a new iterate.
            ! The time limit has not been reached, and we test whether
            ! the following two stopping tests are satisfied:

            ! 1) Terminate if the total number of f and g evaluations
            ! exceeds 900.

            if (isave(34) >= 900) &
            task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'

            ! 2) Terminate if  |proj g|/(1+|f|) < 1.0d-10.

            if (dsave(13) <= 1.0e-10_wp*(1.0_wp + abs(f))) &
            task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'

            ! We wish to print the following information at each iteration:
            !   1) the current iteration number, isave(30),
            !   2) the total number of f and g evaluations, isave(34),
            !   3) the value of the objective function f,
            !   4) the norm of the projected gradient,  dsve(13)
            !
            ! See the comments at the end of driver1 for a description
            ! of the variables isave and dsave.

            write (output_unit,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate' &
                   ,isave(30),'nfg =',isave(34),'f =',f,'|proj g| =',dsave(13)

            ! If the run is to be terminated, we print also the information
            ! contained in task as well as the final value of x.

            if (task(1:4) == 'STOP') then
               write (output_unit,*) task
               write (output_unit,*) 'Final X='
               write (output_unit,'((1x,1p, 6(1x,d11.4)))') (x(i),i = 1,n)
            endif

          endif
        end if
      end do

      ! If task is neither FG nor NEW_X we terminate execution.

      end program driver3