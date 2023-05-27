!*******************************************************************************
!>
!  This is a modernized version L-BFGS-B.
!
!  It is based on the the modified version of L-BFGS-B described
!  in the following paper:
!
!  * Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778:
!    L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
!    Optimization"  (2011). ACM Transactions on Mathematical Software,
!    Volume 38, Issue 1, Article No.: 7, pp 1-4
!
!  The paper describes an improvement and a correction to Algorithm 778.
!  It is shown that the performance of the algorithm can be improved
!  significantly by making a relatively simple modication to the subspace
!  minimization phase. The correction concerns an error caused by the use
!  of routine dpmeps to estimate machine precision.
!
!### License
!  L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License"
!  or "3-clause license")
!  Please read attached file License.txt
!
!### Version
!  * Based on: L-BFGS-B (version 3.0).  April 25, 2011
!  * Refactored and modernized by Jacob Williams, 2023
!
!### Original Authors
!  * J. Nocedal  Department of Electrical Engineering and
!    Computer Science.
!    Northwestern University. Evanston, IL. USA
!  * J.L Morales  Departamento de Matematicas,
!    Instituto Tecnologico Autonomo de Mexico
!    Mexico D.F. Mexico.

    module lbfgsb_module

    use iso_fortran_env, only: wp => real64 ! for now
    use linpack
    use blas

    implicit none

    private

    ! constants:
    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: one  = 1.0_wp
    real(wp),parameter :: two   = 2.0_wp
    real(wp),parameter :: three = 3.0_wp

    public :: setulb

    contains
!*******************************************************************************

    !TODO make a high-level wrapper so the user doesn't have to call the
    !     reverse communication routine directly.

!*******************************************************************************
!>
!  The main routine.
!
!  This subroutine partitions the working arrays `wa` and `iwa`, and
!  then uses the limited memory BFGS method to solve the bound
!  constrained optimization problem by calling [[mainlb]].
!  (The direct method will be used in the subspace minimization.)
!
!### References
!  1. R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, "A limited
!     memory algorithm for bound constrained optimization",
!     SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!  2. C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, "L-BFGS-B: a
!     limited memory FORTRAN code for solving bound constrained
!     optimization problems", Tech. Report, NAM-11, EECS Department,
!     Northwestern University, 1994.
!
!### History
!  * NEOS, November 1994. (Latest revision June 1996.)
!    Optimization Technology Center.
!    Argonne National Laboratory and Northwestern University.
!    Written by Ciyou Zhu in collaboration with
!    R.H. Byrd, P. Lu-Chen and J. Nocedal.

      subroutine setulb(n,m,x,l,u,Nbd,f,g,Factr,Pgtol,Wa,Iwa,Task, &
                        Iprint,Csave,Lsave,Isave,Dsave)
      implicit none

      integer,intent(in) :: n !! the dimension of the problem (the number of variables).
      integer,intent(in) :: m !! the maximum number of variable metric corrections
                              !! used to define the limited memory matrix.
                              !! Values of `m < 3`  are
                              !! not recommended, and large values of `m` can result in excessive
                              !! computing time. The range  `3 <= m <= 20` is recommended.
      real(wp),intent(inout) :: x(n) !! On initial entry
                                     !! it must be set by the user to the values of the initial
                                     !! estimate of the solution vector.  Upon successful exit, it
                                     !! contains the values of the variables at the best point
                                     !! found (usually an approximate solution).
      real(wp),intent(in) :: l(n) !! the lower bound on `x`. If
                                  !! the `i`-th variable has no lower bound,
                                  !! `l(i)` need not be defined.
      real(wp),intent(in) :: u(n) !! the upper bound on `x`. If
                                  !! the `i`-th variable has no upper bound,
                                  !! `u(i)` need not be defined.
      integer,intent(in) :: Nbd(n) !! nbd represents the type of bounds imposed on the
                                   !! variables, and must be specified as follows:
                                   !!
                                   !!  * `nbd(i)=0` if `x(i)` is unbounded
                                   !!  * `nbd(i)=1` if `x(i)` has only a lower bound
                                   !!  * `nbd(i)=2` if `x(i)` has both lower and upper bounds
                                   !!  * `nbd(i)=3` if `x(i)` has only an upper bound
      real(wp),intent(inout) :: f !! On first entry `f` is unspecified.
                                  !! On final exit `f` is the value of the function at x.
                                  !! If the [[setulb]] returns
                                  !! with `task(1:2)= 'FG'`, then `f` must be set by the user to
                                  !! contain the value of the function at the point `x`.
      real(wp),intent(inout) :: g(n) !! On first entry `g` is unspecified.
                                     !! On final exit `g` is the value of the gradient at `x`.
                                     !! If the [[setulb]]
                                     !! returns with `task(1:2)= 'FG'`, then `g` must be set by the user to
                                     !! contain the components of the gradient at the point `x`.
      real(wp),intent(in) :: Factr !! A tolerance in the termination test for the algorithm.
                                   !! The iteration will stop when:
                                   !!
                                   !! `(f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch`
                                   !!
                                   !! where `epsmch` is the machine precision.
                                   !! Typical values for `factr` on a computer
                                   !! with 15 digits of accuracy in real(wp) are:
                                   !!
                                   !!  * `factr=1.0e+12` for low accuracy
                                   !!  * `factr=1.0e+7`  for moderate accuracy
                                   !!  * `factr=1.0e+1`  for extremely high accuracy
                                   !!
                                   !! The user can suppress this termination test by setting `factr=0`.
      real(wp),intent(in) :: Pgtol !! The iteration will stop when:
                                   !!
                                   !! `max{|proj g_i | i = 1, ..., n} <= pgtol`
                                   !!
                                   !! where `pg_i` is the `i`th component of the projected gradient.
                                   !! The user can suppress this termination test by setting `pgtol=0`.
      real(wp) :: Wa(2*m*n+5*n+11*m*m+8*m) !! working array of length `(2mmax + 5)nmax + 11mmax^2 + 8mmax`.
                                           !! This array must not be altered by the user.
      integer :: Iwa(3*n) !! an integer working array of length `3nmax`.
                          !! This array must not be altered by the user.
      character(len=60),intent(inout) :: Task !! Indicates the current job when entering and quitting this subroutine:
                                              !!
                                              !!  * On first entry, it must be set to 'START'.
                                              !!  * On a return with task(1:2)='FG', the user must evaluate the
                                              !!    function f and gradient g at the returned value of x.
                                              !!  * On a return with task(1:5)='NEW_X', an iteration of the
                                              !!    algorithm has concluded, and f and g contain f(x) and g(x)
                                              !!    respectively.  The user can decide whether to continue or stop
                                              !!    the iteration.
                                              !!  * When task(1:4)='CONV', the termination test in L-BFGS-B has been
                                              !!    satisfied;
                                              !!  * When task(1:4)='ABNO', the routine has terminated abnormally
                                              !!    without being able to satisfy the termination conditions,
                                              !!    x contains the best approximation found,
                                              !!    f and g contain f(x) and g(x) respectively;
                                              !!  * When task(1:5)='ERROR', the routine has detected an error in the
                                              !!    input parameters;
                                              !!
                                              !!  On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
                                              !!  contains additional information that the user can print.
                                              !!
                                              !! This array should not be altered unless the user wants to
                                              !! stop the run for some reason.  See driver2 or driver3
                                              !! for a detailed explanation on how to stop the run
                                              !! by assigning task(1:4)='STOP' in the driver.
      integer,intent(in) :: Iprint !! Controls the frequency and type of output generated:
                                   !!
                                   !!  * `iprint<0   ` no output is generated
                                   !!  * `iprint=0   ` print only one line at the last iteration
                                   !!  * `0<iprint<99` print also `f` and `|proj g|` every `iprint` iterations
                                   !!  * `iprint=99  ` print details of every iteration except `n`-vectors
                                   !!  * `iprint=100 ` print also the changes of active set and final `x`
                                   !!  * `iprint>100 ` print details of every iteration including `x` and `g`
                                   !!
                                   !! When `iprint > 0`, the file `iterate.dat` will be created to
                                   !! summarize the iteration.
      character(len=60) :: Csave !! working string
      logical :: Lsave(4) !! A logical working array of dimension 4.
                          !! On exit with `task = 'NEW_X'`, the following information is available:
                          !!
                          !!  * If `lsave(1) = .true.` then the initial `x` did not satisfy the bounds
                          !!    and `x` has been replaced by its projection in the feasible set
                          !!  * If `lsave(2) = .true.` then the problem is constrained
                          !!  * If `lsave(3) = .true.` then each variable has upper and lower bounds
      integer :: Isave(44) !! An integer working array of dimension 44.
                           !! On exit with 'task' = NEW_X, the following information is available:
                           !!
                           !!  * isave(22) = the total number of intervals explored in the
                           !!    search of Cauchy points;
                           !!  * isave(26) = the total number of skipped BFGS updates before
                           !!    the current iteration;
                           !!  * isave(30) = the number of current iteration;
                           !!  * isave(31) = the total number of BFGS updates prior the current
                           !!    iteration;
                           !!  * isave(33) = the number of intervals explored in the search of
                           !!    Cauchy point in the current iteration;
                           !!  * isave(34) = the total number of function and gradient
                           !!    evaluations;
                           !!  * isave(36) = the number of function value or gradient
                           !!    evaluations in the current iteration;
                           !!  * if isave(37) = 0  then the subspace argmin is within the box;
                           !!  * if isave(37) = 1  then the subspace argmin is beyond the box;
                           !!  * isave(38) = the number of free variables in the current
                           !!    iteration;
                           !!  * isave(39) = the number of active constraints in the current
                           !!    iteration;
                           !!  * n + 1 - isave(40) = the number of variables leaving the set of
                           !!    active constraints in the current iteration;
                           !!  * isave(41) = the number of variables entering the set of active
                           !!    constraints in the current iteration.
      real(wp) :: Dsave(29) !! A real(wp) working array of dimension 29.
                            !! On exit with 'task' = NEW_X, the following information is available:
                            !!
                            !!  * dsave(1) = current 'theta' in the BFGS matrix;
                            !!  * dsave(2) = f(x) in the previous iteration;
                            !!  * dsave(3) = factr*epsmch;
                            !!  * dsave(4) = 2-norm of the line search direction vector;
                            !!  * dsave(5) = the machine precision epsmch generated by the code;
                            !!  * dsave(7) = the accumulated time spent on searching for
                            !!    Cauchy points;
                            !!  * dsave(8) = the accumulated time spent on
                            !!    subspace minimization;
                            !!  * dsave(9) = the accumulated time spent on line search;
                            !!  * dsave(11) = the slope of the line search function at
                            !!    the current point of line search;
                            !!  * dsave(12) = the maximum relative step length imposed in
                            !!    line search;
                            !!  * dsave(13) = the infinity norm of the projected gradient;
                            !!  * dsave(14) = the relative step length in the line search;
                            !!  * dsave(15) = the slope of the line search function at
                            !!    the starting point of the line search;
                            !!  * dsave(16) = the square of the 2-norm of the line search
                            !!    direction vector.

    integer :: lws , lr , lz , lt , ld , lxp , lwa , lwy , lsy , lss , &
               lwt , lwn , lsnd

      if ( Task=='START' ) then
         Isave(1) = m*n
         Isave(2) = m**2
         Isave(3) = 4*m**2
         Isave(4) = 1                       ! ws      m*n
         Isave(5) = Isave(4) + Isave(1)     ! wy      m*n
         Isave(6) = Isave(5) + Isave(1)     ! wsy     m**2
         Isave(7) = Isave(6) + Isave(2)     ! wss     m**2
         Isave(8) = Isave(7) + Isave(2)     ! wt      m**2
         Isave(9) = Isave(8) + Isave(2)     ! wn      4*m**2
         Isave(10) = Isave(9) + Isave(3)    ! wsnd    4*m**2
         Isave(11) = Isave(10) + Isave(3)   ! wz      n
         Isave(12) = Isave(11) + n          ! wr      n
         Isave(13) = Isave(12) + n          ! wd      n
         Isave(14) = Isave(13) + n          ! wt      n
         Isave(15) = Isave(14) + n          ! wxp     n
         Isave(16) = Isave(15) + n          ! wa      8*m
      endif
      lws = Isave(4)
      lwy = Isave(5)
      lsy = Isave(6)
      lss = Isave(7)
      lwt = Isave(8)
      lwn = Isave(9)
      lsnd = Isave(10)
      lz = Isave(11)
      lr = Isave(12)
      ld = Isave(13)
      lt = Isave(14)
      lxp = Isave(15)
      lwa = Isave(16)

      call mainlb(n,m,x,l,u,Nbd,f,g,Factr,Pgtol,Wa(lws),Wa(lwy),Wa(lsy),&
                & Wa(lss),Wa(lwt),Wa(lwn),Wa(lsnd),Wa(lz),Wa(lr),Wa(ld),&
                & Wa(lt),Wa(lxp),Wa(lwa),Iwa(1),Iwa(n+1),Iwa(2*n+1),    &
                & Task,Iprint,Csave,Lsave,Isave(22),Dsave)

      end subroutine setulb
!*******************************************************************************

!*******************************************************************************
!>
!  This subroutine solves bound constrained optimization problems by
!  using the compact formula of the limited memory BFGS updates.
!
!### References
!  1. R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, "A limited
!     memory algorithm for bound constrained optimization",
!     SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!  2. C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, "L-BFGS-B: FORTRAN
!     Subroutines for Large Scale Bound Constrained Optimization"
!     Tech. Report, NAM-11, EECS Department, Northwestern University,
!     1994.
!  3. R. Byrd, J. Nocedal and R. Schnabel "Representations of
!     Quasi-Newton Matrices and their use in Limited Memory Methods",
!     Mathematical Programming 63 (1994), no. 4, pp. 129-156.
!
!### History
!  * NEOS, November 1994. (Latest revision June 1996.)
!    Optimization Technology Center.
!    Argonne National Laboratory and Northwestern University.
!    Written by Ciyou Zhu in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

      subroutine mainlb(n,m,x,l,u,Nbd,f,g,Factr,Pgtol,Ws,Wy,Sy,Ss,Wt,Wn,&
                        Snd,z,r,d,t,Xp,Wa,Index,Iwhere,Indx2,Task,      &
                        Iprint,Csave,Lsave,Isave,Dsave)
      implicit none

      character*60 Task , Csave
      logical Lsave(4)
      integer n , m , Iprint , Nbd(n) , Index(n) , Iwhere(n) , Indx2(n) &
            & , Isave(23)

      real(wp) f , Factr , Pgtol , x(n) , l(n) , u(n) , g(n) ,  &
                     & z(n) , r(n) , d(n) , t(n) , Xp(n) , Wa(8*m) ,    &
                     & Ws(n,m) , Wy(n,m) , Sy(m,m) , Ss(m,m) , Wt(m,m) ,&
                     & Wn(2*m,2*m) , Snd(2*m,2*m) , Dsave(29)

!
!     n is an integer variable.
!       On entry n is the number of variables.
!       On exit n is unchanged.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric
!          corrections allowed in the limited memory matrix.
!       On exit m is unchanged.
!
!     x is a real(wp) array of dimension n.
!       On entry x is an approximation to the solution.
!       On exit x is the current approximation.
!
!     l is a real(wp) array of dimension n.
!       On entry l is the lower bound of x.
!       On exit l is unchanged.
!
!     u is a real(wp) array of dimension n.
!       On entry u is the upper bound of x.
!       On exit u is unchanged.
!
!     nbd is an integer array of dimension n.
!       On entry nbd represents the type of bounds imposed on the
!         variables, and must be specified as follows:
!         nbd(i)=0 if x(i) is unbounded,
!                1 if x(i) has only a lower bound,
!                2 if x(i) has both lower and upper bounds,
!                3 if x(i) has only an upper bound.
!       On exit nbd is unchanged.
!
!     f is a real(wp) variable.
!       On first entry f is unspecified.
!       On final exit f is the value of the function at x.
!
!     g is a real(wp) array of dimension n.
!       On first entry g is unspecified.
!       On final exit g is the value of the gradient at x.
!
!     factr is a real(wp) variable.
!       On entry factr >= 0 is specified by the user.  The iteration
!         will stop when
!
!         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
!
!         where epsmch is the machine precision, which is automatically
!         generated by the code.
!       On exit factr is unchanged.
!
!     pgtol is a real(wp) variable.
!       On entry pgtol >= 0 is specified by the user.  The iteration
!         will stop when
!
!                 max{|proj g_i | i = 1, ..., n} <= pgtol
!
!         where pg_i is the ith component of the projected gradient.
!       On exit pgtol is unchanged.
!
!     ws, wy, sy, and wt are real(wp) working arrays used to
!       store the following information defining the limited memory
!          BFGS matrix:
!          ws, of dimension n x m, stores S, the matrix of s-vectors;
!          wy, of dimension n x m, stores Y, the matrix of y-vectors;
!          sy, of dimension m x m, stores S'Y;
!          ss, of dimension m x m, stores S'S;
!          yy, of dimension m x m, stores Y'Y;
!          wt, of dimension m x m, stores the Cholesky factorization
!                                  of (theta*S'S+LD^(-1)L'); see eq.
!                                  (2.26) in [3].
!
!     wn is a real(wp) working array of dimension 2m x 2m
!       used to store the LEL^T factorization of the indefinite matrix
!                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!
!       where     E = [-I  0]
!                     [ 0  I]
!
!     snd is a real(wp) working array of dimension 2m x 2m
!       used to store the lower triangular part of
!                 N = [Y' ZZ'Y   L_a'+R_z']
!                     [L_a +R_z  S'AA'S   ]
!
!     z(n),r(n),d(n),t(n), xp(n),wa(8*m) are real(wp) working arrays.
!       z  is used at different times to store the Cauchy point and
!          the Newton point.
!       xp is used to safeguard the projected Newton direction
!
!     sg(m),sgo(m),yg(m),ygo(m) are real(wp) working arrays.
!
!     index is an integer working array of dimension n.
!       In subroutine freev, index is used to store the free and fixed
!          variables at the Generalized Cauchy Point (GCP).
!
!     iwhere is an integer working array of dimension n used to record
!       the status of the vector x for GCP computation.
!       iwhere(i)=0 or -3 if x(i) is free and has bounds,
!                 1       if x(i) is fixed at l(i), and l(i) /= u(i)
!                 2       if x(i) is fixed at u(i), and u(i) /= l(i)
!                 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
!                -1       if x(i) is always free, i.e., no bounds on it.
!
!     indx2 is an integer working array of dimension n.
!       Within subroutine cauchy, indx2 corresponds to the array iorder.
!       In subroutine freev, a list of variables entering and leaving
!       the free set is stored in indx2, and it is passed on to
!       subroutine formk with this information.
!
!     task is a working string of characters of length 60 indicating
!       the current job when entering and leaving this subroutine.
!
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     csave is a working string of characters of length 60.
!
!     lsave is a logical working array of dimension 4.
!
!     isave is an integer working array of dimension 23.
!
!     dsave is a real(wp) working array of dimension 29.


      logical prjctd , cnstnd , boxed , updatd , wrk
      character*3 word
      integer i , k , nintol , itfile , iback , nskip , head , col ,    &
            & iter , itail , iupdat , nseg , nfgv , info , ifun ,       &
            & iword , nfree , nact , ileave , nenter
      real(wp) theta , fold , dr , rr , tol , xstep ,    &
                     & sbgnrm , ddum , dnorm , dtd , epsmch , cpu1 ,    &
                     & cpu2 , cachyt , sbtime , lnscht , time1 , time2 ,&
                     & gd , gdold , stp , stpmx , time

      if ( Task=='START' ) then

         epsmch = epsilon(one)

         call cpu_time(time1)

!        Initialize counters and scalars when task='START'.

!           for the limited memory BFGS matrices:
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         iback = 0
         itail = 0
         iword = 0
         nact = 0
         ileave = 0
         nenter = 0
         fold = zero
         dnorm = zero
         cpu1 = zero
         gd = zero
         stpmx = zero
         sbgnrm = zero
         stp = zero
         gdold = zero
         dtd = zero

!           for operation counts:
         iter = 0
         nfgv = 0
         nseg = 0
         nintol = 0
         nskip = 0
         nfree = n
         ifun = 0
!           for stopping tolerance:
         tol = Factr*epsmch

!           for measuring running time:
         cachyt = 0
         sbtime = 0
         lnscht = 0

!           'word' records the status of subspace solutions.
         word = '---'

!           'info' records the termination information.
         info = 0

         itfile = 8
!                                open a summary file 'iterate.dat'
         if ( Iprint>=1 ) open (8,file='iterate.dat',status='unknown')

!        Check the input arguments for errors.

         call errclb(n,m,Factr,l,u,Nbd,Task,info,k)
         if ( Task(1:5)=='ERROR' ) then
            call prn3lb(n,x,f,Task,Iprint,info,itfile,iter,nfgv,nintol, &
                      & nskip,nact,sbgnrm,zero,nseg,word,iback,stp,     &
                      & xstep,k,cachyt,sbtime,lnscht)
            return
         endif

         call prn1lb(n,m,l,u,x,Iprint,itfile,epsmch)

!        Initialize iwhere & project x onto the feasible set.

         call active(n,l,u,Nbd,x,Iwhere,Iprint,prjctd,cnstnd,boxed)

!        The end of the initialization.

      else
!          restore local variables.

         prjctd = Lsave(1)
         cnstnd = Lsave(2)
         boxed = Lsave(3)
         updatd = Lsave(4)

         nintol = Isave(1)
         itfile = Isave(3)
         iback = Isave(4)
         nskip = Isave(5)
         head = Isave(6)
         col = Isave(7)
         itail = Isave(8)
         iter = Isave(9)
         iupdat = Isave(10)
         nseg = Isave(12)
         nfgv = Isave(13)
         info = Isave(14)
         ifun = Isave(15)
         iword = Isave(16)
         nfree = Isave(17)
         nact = Isave(18)
         ileave = Isave(19)
         nenter = Isave(20)

         theta = Dsave(1)
         fold = Dsave(2)
         tol = Dsave(3)
         dnorm = Dsave(4)
         epsmch = Dsave(5)
         cpu1 = Dsave(6)
         cachyt = Dsave(7)
         sbtime = Dsave(8)
         lnscht = Dsave(9)
         time1 = Dsave(10)
         gd = Dsave(11)
         stpmx = Dsave(12)
         sbgnrm = Dsave(13)
         stp = Dsave(14)
         gdold = Dsave(15)
         dtd = Dsave(16)

!        After returning from the driver go to the point where execution
!        is to resume.

         if ( Task(1:5)=='FG_LN' ) goto 600
         if ( Task(1:5)=='NEW_X' ) goto 700
         if ( Task(1:5)=='FG_ST' ) goto 100
         if ( Task(1:4)=='STOP' ) then
            if ( Task(7:9)=='CPU' ) then
!                                          restore the previous iterate.
               call dcopy(n,t,1,x,1)
               call dcopy(n,r,1,g,1)
               f = fold
            endif
            goto 900
         endif
      endif

!     Compute f0 and g0.

      Task = 'FG_START'
!          return to the driver to calculate f and g; reenter at 111.
      goto 1000
 100  continue
      nfgv = 1

!     Compute the infinity norm of the (-) projected gradient.

      call projgr(n,l,u,Nbd,x,g,sbgnrm)

      if ( Iprint>=1 ) then
         write (6,99001) iter , f , sbgnrm
99001    format (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,&
               & d12.5)
         write (itfile,99002) iter , nfgv , sbgnrm , f
99002    format (2(1x,i4),5x,'-',5x,'-',3x,'-',5x,'-',5x,'-',8x,'-',3x, &
               & 1p,2(1x,d10.3))
      endif
      if ( sbgnrm<=Pgtol ) then
!                                terminate the algorithm.
         Task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         goto 900
      endif

! ----------------- the beginning of the loop --------------------------

 200  continue
      if ( Iprint>=99 ) write (6,99003) iter + 1

99003 format (//,'ITERATION ',i5)
      iword = -1
!
      if ( .not.cnstnd .and. col>0 ) then
!                                            skip the search for GCP.
         call dcopy(n,x,1,z,1)
         wrk = updatd
         nseg = 0
         goto 300
      endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Compute the Generalized Cauchy Point (GCP).
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call cpu_time(cpu1)
      call cauchy(n,x,l,u,Nbd,g,Indx2,Iwhere,t,d,z,m,Wy,Ws,Sy,Wt,theta, &
                & col,head,Wa(1),Wa(2*m+1),Wa(4*m+1),Wa(6*m+1),nseg,    &
                & Iprint,sbgnrm,info,epsmch)
      if ( info/=0 ) then
!         singular triangular system detected; refresh the lbfgs memory.
         if ( Iprint>=1 ) write (6,99008)
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         call cpu_time(cpu2)
         cachyt = cachyt + cpu2 - cpu1
         goto 200
      endif
      call cpu_time(cpu2)
      cachyt = cachyt + cpu2 - cpu1
      nintol = nintol + nseg

!     Count the entering and leaving variables for iter > 0;
!     find the index set of free and active variables at the GCP.

      call freev(n,nfree,Index,nenter,ileave,Indx2,Iwhere,wrk,updatd,   &
               & cnstnd,Iprint,iter)
      nact = n - nfree

 300  continue

!     If there are no free variables or B=theta*I, then
!                                        skip the subspace minimization.

      if ( nfree==0 .or. col==0 ) goto 500

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Subspace minimization.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call cpu_time(cpu1)

!     Form  the LEL^T factorization of the indefinite
!       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!       where     E = [-I  0]
!                     [ 0  I]

      if ( wrk ) call formk(n,nfree,Index,nenter,ileave,Indx2,iupdat,   &
                          & updatd,Wn,Snd,m,Ws,Wy,Sy,theta,col,head,    &
                          & info)
      if ( info/=0 ) then
!          nonpositive definiteness in Cholesky factorization;
!          refresh the lbfgs memory and restart the iteration.
         if ( Iprint>=1 ) write (6,99004)
99004    format (/,                                                     &
        &' Nonpositive definiteness in Cholesky factorization in formk;'&
       & ,/,'   refresh the lbfgs memory and restart the iteration.')
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         call cpu_time(cpu2)
         sbtime = sbtime + cpu2 - cpu1
         goto 200
      endif

!        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
!                                                   from 'cauchy').
      call cmprlb(n,m,x,g,Ws,Wy,Sy,Wt,z,r,Wa,Index,theta,col,head,nfree,&
                & cnstnd,info)
      if ( info/=0 ) goto 400

!-jlm-jn   call the direct method.

      call subsm(n,m,nfree,Index,l,u,Nbd,z,r,Xp,Ws,Wy,theta,x,g,col,    &
               & head,iword,Wa,Wn,Iprint,info)
 400  continue
      if ( info/=0 ) then
!          singular triangular system detected;
!          refresh the lbfgs memory and restart the iteration.
         if ( Iprint>=1 ) write (6,99008)
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         call cpu_time(cpu2)
         sbtime = sbtime + cpu2 - cpu1
         goto 200
      endif

      call cpu_time(cpu2)
      sbtime = sbtime + cpu2 - cpu1
 500  continue

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Line search and optimality tests.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     Generate the search direction d:=z-x.

      do i = 1 , n
         d(i) = z(i) - x(i)
      enddo
      call cpu_time(cpu1)
 600  continue
      call lnsrlb(n,l,u,Nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm,dtd,  &
                & xstep,stpmx,iter,ifun,iback,nfgv,info,Task,boxed,     &
                & cnstnd,Csave,Isave(22),Dsave(17))
      if ( info/=0 .or. iback>=20 ) then
!          restore the previous iterate.
         call dcopy(n,t,1,x,1)
         call dcopy(n,r,1,g,1)
         f = fold
         if ( col==0 ) then
!             abnormal termination.
            if ( info==0 ) then
               info = -9
!                restore the actual number of f and g evaluations etc.
               nfgv = nfgv - 1
               ifun = ifun - 1
               iback = iback - 1
            endif
            Task = 'ABNORMAL_TERMINATION_IN_LNSRCH'
            iter = iter + 1
            goto 900
         else
!             refresh the lbfgs memory and restart the iteration.
            if ( Iprint>=1 ) write (6,99005)
99005       format (/,' Bad direction in the line search;',/,           &
               &'   refresh the lbfgs memory and restart the iteration.'&
              & )
            if ( info==0 ) nfgv = nfgv - 1
            info = 0
            col = 0
            head = 1
            theta = one
            iupdat = 0
            updatd = .false.
            Task = 'RESTART_FROM_LNSRCH'
            call cpu_time(cpu2)
            lnscht = lnscht + cpu2 - cpu1
            goto 200
         endif
      elseif ( Task(1:5)=='FG_LN' ) then
!          return to the driver for calculating f and g; reenter at 666.
         goto 1000
      else
!          calculate and print out the quantities related to the new X.
         call cpu_time(cpu2)
         lnscht = lnscht + cpu2 - cpu1
         iter = iter + 1

!        Compute the infinity norm of the projected (-)gradient.

         call projgr(n,l,u,Nbd,x,g,sbgnrm)

!        Print iteration information.

         call prn2lb(n,x,f,g,Iprint,itfile,iter,nfgv,nact,sbgnrm,nseg,  &
                   & word,iword,iback,stp,xstep)
         goto 1000
      endif
 700  continue

!     Test for termination.

      if ( sbgnrm<=Pgtol ) then
!                                terminate the algorithm.
         Task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         goto 900
      endif

      ddum = max(abs(fold),abs(f),one)
      if ( (fold-f)<=tol*ddum ) then
!                                        terminate the algorithm.
         Task = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
         if ( iback>=10 ) info = -5
!           i.e., to issue a warning if iback>10 in the line search.
         goto 900
      endif

!     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.

      do i = 1 , n
         r(i) = g(i) - r(i)
      enddo
      rr = ddot(n,r,1,r,1)
      if ( stp==one ) then
         dr = gd - gdold
         ddum = -gdold
      else
         dr = (gd-gdold)*stp
         call dscal(n,stp,d,1)
         ddum = -gdold*stp
      endif

      if ( dr<=epsmch*ddum ) then
!                            skip the L-BFGS update.
         nskip = nskip + 1
         updatd = .false.
         if ( Iprint>=1 ) write (6,99006) dr , ddum
99006    format ('  ys=',1p,e10.3,'  -gs=',1p,e10.3,                    &
                &' BFGS update SKIPPED')
         goto 800
      endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Update the L-BFGS matrix.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      updatd = .true.
      iupdat = iupdat + 1

!     Update matrices WS and WY and form the middle matrix in B.

      call matupd(n,m,Ws,Wy,Sy,Ss,d,r,itail,iupdat,col,head,theta,rr,dr,&
                & stp,dtd)

!     Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
!        Store T in the upper triangular of the array wt;
!        Cholesky factorize T to J*J' with
!           J' stored in the upper triangular of wt.

      call formt(m,Wt,Sy,Ss,col,theta,info)

      if ( info/=0 ) then
!          nonpositive definiteness in Cholesky factorization;
!          refresh the lbfgs memory and restart the iteration.
         if ( Iprint>=1 ) write (6,99007)
99007    format (/,                                                     &
        &' Nonpositive definiteness in Cholesky factorization in formt;'&
       & ,/,'   refresh the lbfgs memory and restart the iteration.')
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         goto 200
      endif

!     Now the inverse of the middle matrix in B is

!       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
!       [ -L*D^(-1/2)   J ] [  0        J'          ]

 800  continue

! -------------------- the end of the loop -----------------------------

      goto 200
 900  continue
      call cpu_time(time2)
      time = time2 - time1
      call prn3lb(n,x,f,Task,Iprint,info,itfile,iter,nfgv,nintol,nskip, &
                & nact,sbgnrm,time,nseg,word,iback,stp,xstep,k,cachyt,  &
                & sbtime,lnscht)
 1000 continue

!     Save local variables.

      Lsave(1) = prjctd
      Lsave(2) = cnstnd
      Lsave(3) = boxed
      Lsave(4) = updatd

      Isave(1) = nintol
      Isave(3) = itfile
      Isave(4) = iback
      Isave(5) = nskip
      Isave(6) = head
      Isave(7) = col
      Isave(8) = itail
      Isave(9) = iter
      Isave(10) = iupdat
      Isave(12) = nseg
      Isave(13) = nfgv
      Isave(14) = info
      Isave(15) = ifun
      Isave(16) = iword
      Isave(17) = nfree
      Isave(18) = nact
      Isave(19) = ileave
      Isave(20) = nenter

      Dsave(1) = theta
      Dsave(2) = fold
      Dsave(3) = tol
      Dsave(4) = dnorm
      Dsave(5) = epsmch
      Dsave(6) = cpu1
      Dsave(7) = cachyt
      Dsave(8) = sbtime
      Dsave(9) = lnscht
      Dsave(10) = time1
      Dsave(11) = gd
      Dsave(12) = stpmx
      Dsave(13) = sbgnrm
      Dsave(14) = stp
      Dsave(15) = gdold
      Dsave(16) = dtd

      continue
99008 format (/,' Singular triangular system detected;',/,              &
             &'   refresh the lbfgs memory and restart the iteration.')

      end subroutine mainlb
!*******************************************************************************

!*******************************************************************************
!>
!  This subroutine initializes iwhere and projects the initial x to
!  the feasible set if necessary.

      subroutine active(n,l,u,Nbd,x,Iwhere,Iprint,Prjctd,Cnstnd,Boxed)
      implicit none

      logical Prjctd , Cnstnd , Boxed
      integer n , Iprint , Nbd(n) , Iwhere(n)
      real(wp) x(n) , l(n) , u(n)

!
!     iwhere is an integer array of dimension n.
!       On entry iwhere is unspecified.
!       On exit iwhere(i)=-1  if x(i) has no bounds
!                         3   if l(i)=u(i)
!                         0   otherwise.
!       In cauchy, iwhere is given finer gradations.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer nbdd , i

!     Initialize nbdd, prjctd, cnstnd and boxed.

      nbdd = 0
      Prjctd = .false.
      Cnstnd = .false.
      Boxed = .true.

!     Project the initial x to the feasible set if necessary.

      do i = 1 , n
         if ( Nbd(i)>0 ) then
            if ( Nbd(i)<=2 .and. x(i)<=l(i) ) then
               if ( x(i)<l(i) ) then
                  Prjctd = .true.
                  x(i) = l(i)
               endif
               nbdd = nbdd + 1
            elseif ( Nbd(i)>=2 .and. x(i)>=u(i) ) then
               if ( x(i)>u(i) ) then
                  Prjctd = .true.
                  x(i) = u(i)
               endif
               nbdd = nbdd + 1
            endif
         endif
      enddo

!     Initialize iwhere and assign values to cnstnd and boxed.

      do i = 1 , n
         if ( Nbd(i)/=2 ) Boxed = .false.
         if ( Nbd(i)==0 ) then
!                                this variable is always free
            Iwhere(i) = -1

!           otherwise set x(i)=mid(x(i), u(i), l(i)).
         else
            Cnstnd = .true.
            if ( Nbd(i)==2 .and. u(i)-l(i)<=zero ) then
!                   this variable is always fixed
               Iwhere(i) = 3
            else
               Iwhere(i) = 0
            endif
         endif
      enddo

      if ( Iprint>=0 ) then
         if ( Prjctd ) write (6,*)                                      &
           &'The initial X is infeasible.  Restart with its projection.'
         if ( .not.Cnstnd ) write (6,*) 'This problem is unconstrained.'
      endif

      if ( Iprint>0 ) write (6,99001) nbdd

99001 format (/,'At X0 ',i9,' variables are exactly at the bounds')

      end subroutine active

      subroutine bmv(m,Sy,Wt,Col,v,p,Info)
      implicit none

      integer m , Col , Info
      real(wp) Sy(m,m) , Wt(m,m) , v(2*Col) , p(2*Col)

!     ************
!
!     Subroutine bmv
!
!     This subroutine computes the product of the 2m x 2m middle matrix
!       in the compact L-BFGS formula of B and a 2m vector v;
!       it returns the product in p.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     sy is a real(wp) array of dimension m x m.
!       On entry sy specifies the matrix S'Y.
!       On exit sy is unchanged.
!
!     wt is a real(wp) array of dimension m x m.
!       On entry wt specifies the upper triangular matrix J' which is
!         the Cholesky factor of (thetaS'S+LD^(-1)L').
!       On exit wt is unchanged.
!
!     col is an integer variable.
!       On entry col specifies the number of s-vectors (or y-vectors)
!         stored in the compact L-BFGS formula.
!       On exit col is unchanged.
!
!     v is a real(wp) array of dimension 2col.
!       On entry v specifies vector v.
!       On exit v is unchanged.
!
!     p is a real(wp) array of dimension 2col.
!       On entry p is unspecified.
!       On exit p is the product Mv.
!
!     info is an integer variable.
!       On entry info is unspecified.
!       On exit info = 0       for normal return,
!                    = nonzero for abnormal return when the system
!                                to be solved by dtrsl is singular.
!
!     Subprograms called:
!
!       Linpack ... dtrsl.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer i , k , i2
      real(wp) sum

      if ( Col==0 ) return

!     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
!                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].

!       solve Jp2=v2+LD^(-1)v1.
      p(Col+1) = v(Col+1)
      do i = 2 , Col
         i2 = Col + i
         sum = 0.0d0
         do k = 1 , i - 1
            sum = sum + Sy(i,k)*v(k)/Sy(k,k)
         enddo
         p(i2) = v(i2) + sum
      enddo
!     Solve the triangular system
      call dtrsl(Wt,m,Col,p(Col+1),11,Info)
      if ( Info/=0 ) return

!       solve D^(1/2)p1=v1.
      do i = 1 , Col
         p(i) = v(i)/sqrt(Sy(i,i))
      enddo

!     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
!                    [  0         J'           ] [ p2 ]   [ p2 ].

!       solve J^Tp2=p2.
      call dtrsl(Wt,m,Col,p(Col+1),01,Info)
      if ( Info/=0 ) return

!       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
!                 =-D^(-1/2)p1+D^(-1)L'p2.
      do i = 1 , Col
         p(i) = -p(i)/sqrt(Sy(i,i))
      enddo
      do i = 1 , Col
         sum = 0.d0
         do k = i + 1 , Col
            sum = sum + Sy(k,i)*p(Col+k)/Sy(i,i)
         enddo
         p(i) = p(i) + sum
      enddo

      end subroutine bmv

      subroutine cauchy(n,x,l,u,Nbd,g,Iorder,Iwhere,t,d,Xcp,m,Wy,Ws,Sy, &
                      & Wt,Theta,Col,Head,p,c,Wbp,v,Nseg,Iprint,Sbgnrm, &
                      & Info,Epsmch)
      implicit none

      integer n , m , Head , Col , Nseg , Iprint , Info , Nbd(n) ,      &
            & Iorder(n) , Iwhere(n)
      real(wp) Theta , Epsmch , x(n) , l(n) , u(n) , g(n) ,     &
                     & t(n) , d(n) , Xcp(n) , Wy(n,Col) , Ws(n,Col) ,   &
                     & Sy(m,m) , Wt(m,m) , p(2*m) , c(2*m) , Wbp(2*m) , &
                     & v(2*m)

!     ************
!
!     Subroutine cauchy
!
!     For given x, l, u, g (with sbgnrm > 0), and a limited memory
!       BFGS matrix B defined in terms of matrices WY, WS, WT, and
!       scalars head, col, and theta, this subroutine computes the
!       generalized Cauchy point (GCP), defined as the first local
!       minimizer of the quadratic
!
!                  Q(x + s) = g's + 1/2 s'Bs
!
!       along the projected gradient direction P(x-tg,l,u).
!       The routine returns the GCP in xcp.
!
!     n is an integer variable.
!       On entry n is the dimension of the problem.
!       On exit n is unchanged.
!
!     x is a real(wp) array of dimension n.
!       On entry x is the starting point for the GCP computation.
!       On exit x is unchanged.
!
!     l is a real(wp) array of dimension n.
!       On entry l is the lower bound of x.
!       On exit l is unchanged.
!
!     u is a real(wp) array of dimension n.
!       On entry u is the upper bound of x.
!       On exit u is unchanged.
!
!     nbd is an integer array of dimension n.
!       On entry nbd represents the type of bounds imposed on the
!         variables, and must be specified as follows:
!         nbd(i)=0 if x(i) is unbounded,
!                1 if x(i) has only a lower bound,
!                2 if x(i) has both lower and upper bounds, and
!                3 if x(i) has only an upper bound.
!       On exit nbd is unchanged.
!
!     g is a real(wp) array of dimension n.
!       On entry g is the gradient of f(x).  g must be a nonzero vector.
!       On exit g is unchanged.
!
!     iorder is an integer working array of dimension n.
!       iorder will be used to store the breakpoints in the piecewise
!       linear path and free variables encountered. On exit,
!         iorder(1),...,iorder(nleft) are indices of breakpoints
!                                which have not been encountered;
!         iorder(nleft+1),...,iorder(nbreak) are indices of
!                                     encountered breakpoints; and
!         iorder(nfree),...,iorder(n) are indices of variables which
!                 have no bound constraits along the search direction.
!
!     iwhere is an integer array of dimension n.
!       On entry iwhere indicates only the permanently fixed (iwhere=3)
!       or free (iwhere= -1) components of x.
!       On exit iwhere records the status of the current x variables.
!       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
!                 0   if x(i) is free and has bounds, and is moved
!                 1   if x(i) is fixed at l(i), and l(i) /= u(i)
!                 2   if x(i) is fixed at u(i), and u(i) /= l(i)
!                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
!                 -1  if x(i) is always free, i.e., it has no bounds.
!
!     t is a real(wp) working array of dimension n.
!       t will be used to store the break points.
!
!     d is a real(wp) array of dimension n used to store
!       the Cauchy direction P(x-tg)-x.
!
!     xcp is a real(wp) array of dimension n used to return the
!       GCP on exit.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     ws, wy, sy, and wt are real(wp) arrays.
!       On entry they store information that defines the
!                             limited memory BFGS matrix:
!         ws(n,m) stores S, a set of s-vectors;
!         wy(n,m) stores Y, a set of y-vectors;
!         sy(m,m) stores S'Y;
!         wt(m,m) stores the
!                 Cholesky factorization of (theta*S'S+LD^(-1)L').
!       On exit these arrays are unchanged.
!
!     theta is a real(wp) variable.
!       On entry theta is the scaling factor specifying B_0 = theta I.
!       On exit theta is unchanged.
!
!     col is an integer variable.
!       On entry col is the actual number of variable metric
!         corrections stored so far.
!       On exit col is unchanged.
!
!     head is an integer variable.
!       On entry head is the location of the first s-vector (or y-vector)
!         in S (or Y).
!       On exit col is unchanged.
!
!     p is a real(wp) working array of dimension 2m.
!       p will be used to store the vector p = W^(T)d.
!
!     c is a real(wp) working array of dimension 2m.
!       c will be used to store the vector c = W^(T)(xcp-x).
!
!     wbp is a real(wp) working array of dimension 2m.
!       wbp will be used to store the row of W corresponding
!         to a breakpoint.
!
!     v is a real(wp) working array of dimension 2m.
!
!     nseg is an integer variable.
!       On exit nseg records the number of quadratic segments explored
!         in searching for the GCP.
!
!     sg and yg are real(wp) arrays of dimension m.
!       On entry sg  and yg store S'g and Y'g correspondingly.
!       On exit they are unchanged.
!
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     sbgnrm is a real(wp) variable.
!       On entry sbgnrm is the norm of the projected gradient at x.
!       On exit sbgnrm is unchanged.
!
!     info is an integer variable.
!       On entry info is 0.
!       On exit info = 0       for normal return,
!                    = nonzero for abnormal return when the the system
!                              used in routine bmv is singular.
!
!     Subprograms called:
!
!       L-BFGS-B Library ... hpsolb, bmv.
!
!       Linpack ... dscal dcopy, daxpy.
!
!
!     References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, "A limited
!       memory algorithm for bound constrained optimization",
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, "L-BFGS-B: FORTRAN
!       Subroutines for Large Scale Bound Constrained Optimization"
!       Tech. Report, NAM-11, EECS Department, Northwestern University,
!       1994.
!
!       (Postscript files of these papers are available via anonymous
!        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      logical xlower , xupper , bnded
      integer i , j , col2 , nfree , nbreak , pointr , ibp , nleft ,    &
            & ibkmin , iter
      real(wp) f1 , f2 , dt , dtm , tsum , dibp , zibp , dibp2 ,&
                     & bkmin , tu , tl , wmc , wmp , wmw , tj ,  &
                     & tj0 , neggi , Sbgnrm , f2_org

!     Check the status of the variables, reset iwhere(i) if necessary;
!       compute the Cauchy direction d and the breakpoints t; initialize
!       the derivative f1 and the vector p = W'd (for theta = 1).

      if ( Sbgnrm<=zero ) then
         if ( Iprint>=0 ) write (6,*) 'Subgnorm = 0.  GCP = X.'
         call dcopy(n,x,1,Xcp,1)
         return
      endif
      bnded = .true.
      nfree = n + 1
      nbreak = 0
      ibkmin = 0
      bkmin = zero
      col2 = 2*Col
      f1 = zero
      if ( Iprint>=99 ) write (6,99001)
99001 format (/,'---------------- CAUCHY entered-------------------')

!     We set p to zero and build it up as we determine d.

      do i = 1 , col2
         p(i) = zero
      enddo

!     In the following loop we determine for each variable its bound
!        status and its breakpoint, and update p accordingly.
!        Smallest breakpoint is identified.

      do i = 1 , n
         neggi = -g(i)
         if ( Iwhere(i)/=3 .and. Iwhere(i)/=-1 ) then
!             if x(i) is not a constant and has bounds,
!             compute the difference between x(i) and its bounds.
            if ( Nbd(i)<=2 ) tl = x(i) - l(i)
            if ( Nbd(i)>=2 ) tu = u(i) - x(i)

!           If a variable is close enough to a bound
!             we treat it as at bound.
            xlower = Nbd(i)<=2 .and. tl<=zero
            xupper = Nbd(i)>=2 .and. tu<=zero

!              reset iwhere(i).
            Iwhere(i) = 0
            if ( xlower ) then
               if ( neggi<=zero ) Iwhere(i) = 1
            elseif ( xupper ) then
               if ( neggi>=zero ) Iwhere(i) = 2
            else
               if ( abs(neggi)<=zero ) Iwhere(i) = -3
            endif
         endif
         pointr = Head
         if ( Iwhere(i)/=0 .and. Iwhere(i)/=-1 ) then
            d(i) = zero
         else
            d(i) = neggi
            f1 = f1 - neggi*neggi
!             calculate p := p - W'e_i* (g_i).
            do j = 1 , Col
               p(j) = p(j) + Wy(i,pointr)*neggi
               p(Col+j) = p(Col+j) + Ws(i,pointr)*neggi
               pointr = mod(pointr,m) + 1
            enddo
            if ( Nbd(i)<=2 .and. Nbd(i)/=0 .and. neggi<zero )    &
               & then
!                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               Iorder(nbreak) = i
               t(nbreak) = tl/(-neggi)
               if ( nbreak==1 .or. t(nbreak)<bkmin ) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            elseif ( Nbd(i)>=2 .and. neggi>zero ) then
!                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               Iorder(nbreak) = i
               t(nbreak) = tu/neggi
               if ( nbreak==1 .or. t(nbreak)<bkmin ) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else
!                x(i) + d(i) is not bounded.
               nfree = nfree - 1
               Iorder(nfree) = i
               if ( abs(neggi)>zero ) bnded = .false.
            endif
         endif
      enddo

!     The indices of the nonzero components of d are now stored
!       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
!       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.

!                   complete the initialization of p for theta not= one.
      if ( Theta/=one ) call dscal(Col,Theta,p(Col+1),1)

!     Initialize GCP xcp = x.

      call dcopy(n,x,1,Xcp,1)

      if ( nbreak==0 .and. nfree==n+1 ) then
!                  is a zero vector, return with the initial xcp as GCP.
         if ( Iprint>100 ) write (6,99006) (Xcp(i),i=1,n)
         return
      endif

!     Initialize c = W'(xcp - x) = 0.

      do j = 1 , col2
         c(j) = zero
      enddo

!     Initialize derivative f2.

      f2 = -Theta*f1
      f2_org = f2
      if ( Col>0 ) then
         call bmv(m,Sy,Wt,Col,p,v,Info)
         if ( Info/=0 ) return
         f2 = f2 - ddot(col2,v,1,p,1)
      endif
      dtm = -f1/f2
      tsum = zero
      Nseg = 1
      if ( Iprint>=99 ) write (6,*) 'There are ' , nbreak ,           &
                                     &'  breakpoints '

!     If there are no breakpoints, locate the GCP and return.

      if ( nbreak==0 ) goto 200

      nleft = nbreak
      iter = 1


      tj = zero

!------------------- the beginning of the loop -------------------------

 100  continue

!     Find the next smallest breakpoint;
!       compute dt = t(nleft) - t(nleft + 1).

      tj0 = tj
      if ( iter==1 ) then
!         Since we already have the smallest breakpoint we need not do
!         heapsort yet. Often only one breakpoint is used and the
!         cost of heapsort is avoided.
         tj = bkmin
         ibp = Iorder(ibkmin)
      else
         if ( iter==2 ) then
!             Replace the already used smallest breakpoint with the
!             breakpoint numbered nbreak > nlast, before heapsort call.
            if ( ibkmin/=nbreak ) then
               t(ibkmin) = t(nbreak)
               Iorder(ibkmin) = Iorder(nbreak)
            endif
!        Update heap structure of breakpoints
!           (if iter=2, initialize heap).
         endif
         call hpsolb(nleft,t,Iorder,iter-2)
         tj = t(nleft)
         ibp = Iorder(nleft)
      endif

      dt = tj - tj0

      if ( dt/=zero .and. Iprint>=100 ) then
         write (6,99002) Nseg , f1 , f2
99002    format (/,'Piece    ',i3,' --f1, f2 at start point ',1p,       &
               & 2(1x,d11.4))
         write (6,99003) dt
99003    format ('Distance to the next break point =  ',1p,d11.4)
         write (6,99007) dtm
      endif

!     If a minimizer is within this interval, locate the GCP and return.

      if ( dtm<dt ) goto 200

!     Otherwise fix one variable and
!       reset the corresponding component of d to zero.

      tsum = tsum + dt
      nleft = nleft - 1
      iter = iter + 1
      dibp = d(ibp)
      d(ibp) = zero
      if ( dibp>zero ) then
         zibp = u(ibp) - x(ibp)
         Xcp(ibp) = u(ibp)
         Iwhere(ibp) = 2
      else
         zibp = l(ibp) - x(ibp)
         Xcp(ibp) = l(ibp)
         Iwhere(ibp) = 1
      endif
      if ( Iprint>=100 ) write (6,*) 'Variable  ' , ibp ,             &
                                      &'  is fixed.'
      if ( nleft==0 .and. nbreak==n ) then
!                                             all n variables are fixed,
!                                                return with xcp as GCP.
         dtm = dt
         goto 300
      endif

!     Update the derivative information.

      Nseg = Nseg + 1
      dibp2 = dibp**2

!     Update f1 and f2.

!        temporarily set f1 and f2 for col=0.
      f1 = f1 + dt*f2 + dibp2 - Theta*dibp*zibp
      f2 = f2 - Theta*dibp2

      if ( Col>0 ) then
!                          update c = c + dt*p.
         call daxpy(col2,dt,p,1,c,1)

!           choose wbp,
!           the row of W corresponding to the breakpoint encountered.
         pointr = Head
         do j = 1 , Col
            Wbp(j) = Wy(ibp,pointr)
            Wbp(Col+j) = Theta*Ws(ibp,pointr)
            pointr = mod(pointr,m) + 1
         enddo

!           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
         call bmv(m,Sy,Wt,Col,Wbp,v,Info)
         if ( Info/=0 ) return
         wmc = ddot(col2,c,1,v,1)
         wmp = ddot(col2,p,1,v,1)
         wmw = ddot(col2,Wbp,1,v,1)

!           update p = p - dibp*wbp.
         call daxpy(col2,-dibp,Wbp,1,p,1)

!           complete updating f1 and f2 while col > 0.
         f1 = f1 + dibp*wmc
         f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
      endif

      f2 = max(Epsmch*f2_org,f2)
      if ( nleft>0 ) then
         dtm = -f1/f2
         goto 100
!                 to repeat the loop for unsearched intervals.
      elseif ( bnded ) then
         f1 = zero
         f2 = zero
         dtm = zero
      else
         dtm = -f1/f2
      endif

!------------------- the end of the loop -------------------------------

 200  continue
      if ( Iprint>=99 ) then
         write (6,*)
         write (6,*) 'GCP found in this segment'
         write (6,99004) Nseg , f1 , f2
99004    format ('Piece    ',i3,' --f1, f2 at start point ',1p,         &
               & 2(1x,d11.4))
         write (6,99007) dtm
      endif
      if ( dtm<=zero ) dtm = zero
      tsum = tsum + dtm

!     Move free variables (i.e., the ones w/o breakpoints) and
!       the variables whose breakpoints haven't been reached.

      call daxpy(n,tsum,d,1,Xcp,1)

 300  continue

!     Update c = c + dtm*p = W'(x^c - x)
!       which will be used in computing r = Z'(B(x^c - x) + g).

      if ( Col>0 ) call daxpy(col2,dtm,p,1,c,1)
      if ( Iprint>100 ) write (6,99006) (Xcp(i),i=1,n)
      if ( Iprint>=99 ) write (6,99005)
99005 format (/,'---------------- exit CAUCHY----------------------',/)

      continue

99006 format ('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))
99007 format ('Distance to the stationary point =  ',1p,d11.4)

      end subroutine cauchy

      subroutine cmprlb(n,m,x,g,Ws,Wy,Sy,Wt,z,r,Wa,Index,Theta,Col,Head,&
                      & Nfree,Cnstnd,Info)
      implicit none

      logical Cnstnd
      integer n , m , Col , Head , Nfree , Info , Index(n)
      real(wp) Theta , x(n) , g(n) , z(n) , r(n) , Wa(4*m) ,    &
                     & Ws(n,m) , Wy(n,m) , Sy(m,m) , Wt(m,m)

!     ************
!
!     Subroutine cmprlb
!
!       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
!         wa(2m+1)=W'(xcp-x) from subroutine cauchy.
!
!     Subprograms called:
!
!       L-BFGS-B Library ... bmv.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer i , j , k , pointr
      real(wp) a1 , a2

      if ( .not.Cnstnd .and. Col>0 ) then
         do i = 1 , n
            r(i) = -g(i)
         enddo
      else
         do i = 1 , Nfree
            k = Index(i)
            r(i) = -Theta*(z(k)-x(k)) - g(k)
         enddo
         call bmv(m,Sy,Wt,Col,Wa(2*m+1),Wa(1),Info)
         if ( Info/=0 ) then
            Info = -8
            return
         endif
         pointr = Head
         do j = 1 , Col
            a1 = Wa(j)
            a2 = Theta*Wa(Col+j)
            do i = 1 , Nfree
               k = Index(i)
               r(i) = r(i) + Wy(k,pointr)*a1 + Ws(k,pointr)*a2
            enddo
            pointr = mod(pointr,m) + 1
         enddo
      endif

      end subroutine cmprlb

      subroutine errclb(n,m,Factr,l,u,Nbd,Task,Info,k)
      implicit none

      character*60 Task
      integer n , m , Info , k , Nbd(n)
      real(wp) Factr , l(n) , u(n)

!     ************
!
!     Subroutine errclb
!
!     This subroutine checks the validity of the input data.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer i

!     Check the input arguments for errors.

      if ( n<=0 ) Task = 'ERROR: N <= 0'
      if ( m<=0 ) Task = 'ERROR: M <= 0'
      if ( Factr<zero ) Task = 'ERROR: FACTR < 0'

!     Check the validity of the arrays nbd(i), u(i), and l(i).

      do i = 1 , n
         if ( Nbd(i)<0 .or. Nbd(i)>3 ) then
!                                                   return
            Task = 'ERROR: INVALID NBD'
            Info = -6
            k = i
         endif
         if ( Nbd(i)==2 ) then
            if ( l(i)>u(i) ) then
!                                    return
               Task = 'ERROR: NO FEASIBLE SOLUTION'
               Info = -7
               k = i
            endif
         endif
      enddo

      end subroutine errclb

      subroutine formk(n,Nsub,Ind,Nenter,Ileave,Indx2,Iupdat,Updatd,Wn, &
                     & Wn1,m,Ws,Wy,Sy,Theta,Col,Head,Info)
      implicit none

      integer n , Nsub , m , Col , Head , Nenter , Ileave , Iupdat ,    &
            & Info , Ind(n) , Indx2(n)
      real(wp) Theta , Wn(2*m,2*m) , Wn1(2*m,2*m) , Ws(n,m) ,   &
                     & Wy(n,m) , Sy(m,m)
      logical Updatd

!     ************
!
!     Subroutine formk
!
!     This subroutine forms  the LEL^T factorization of the indefinite
!
!       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!                                                    where E = [-I  0]
!                                                              [ 0  I]
!     The matrix K can be shown to be equal to the matrix M^[-1]N
!       occurring in section 5.1 of [1], as well as to the matrix
!       Mbar^[-1] Nbar in section 5.3.
!
!     n is an integer variable.
!       On entry n is the dimension of the problem.
!       On exit n is unchanged.
!
!     nsub is an integer variable
!       On entry nsub is the number of subspace variables in free set.
!       On exit nsub is not changed.
!
!     ind is an integer array of dimension nsub.
!       On entry ind specifies the indices of subspace variables.
!       On exit ind is unchanged.
!
!     nenter is an integer variable.
!       On entry nenter is the number of variables entering the
!         free set.
!       On exit nenter is unchanged.
!
!     ileave is an integer variable.
!       On entry indx2(ileave),...,indx2(n) are the variables leaving
!         the free set.
!       On exit ileave is unchanged.
!
!     indx2 is an integer array of dimension n.
!       On entry indx2(1),...,indx2(nenter) are the variables entering
!         the free set, while indx2(ileave),...,indx2(n) are the
!         variables leaving the free set.
!       On exit indx2 is unchanged.
!
!     iupdat is an integer variable.
!       On entry iupdat is the total number of BFGS updates made so far.
!       On exit iupdat is unchanged.
!
!     updatd is a logical variable.
!       On entry 'updatd' is true if the L-BFGS matrix is updatd.
!       On exit 'updatd' is unchanged.
!
!     wn is a real(wp) array of dimension 2m x 2m.
!       On entry wn is unspecified.
!       On exit the upper triangle of wn stores the LEL^T factorization
!         of the 2*col x 2*col indefinite matrix
!                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!
!     wn1 is a real(wp) array of dimension 2m x 2m.
!       On entry wn1 stores the lower triangular part of
!                     [Y' ZZ'Y   L_a'+R_z']
!                     [L_a+R_z   S'AA'S   ]
!         in the previous iteration.
!       On exit wn1 stores the corresponding updated matrices.
!       The purpose of wn1 is just to store these inner products
!       so they can be easily updated and inserted into wn.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     ws, wy, sy, and wtyy are real(wp) arrays;
!     theta is a real(wp) variable;
!     col is an integer variable;
!     head is an integer variable.
!       On entry they store the information defining the
!                                          limited memory BFGS matrix:
!         ws(n,m) stores S, a set of s-vectors;
!         wy(n,m) stores Y, a set of y-vectors;
!         sy(m,m) stores S'Y;
!         wtyy(m,m) stores the Cholesky factorization
!                                   of (theta*S'S+LD^(-1)L')
!         theta is the scaling factor specifying B_0 = theta I;
!         col is the number of variable metric corrections stored;
!         head is the location of the 1st s- (or y-) vector in S (or Y).
!       On exit they are unchanged.
!
!     info is an integer variable.
!       On entry info is unspecified.
!       On exit info =  0 for normal return;
!                    = -1 when the 1st Cholesky factorization failed;
!                    = -2 when the 2st Cholesky factorization failed.
!
!     Subprograms called:
!
!       Linpack ... dcopy, dpofa, dtrsl.
!
!
!     References:
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, "A limited
!       memory algorithm for bound constrained optimization",
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, "L-BFGS-B: a
!       limited memory FORTRAN code for solving bound constrained
!       optimization problems", Tech. Report, NAM-11, EECS Department,
!       Northwestern University, 1994.
!
!       (Postscript files of these papers are available via anonymous
!        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer m2 , ipntr , jpntr , iy , is , jy , js , is1 , js1 , k1 , &
            & i , k , col2 , pbegin , pend , dbegin , dend , upcl
      real(wp) temp1 , temp2 , temp3 , temp4

!     Form the lower triangular part of
!               WN1 = [Y' ZZ'Y   L_a'+R_z']
!                     [L_a+R_z   S'AA'S   ]
!        where L_a is the strictly lower triangular part of S'AA'Y
!              R_z is the upper triangular part of S'ZZ'Y.

      if ( Updatd ) then
         if ( Iupdat>m ) then
!                                 shift old part of WN1.
            do jy = 1 , m - 1
               js = m + jy
               call dcopy(m-jy,Wn1(jy+1,jy+1),1,Wn1(jy,jy),1)
               call dcopy(m-jy,Wn1(js+1,js+1),1,Wn1(js,js),1)
               call dcopy(m-1,Wn1(m+2,jy+1),1,Wn1(m+1,jy),1)
            enddo
         endif

!          put new rows in blocks (1,1), (2,1) and (2,2).
         pbegin = 1
         pend = Nsub
         dbegin = Nsub + 1
         dend = n
         iy = Col
         is = m + Col
         ipntr = Head + Col - 1
         if ( ipntr>m ) ipntr = ipntr - m
         jpntr = Head
         do jy = 1 , Col
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
!             compute element jy of row 'col' of Y'ZZ'Y
            do k = pbegin , pend
               k1 = Ind(k)
               temp1 = temp1 + Wy(k1,ipntr)*Wy(k1,jpntr)
            enddo
!             compute elements jy of row 'col' of L_a and S'AA'S
            do k = dbegin , dend
               k1 = Ind(k)
               temp2 = temp2 + Ws(k1,ipntr)*Ws(k1,jpntr)
               temp3 = temp3 + Ws(k1,ipntr)*Wy(k1,jpntr)
            enddo
            Wn1(iy,jy) = temp1
            Wn1(is,js) = temp2
            Wn1(is,jy) = temp3
            jpntr = mod(jpntr,m) + 1
         enddo

!          put new column in block (2,1).
         jy = Col
         jpntr = Head + Col - 1
         if ( jpntr>m ) jpntr = jpntr - m
         ipntr = Head
         do i = 1 , Col
            is = m + i
            temp3 = zero
!             compute element i of column 'col' of R_z
            do k = pbegin , pend
               k1 = Ind(k)
               temp3 = temp3 + Ws(k1,ipntr)*Wy(k1,jpntr)
            enddo
            ipntr = mod(ipntr,m) + 1
            Wn1(is,jy) = temp3
         enddo
         upcl = Col - 1
      else
         upcl = Col
      endif

!       modify the old parts in blocks (1,1) and (2,2) due to changes
!       in the set of free variables.
      ipntr = Head
      do iy = 1 , upcl
         is = m + iy
         jpntr = Head
         do jy = 1 , iy
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
            temp4 = zero
            do k = 1 , Nenter
               k1 = Indx2(k)
               temp1 = temp1 + Wy(k1,ipntr)*Wy(k1,jpntr)
               temp2 = temp2 + Ws(k1,ipntr)*Ws(k1,jpntr)
            enddo
            do k = Ileave , n
               k1 = Indx2(k)
               temp3 = temp3 + Wy(k1,ipntr)*Wy(k1,jpntr)
               temp4 = temp4 + Ws(k1,ipntr)*Ws(k1,jpntr)
            enddo
            Wn1(iy,jy) = Wn1(iy,jy) + temp1 - temp3
            Wn1(is,js) = Wn1(is,js) - temp2 + temp4
            jpntr = mod(jpntr,m) + 1
         enddo
         ipntr = mod(ipntr,m) + 1
      enddo

!       modify the old parts in block (2,1).
      ipntr = Head
      do is = m + 1 , m + upcl
         jpntr = Head
         do jy = 1 , upcl
            temp1 = zero
            temp3 = zero
            do k = 1 , Nenter
               k1 = Indx2(k)
               temp1 = temp1 + Ws(k1,ipntr)*Wy(k1,jpntr)
            enddo
            do k = Ileave , n
               k1 = Indx2(k)
               temp3 = temp3 + Ws(k1,ipntr)*Wy(k1,jpntr)
            enddo
            if ( is<=jy+m ) then
               Wn1(is,jy) = Wn1(is,jy) + temp1 - temp3
            else
               Wn1(is,jy) = Wn1(is,jy) - temp1 + temp3
            endif
            jpntr = mod(jpntr,m) + 1
         enddo
         ipntr = mod(ipntr,m) + 1
      enddo

!     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
!                                     [-L_a +R_z        S'AA'S*theta]

      m2 = 2*m
      do iy = 1 , Col
         is = Col + iy
         is1 = m + iy
         do jy = 1 , iy
            js = Col + jy
            js1 = m + jy
            Wn(jy,iy) = Wn1(iy,jy)/Theta
            Wn(js,is) = Wn1(is1,js1)*Theta
         enddo
         do jy = 1 , iy - 1
            Wn(jy,is) = -Wn1(is1,jy)
         enddo
         do jy = iy , Col
            Wn(jy,is) = Wn1(is1,jy)
         enddo
         Wn(iy,iy) = Wn(iy,iy) + Sy(iy,iy)
      enddo

!     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
!                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ]

!        first Cholesky factor (1,1) block of wn to get LL'
!                          with L' stored in the upper triangle of wn.
      call dpofa(Wn,m2,Col,Info)
      if ( Info/=0 ) then
         Info = -1
         return
      endif
!        then form L^-1(-L_a'+R_z') in the (1,2) block.
      col2 = 2*Col
      do js = Col + 1 , col2
         call dtrsl(Wn,m2,Col,Wn(1,js),11,Info)
      enddo

!     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
!        upper triangle of (2,2) block of wn.


      do is = Col + 1 , col2
         do js = is , col2
            Wn(is,js) = Wn(is,js) + ddot(Col,Wn(1,is),1,Wn(1,js),1)
         enddo
      enddo

!     Cholesky factorization of (2,2) block of wn.

      call dpofa(Wn(Col+1,Col+1),m2,Col,Info)
      if ( Info/=0 ) then
         Info = -2
         return
      endif

      end subroutine formk

      subroutine formt(m,Wt,Sy,Ss,Col,Theta,Info)
      implicit none

      integer m , Col , Info
      real(wp) Theta , Wt(m,m) , Sy(m,m) , Ss(m,m)

!     ************
!
!     Subroutine formt
!
!       This subroutine forms the upper half of the pos. def. and symm.
!         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
!         of the array wt, and performs the Cholesky factorization of T
!         to produce J*J', with J' stored in the upper triangle of wt.
!
!     Subprograms called:
!
!       Linpack ... dpofa.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer i , j , k , k1
      real(wp) ddum


!     Form the upper half of  T = theta*SS + L*D^(-1)*L',
!        store T in the upper triangle of the array wt.

      do j = 1 , Col
         Wt(1,j) = Theta*Ss(1,j)
      enddo
      do i = 2 , Col
         do j = i , Col
            k1 = min(i,j) - 1
            ddum = zero
            do k = 1 , k1
               ddum = ddum + Sy(i,k)*Sy(j,k)/Sy(k,k)
            enddo
            Wt(i,j) = ddum + Theta*Ss(i,j)
         enddo
      enddo

!     Cholesky factorize T to J*J' with
!        J' stored in the upper triangle of wt.

      call dpofa(Wt,m,Col,Info)
      if ( Info/=0 ) Info = -3

      end subroutine formt

      subroutine freev(n,Nfree,Index,Nenter,Ileave,Indx2,Iwhere,Wrk,    &
                     & Updatd,Cnstnd,Iprint,Iter)
      implicit none

      integer n , Nfree , Nenter , Ileave , Iprint , Iter , Index(n) ,  &
            & Indx2(n) , Iwhere(n)
      logical Wrk , Updatd , Cnstnd

!     ************
!
!     Subroutine freev
!
!     This subroutine counts the entering and leaving variables when
!       iter > 0, and finds the index set of free and active variables
!       at the GCP.
!
!     cnstnd is a logical variable indicating whether bounds are present
!
!     index is an integer array of dimension n
!       for i=1,...,nfree, index(i) are the indices of free variables
!       for i=nfree+1,...,n, index(i) are the indices of bound variables
!       On entry after the first iteration, index gives
!         the free variables at the previous iteration.
!       On exit it gives the free variables based on the determination
!         in cauchy using the array iwhere.
!
!     indx2 is an integer array of dimension n
!       On entry indx2 is unspecified.
!       On exit with iter>0, indx2 indicates which variables
!          have changed status since the previous iteration.
!       For i= 1,...,nenter, indx2(i) have changed from bound to free.
!       For i= ileave+1,...,n, indx2(i) have changed from free to bound.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer iact , i , k

      Nenter = 0
      Ileave = n + 1
      if ( Iter>0 .and. Cnstnd ) then
!                           count the entering and leaving variables.
         do i = 1 , Nfree
            k = Index(i)

!            write(6,*) ' k  = index(i) ', k
!            write(6,*) ' index = ', i

            if ( Iwhere(k)>0 ) then
               Ileave = Ileave - 1
               Indx2(Ileave) = k
               if ( Iprint>=100 ) write (6,*) 'Variable ' , k ,       &
                   &' leaves the set of free variables'
            endif
         enddo
         do i = 1 + Nfree , n
            k = Index(i)
            if ( Iwhere(k)<=0 ) then
               Nenter = Nenter + 1
               Indx2(Nenter) = k
               if ( Iprint>=100 ) write (6,*) 'Variable ' , k ,       &
                   &' enters the set of free variables'
            endif
         enddo
         if ( Iprint>=99 ) write (6,*) n + 1 - Ileave ,               &
                                   &' variables leave; ' , Nenter ,     &
                                   &' variables enter'
      endif
      Wrk = (Ileave<n+1) .or. (Nenter>0) .or. Updatd

!     Find the index set of free and active variables at the GCP.

      Nfree = 0
      iact = n + 1
      do i = 1 , n
         if ( Iwhere(i)<=0 ) then
            Nfree = Nfree + 1
            Index(Nfree) = i
         else
            iact = iact - 1
            Index(iact) = i
         endif
      enddo
      if ( Iprint>=99 ) write (6,*) Nfree ,                           &
                                     &' variables are free at GCP ' ,   &
                                    & Iter + 1

      end subroutine freev

      subroutine hpsolb(n,t,Iorder,Iheap)
      implicit none

      integer Iheap , n , Iorder(n)
      real(wp) t(n)

!     ************
!
!     Subroutine hpsolb
!
!     This subroutine sorts out the least element of t, and puts the
!       remaining elements of t in a heap.
!
!     n is an integer variable.
!       On entry n is the dimension of the arrays t and iorder.
!       On exit n is unchanged.
!
!     t is a real(wp) array of dimension n.
!       On entry t stores the elements to be sorted,
!       On exit t(n) stores the least elements of t, and t(1) to t(n-1)
!         stores the remaining elements in the form of a heap.
!
!     iorder is an integer array of dimension n.
!       On entry iorder(i) is the index of t(i).
!       On exit iorder(i) is still the index of t(i), but iorder may be
!         permuted in accordance with t.
!
!     iheap is an integer variable specifying the task.
!       On entry iheap should be set as follows:
!         iheap == 0 if t(1) to t(n) is not in the form of a heap,
!         iheap /= 0 if otherwise.
!       On exit iheap is unchanged.
!
!
!     References:
!       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!     ************

      integer i , j , k , indxin , indxou
      real(wp) ddum , out

      if ( Iheap==0 ) then

!        Rearrange the elements t(1) to t(n) to form a heap.

         do k = 2 , n
            ddum = t(k)
            indxin = Iorder(k)

!           Add ddum to the heap.
            i = k
            do
               if ( i>1 ) then
                  j = i/2
                  if ( ddum<t(j) ) then
                     t(i) = t(j)
                     Iorder(i) = Iorder(j)
                     i = j
                     cycle
                  endif
               endif
               exit
            end do
            t(i) = ddum
            Iorder(i) = indxin
         enddo
      endif

!     Assign to 'out' the value of t(1), the least member of the heap,
!        and rearrange the remaining members to form a heap as
!        elements 1 to n-1 of t.

      if ( n>1 ) then
         i = 1
         out = t(1)
         indxou = Iorder(1)
         ddum = t(n)
         indxin = Iorder(n)

!        Restore the heap
         do
            j = i + i
            if ( j<=n-1 ) then
               if ( t(j+1)<t(j) ) j = j + 1
               if ( t(j)<ddum ) then
                  t(i) = t(j)
                  Iorder(i) = Iorder(j)
                  i = j
                  cycle
               endif
            endif
            exit
         end do
         t(i) = ddum
         Iorder(i) = indxin

!     Put the least member in t(n).

         t(n) = out
         Iorder(n) = indxou
      endif

      end subroutine hpsolb

      subroutine lnsrlb(n,l,u,Nbd,x,f,Fold,Gd,Gdold,g,d,r,t,z,Stp,Dnorm,&
                      & Dtd,Xstep,Stpmx,Iter,Ifun,Iback,Nfgv,Info,Task, &
                      & Boxed,Cnstnd,Csave,Isave,Dsave)
      implicit none

      character*60 Task , Csave
      logical Boxed , Cnstnd
      integer n , Iter , Ifun , Iback , Nfgv , Info , Nbd(n) , Isave(2)
      real(wp) f , Fold , Gd , Gdold , Stp , Dnorm , Dtd ,      &
                     & Xstep , Stpmx , x(n) , l(n) , u(n) , g(n) ,      &
                     & d(n) , r(n) , t(n) , z(n) , Dsave(13)
!     **********
!
!     Subroutine lnsrlb
!
!     This subroutine calls subroutine dcsrch from the Minpack2 library
!       to perform the line search.  Subroutine dscrch is safeguarded so
!       that all trial points lie within the feasible region.
!
!     Subprograms called:
!
!       Minpack2 Library ... dcsrch.
!
!       Linpack ... dtrsl, ddot.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     **********

      integer i
      real(wp) a1 , a2

      real(wp),parameter :: big  = 1.0e+10_wp
      real(wp),parameter :: ftol = 1.0e-3_wp
      real(wp),parameter :: gtol = 0.9_wp
      real(wp),parameter :: xtol = 0.1_wp

      if ( Task(1:5)/='FG_LN' ) then

         Dtd = ddot(n,d,1,d,1)
         Dnorm = sqrt(Dtd)

   !     Determine the maximum step length.

         Stpmx = big
         if ( Cnstnd ) then
            if ( Iter==0 ) then
               Stpmx = one
            else
               do i = 1 , n
                  a1 = d(i)
                  if ( Nbd(i)/=0 ) then
                     if ( a1<zero .and. Nbd(i)<=2 ) then
                        a2 = l(i) - x(i)
                        if ( a2>=zero ) then
                           Stpmx = zero
                        elseif ( a1*Stpmx<a2 ) then
                           Stpmx = a2/a1
                        endif
                     elseif ( a1>zero .and. Nbd(i)>=2 ) then
                        a2 = u(i) - x(i)
                        if ( a2<=zero ) then
                           Stpmx = zero
                        elseif ( a1*Stpmx>a2 ) then
                           Stpmx = a2/a1
                        endif
                     endif
                  endif
               enddo
            endif
         endif

         if ( Iter==0 .and. .not.Boxed ) then
            Stp = min(one/Dnorm,Stpmx)
         else
            Stp = one
         endif

         call dcopy(n,x,1,t,1)
         call dcopy(n,g,1,r,1)
         Fold = f
         Ifun = 0
         Iback = 0
         Csave = 'START'

      end if

      Gd = ddot(n,g,1,d,1)
      if ( Ifun==0 ) then
         Gdold = Gd
         if ( Gd>=zero ) then
!                               the directional derivative >=0.
!                               Line search is impossible.
            write (6,*) ' ascent direction in projection gd = ' , Gd
            Info = -4
            return
         endif
      endif

      call dcsrch(f,Gd,Stp,ftol,gtol,xtol,zero,Stpmx,Csave,Isave,Dsave)

      Xstep = Stp*Dnorm
      if ( Csave(1:4)/='CONV' .and. Csave(1:4)/='WARN' ) then
         Task = 'FG_LNSRCH'
         Ifun = Ifun + 1
         Nfgv = Nfgv + 1
         Iback = Ifun - 1
         if ( Stp==one ) then
            call dcopy(n,z,1,x,1)
         else
            do i = 1 , n
               x(i) = Stp*d(i) + t(i)
            enddo
         endif
      else
         Task = 'NEW_X'
      endif

      end subroutine lnsrlb

      subroutine matupd(n,m,Ws,Wy,Sy,Ss,d,r,Itail,Iupdat,Col,Head,Theta,&
                      & Rr,Dr,Stp,Dtd)
      implicit none

      integer n , m , Itail , Iupdat , Col , Head
      real(wp) Theta , Rr , Dr , Stp , Dtd , d(n) , r(n) ,      &
                     & Ws(n,m) , Wy(n,m) , Sy(m,m) , Ss(m,m)

!     ************
!
!     Subroutine matupd
!
!       This subroutine updates matrices WS and WY, and forms the
!         middle matrix in B.
!
!     Subprograms called:
!
!       Linpack ... dcopy, ddot.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer j , pointr

!     Set pointers for matrices WS and WY.

      if ( Iupdat<=m ) then
         Col = Iupdat
         Itail = mod(Head+Iupdat-2,m) + 1
      else
         Itail = mod(Itail,m) + 1
         Head = mod(Head,m) + 1
      endif

!     Update matrices WS and WY.

      call dcopy(n,d,1,Ws(1,Itail),1)
      call dcopy(n,r,1,Wy(1,Itail),1)

!     Set theta=yy/ys.

      Theta = Rr/Dr

!     Form the middle matrix in B.

!        update the upper triangle of SS,
!                                         and the lower triangle of SY:
      if ( Iupdat>m ) then
!                              move old information
         do j = 1 , Col - 1
            call dcopy(j,Ss(2,j+1),1,Ss(1,j),1)
            call dcopy(Col-j,Sy(j+1,j+1),1,Sy(j,j),1)
         enddo
      endif
!        add new information: the last row of SY
!                                             and the last column of SS:
      pointr = Head
      do j = 1 , Col - 1
         Sy(Col,j) = ddot(n,d,1,Wy(1,pointr),1)
         Ss(j,Col) = ddot(n,Ws(1,pointr),1,d,1)
         pointr = mod(pointr,m) + 1
      enddo
      if ( Stp==one ) then
         Ss(Col,Col) = Dtd
      else
         Ss(Col,Col) = Stp*Stp*Dtd
      endif
      Sy(Col,Col) = Dr

      end subroutine matupd

      subroutine prn1lb(n,m,l,u,x,Iprint,Itfile,Epsmch)
      implicit none

      integer n , m , Iprint , Itfile
      real(wp) Epsmch , x(n) , l(n) , u(n)

!     ************
!
!     Subroutine prn1lb
!
!     This subroutine prints the input data, initial point, upper and
!       lower bounds of each variable, machine precision, as well as
!       the headings of the output.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer i

      if ( Iprint>=0 ) then
         write (6,99001) Epsmch
99001    format ('RUNNING THE L-BFGS-B CODE',/,/,'           * * *',/,/,&
                &'Machine precision =',1p,d10.3)
         write (6,*) 'N = ' , n , '    M = ' , m
         if ( Iprint>=1 ) then
            write (Itfile,99002) Epsmch
99002       format ('RUNNING THE L-BFGS-B CODE',/,/,                    &
                   &'it    = iteration number',/,                       &
                   &'nf    = number of function evaluations',/,         &
         &'nseg  = number of segments explored during the Cauchy search'&
        & ,/,                                                           &
      &'nact  = number of active bounds at the generalized Cauchy point'&
     & ,/,                                                              &
      &'sub   = manner in which the subspace minimization terminated:', &
     & /,'        con = converged, bnd = a bound was reached',/,        &
      &'itls  = number of iterations performed in the line search',/,   &
      &'stepl = step length used',/,                                    &
      &'tstep = norm of the displacement (total step)',/,               &
      &'projg = norm of the projected gradient',/,                      &
      &'f     = function value',/,/,'           * * *',/,/,             &
      &'Machine precision =',1p,d10.3)
            write (Itfile,*) 'N = ' , n , '    M = ' , m
            write (Itfile,99003)
99003       format (/,3x,'it',3x,'nf',2x,'nseg',2x,'nact',2x,'sub',2x,  &
                   &'itls',2x,'stepl',4x,'tstep',5x,'projg',8x,'f')
            if ( Iprint>100 ) then
               write (6,99004) 'L =' , (l(i),i=1,n)
               write (6,99004) 'X0 =' , (x(i),i=1,n)
               write (6,99004) 'U =' , (u(i),i=1,n)
            endif
         endif
      endif

      continue

99004 format (/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))

      end subroutine prn1lb

      subroutine prn2lb(n,x,f,g,Iprint,Itfile,Iter,Nfgv,Nact,Sbgnrm,    &
                      & Nseg,Word,Iword,Iback,Stp,Xstep)
      implicit none

      character*3 Word
      integer n , Iprint , Itfile , Iter , Nfgv , Nact , Nseg , Iword , &
            & Iback
      real(wp) f , Sbgnrm , Stp , Xstep , x(n) , g(n)

!     ************
!
!     Subroutine prn2lb
!
!     This subroutine prints out new information after a successful
!       line search.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer i , imod

!           'word' records the status of subspace solutions.
      if ( Iword==0 ) then
!                            the subspace minimization converged.
         Word = 'con'
      elseif ( Iword==1 ) then
!                          the subspace minimization stopped at a bound.
         Word = 'bnd'
      elseif ( Iword==5 ) then
!                             the truncated Newton step has been used.
         Word = 'TNT'
      else
         Word = '---'
      endif
      if ( Iprint>=99 ) then
         write (6,*) 'LINE SEARCH' , Iback , ' times; norm of step = ' ,&
                   & Xstep
         write (6,99003) Iter , f , Sbgnrm
         if ( Iprint>100 ) then
            write (6,99002) 'X =' , (x(i),i=1,n)
            write (6,99002) 'G =' , (g(i),i=1,n)
         endif
      elseif ( Iprint>0 ) then
         imod = mod(Iter,Iprint)
         if ( imod==0 ) write (6,99003) Iter , f , Sbgnrm
      endif
      if ( Iprint>=1 ) write (Itfile,99001) Iter , Nfgv , Nseg ,      &
                              & Nact , Word , Iback , Stp , Xstep ,     &
                              & Sbgnrm , f
99001 format (2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1p,2(1x,d10.3)&
            & )

      continue

99002 format (/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
99003 format (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,   &
            & d12.5)

      end subroutine prn2lb

      subroutine prn3lb(n,x,f,Task,Iprint,Info,Itfile,Iter,Nfgv,Nintol, &
                      & Nskip,Nact,Sbgnrm,Time,Nseg,Word,Iback,Stp,     &
                      & Xstep,k,Cachyt,Sbtime,Lnscht)
      implicit none

      character*60 Task
      character*3 Word
      integer n , Iprint , Info , Itfile , Iter , Nfgv , Nintol ,       &
            & Nskip , Nact , Nseg , Iback , k
      real(wp) f , Sbgnrm , Time , Stp , Xstep , Cachyt ,       &
                     & Sbtime , Lnscht , x(n)

!     ************
!
!     Subroutine prn3lb
!
!     This subroutine prints out information when either a built-in
!       convergence test is satisfied or when an error message is
!       generated.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer i

      if ( Task(1:5)/='ERROR' ) then

         if ( Iprint>=0 ) then
            write (6,99001)
   99001    format (/,'           * * *',/,/,                              &
                  &'Tit   = total number of iterations',/,                &
                  &'Tnf   = total number of function evaluations',/,      &
                  &'Tnint = total number of segments explored during',    &
                  &' Cauchy searches',/,                                  &
                  &'Skip  = number of BFGS updates skipped',/,            &
                  &'Nact  = number of active bounds at final generalized',&
                  &' Cauchy point',/,                                     &
                  &'Projg = norm of the final projected gradient',/,      &
                  &'F     = final function value',/,/,'           * * *')
            write (6,99002)
   99002    format (/,3x,'N',4x,'Tit',5x,'Tnf',2x,'Tnint',2x,'Skip',2x,    &
                  &'Nact',5x,'Projg',8x,'F')
            write (6,99003) n , Iter , Nfgv , Nintol , Nskip , Nact ,      &
                        & Sbgnrm , f
   99003    format (i5,2(1x,i6),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d10.3))
            if ( Iprint>=100 ) then
               write (6,99004) 'X =' , (x(i),i=1,n)

   99004       format (/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
            endif
            if ( Iprint>=1 ) write (6,*) ' F =' , f
         endif

      end if

      if ( Iprint>=0 ) then
         write (6,99008) Task
         if ( Info/=0 ) then
            if ( Info==-1 ) write (6,99009)
            if ( Info==-2 ) write (6,99010)
            if ( Info==-3 ) write (6,99011)
            if ( Info==-4 ) write (6,99012)
            if ( Info==-5 ) write (6,99013)
            if ( Info==-6 ) write (6,*) ' Input nbd(' , k ,           &
                                    &') is invalid.'
            if ( Info==-7 ) write (6,*) ' l(' , k , ') > u(' , k ,    &
                                    &').  No feasible solution.'
            if ( Info==-8 ) write (6,99014)
            if ( Info==-9 ) write (6,99015)
         endif
         if ( Iprint>=1 ) write (6,99005) Cachyt , Sbtime , Lnscht
99005    format (/,' Cauchy                time',1p,e10.3,' seconds.',  &
                &/' Subspace minimization time',1p,e10.3,' seconds.',   &
                &/' Line search           time',1p,e10.3,' seconds.')
         write (6,99007) Time
         if ( Iprint>=1 ) then
            if ( Info==-4 .or. Info==-9 ) then
               write (Itfile,99006) Iter , Nfgv , Nseg , Nact , Word ,  &
                                  & Iback , Stp , Xstep
99006          format (2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6x,  &
                      &'-',10x,'-')
            endif
            write (Itfile,99008) Task
            if ( Info/=0 ) then
               if ( Info==-1 ) write (Itfile,99009)
               if ( Info==-2 ) write (Itfile,99010)
               if ( Info==-3 ) write (Itfile,99011)
               if ( Info==-4 ) write (Itfile,99012)
               if ( Info==-5 ) write (Itfile,99013)
               if ( Info==-8 ) write (Itfile,99014)
               if ( Info==-9 ) write (Itfile,99015)
            endif
            write (Itfile,99007) Time
         endif
      endif

      continue
99007 format (/,' Total User time',1p,e10.3,' seconds.',/)
99008 format (/,a60)
99009 format (/,                                                        &
      &' Matrix in 1st Cholesky factorization in formk is not Pos. Def.'&
     & )
99010 format (/,                                                        &
      &' Matrix in 2st Cholesky factorization in formk is not Pos. Def.'&
     & )
99011 format (/,                                                        &
      &' Matrix in the Cholesky factorization in formt is not Pos. Def.'&
     & )
99012 format (/,' Derivative >= 0, backtracking line search impossible.'&
            & ,/,'   Previous x, f and g restored.',/,                  &
        &' Possible causes: 1 error in function or gradient evaluation;'&
       & ,/,'                  2 rounding errors dominate computation.')
99013 format (/,' Warning:  more than 10 function and gradient',/,      &
             &'   evaluations in the last line search.  Termination',/, &
             &'   may possibly be caused by a bad search direction.')
99014 format (/,' The triangular system is singular.')
99015 format (/,                                                        &
       &' Line search cannot locate an adequate point after 20 function'&
      & ,/,'  and gradient evaluations.  Previous x, f and g restored.',&
      & /,                                                              &
       &' Possible causes: 1 error in function or gradient evaluation;',&
      & /,'                  2 rounding error dominate computation.')

      end subroutine prn3lb

      subroutine projgr(n,l,u,Nbd,x,g,Sbgnrm)
      implicit none

      integer n , Nbd(n)
      real(wp) Sbgnrm , x(n) , l(n) , u(n) , g(n)

!     ************
!
!     Subroutine projgr
!
!     This subroutine computes the infinity norm of the projected
!       gradient.
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer i
      real(wp) gi

      Sbgnrm = zero
      do i = 1 , n
         gi = g(i)
         if ( Nbd(i)/=0 ) then
            if ( gi<zero ) then
               if ( Nbd(i)>=2 ) gi = max((x(i)-u(i)),gi)
            else
               if ( Nbd(i)<=2 ) gi = min((x(i)-l(i)),gi)
            endif
         endif
         Sbgnrm = max(Sbgnrm,abs(gi))
      enddo

      end subroutine projgr

      subroutine subsm(n,m,Nsub,Ind,l,u,Nbd,x,d,Xp,Ws,Wy,Theta,Xx,Gg,   &
                     & Col,Head,Iword,Wv,Wn,Iprint,Info)
      implicit none

      integer n , m , Nsub , Col , Head , Iword , Iprint , Info ,       &
            & Ind(Nsub) , Nbd(n)
      real(wp) Theta , l(n) , u(n) , x(n) , d(n) , Xp(n) ,      &
                     & Xx(n) , Gg(n) , Ws(n,m) , Wy(n,m) , Wv(2*m) ,    &
                     & Wn(2*m,2*m)

!     **********************************************************************
!
!     This routine contains the major changes in the updated version.
!     The changes are described in the accompanying paper
!
!      Jose Luis Morales, Jorge Nocedal
!      "Remark On Algorithm 788: L-BFGS-B: Fortran Subroutines for Large-Scale
!       Bound Constrained Optimization". Decemmber 27, 2010.
!
!             J.L. Morales  Departamento de Matematicas,
!                           Instituto Tecnologico Autonomo de Mexico
!                           Mexico D.F.
!
!             J, Nocedal    Department of Electrical Engineering and
!                           Computer Science.
!                           Northwestern University. Evanston, IL. USA
!
!                           January 17, 2011
!
!      **********************************************************************
!
!
!     Subroutine subsm
!
!     Given xcp, l, u, r, an index set that specifies
!       the active set at xcp, and an l-BFGS matrix B
!       (in terms of WY, WS, SY, WT, head, col, and theta),
!       this subroutine computes an approximate solution
!       of the subspace problem
!
!       (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
!
!             subject to l<=x<=u
!                       x_i=xcp_i for all i in A(xcp)
!
!       along the subspace unconstrained Newton direction
!
!          d = -(Z'BZ)^(-1) r.
!
!       The formula for the Newton direction, given the L-BFGS matrix
!       and the Sherman-Morrison formula, is
!
!          d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
!
!       where
!                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                     [L_a -R_z           theta*S'AA'S ]
!
!     Note that this procedure for computing d differs
!     from that described in [1]. One can show that the matrix K is
!     equal to the matrix M^[-1]N in that paper.
!
!     n is an integer variable.
!       On entry n is the dimension of the problem.
!       On exit n is unchanged.
!
!     m is an integer variable.
!       On entry m is the maximum number of variable metric corrections
!         used to define the limited memory matrix.
!       On exit m is unchanged.
!
!     nsub is an integer variable.
!       On entry nsub is the number of free variables.
!       On exit nsub is unchanged.
!
!     ind is an integer array of dimension nsub.
!       On entry ind specifies the coordinate indices of free variables.
!       On exit ind is unchanged.
!
!     l is a real(wp) array of dimension n.
!       On entry l is the lower bound of x.
!       On exit l is unchanged.
!
!     u is a real(wp) array of dimension n.
!       On entry u is the upper bound of x.
!       On exit u is unchanged.
!
!     nbd is a integer array of dimension n.
!       On entry nbd represents the type of bounds imposed on the
!         variables, and must be specified as follows:
!         nbd(i)=0 if x(i) is unbounded,
!                1 if x(i) has only a lower bound,
!                2 if x(i) has both lower and upper bounds, and
!                3 if x(i) has only an upper bound.
!       On exit nbd is unchanged.
!
!     x is a real(wp) array of dimension n.
!       On entry x specifies the Cauchy point xcp.
!       On exit x(i) is the minimizer of Q over the subspace of
!                                                        free variables.
!
!     d is a real(wp) array of dimension n.
!       On entry d is the reduced gradient of Q at xcp.
!       On exit d is the Newton direction of Q.
!
!    xp is a real(wp) array of dimension n.
!       used to safeguard the projected Newton direction
!
!    xx is a real(wp) array of dimension n
!       On entry it holds the current iterate
!       On output it is unchanged

!    gg is a real(wp) array of dimension n
!       On entry it holds the gradient at the current iterate
!       On output it is unchanged
!
!     ws and wy are real(wp) arrays;
!     theta is a real(wp) variable;
!     col is an integer variable;
!     head is an integer variable.
!       On entry they store the information defining the
!                                          limited memory BFGS matrix:
!         ws(n,m) stores S, a set of s-vectors;
!         wy(n,m) stores Y, a set of y-vectors;
!         theta is the scaling factor specifying B_0 = theta I;
!         col is the number of variable metric corrections stored;
!         head is the location of the 1st s- (or y-) vector in S (or Y).
!       On exit they are unchanged.
!
!     iword is an integer variable.
!       On entry iword is unspecified.
!       On exit iword specifies the status of the subspace solution.
!         iword = 0 if the solution is in the box,
!                 1 if some bound is encountered.
!
!     wv is a real(wp) working array of dimension 2m.
!
!     wn is a real(wp) array of dimension 2m x 2m.
!       On entry the upper triangle of wn stores the LEL^T factorization
!         of the indefinite matrix
!
!              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
!                  [L_a -R_z           theta*S'AA'S ]
!                                                    where E = [-I  0]
!                                                              [ 0  I]
!       On exit wn is unchanged.
!
!     iprint is an INTEGER variable that must be set by the user.
!       It controls the frequency and type of output generated:
!        iprint<0    no output is generated;
!        iprint=0    print only one line at the last iteration;
!        0<iprint<99 print also f and |proj g| every iprint iterations;
!        iprint=99   print details of every iteration except n-vectors;
!        iprint=100  print also the changes of active set and final x;
!        iprint>100  print details of every iteration including x and g;
!       When iprint > 0, the file iterate.dat will be created to
!                        summarize the iteration.
!
!     info is an integer variable.
!       On entry info is unspecified.
!       On exit info = 0       for normal return,
!                    = nonzero for abnormal return
!                                  when the matrix K is ill-conditioned.
!
!     Subprograms called:
!
!       Linpack dtrsl.
!
!
!     References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, "A limited
!       memory algorithm for bound constrained optimization",
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!
!
!                           *  *  *
!
!     NEOS, November 1994. (Latest revision June 1996.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!
!
!     ************

      integer pointr , m2 , col2 , ibd , jy , js , i , j , k
      real(wp) alpha , xk , dk , temp1 , temp2
!
      real(wp) dd_p

      if ( Nsub<=0 ) return
      if ( Iprint>=99 ) write (6,99001)

99001 format (/,'----------------SUBSM entered-----------------',/)

!     Compute wv = W'Zd.

      pointr = Head
      do i = 1 , Col
         temp1 = zero
         temp2 = zero
         do j = 1 , Nsub
            k = Ind(j)
            temp1 = temp1 + Wy(k,pointr)*d(j)
            temp2 = temp2 + Ws(k,pointr)*d(j)
         enddo
         Wv(i) = temp1
         Wv(Col+i) = Theta*temp2
         pointr = mod(pointr,m) + 1
      enddo

!     Compute wv:=K^(-1)wv.

      m2 = 2*m
      col2 = 2*Col
      call dtrsl(Wn,m2,col2,Wv,11,Info)
      if ( Info/=0 ) return
      do i = 1 , Col
         Wv(i) = -Wv(i)
      enddo
      call dtrsl(Wn,m2,col2,Wv,01,Info)
      if ( Info/=0 ) return

!     Compute d = (1/theta)d + (1/theta**2)Z'W wv.

      pointr = Head
      do jy = 1 , Col
         js = Col + jy
         do i = 1 , Nsub
            k = Ind(i)
            d(i) = d(i) + Wy(k,pointr)*Wv(jy)/Theta + Ws(k,pointr)      &
                 & *Wv(js)
         enddo
         pointr = mod(pointr,m) + 1
      enddo

      call dscal(Nsub,one/Theta,d,1)
!
!-----------------------------------------------------------------
!     Let us try the projection, d is the Newton direction

      Iword = 0

      call dcopy(n,x,1,Xp,1)
!
      do i = 1 , Nsub
         k = Ind(i)
         dk = d(i)
         xk = x(k)
         if ( Nbd(k)/=0 ) then
!
            if ( Nbd(k)==1 ) then          ! lower bounds only
               x(k) = max(l(k),xk+dk)
               if ( x(k)==l(k) ) Iword = 1
            else
!
               if ( Nbd(k)==2 ) then       ! upper and lower bounds
                  xk = max(l(k),xk+dk)
                  x(k) = min(u(k),xk)
                  if ( x(k)==l(k) .or. x(k)==u(k) ) Iword = 1
               else
!
                  if ( Nbd(k)==3 ) then    ! upper bounds only
                     x(k) = min(u(k),xk+dk)
                     if ( x(k)==u(k) ) Iword = 1
                  endif
               endif
            endif
!
         else                                ! free variables
            x(k) = xk + dk
         endif
      enddo

      main : block

         if ( Iword==0 ) exit main
   !
   !     check sign of the directional derivative
   !
         dd_p = zero
         do i = 1 , n
            dd_p = dd_p + (x(i)-Xx(i))*Gg(i)
         enddo
         if ( dd_p<=zero ) exit main

         call dcopy(n,Xp,1,x,1)
         write (6,*) ' Positive dir derivative in projection '
         write (6,*) ' Using the backtracking step '
   !
   !-----------------------------------------------------------------
   !
         alpha = one
         temp1 = alpha
         ibd = 0
         do i = 1 , Nsub
            k = Ind(i)
            dk = d(i)
            if ( Nbd(k)/=0 ) then
               if ( dk<zero .and. Nbd(k)<=2 ) then
                  temp2 = l(k) - x(k)
                  if ( temp2>=zero ) then
                     temp1 = zero
                  elseif ( dk*alpha<temp2 ) then
                     temp1 = temp2/dk
                  endif
               elseif ( dk>zero .and. Nbd(k)>=2 ) then
                  temp2 = u(k) - x(k)
                  if ( temp2<=zero ) then
                     temp1 = zero
                  elseif ( dk*alpha>temp2 ) then
                     temp1 = temp2/dk
                  endif
               endif
               if ( temp1<alpha ) then
                  alpha = temp1
                  ibd = i
               endif
            endif
         enddo

         if ( alpha<one ) then
            dk = d(ibd)
            k = Ind(ibd)
            if ( dk>zero ) then
               x(k) = u(k)
               d(ibd) = zero
            elseif ( dk<zero ) then
               x(k) = l(k)
               d(ibd) = zero
            endif
         endif
         do i = 1 , Nsub
            k = Ind(i)
            x(k) = x(k) + alpha*d(i)
         enddo

      end block main

      if ( Iprint>=99 ) write (6,99002)
99002 format (/,'----------------exit SUBSM --------------------',/)

      end subroutine subsm

      subroutine dcsrch(f,g,Stp,Ftol,Gtol,Xtol,Stpmin,Stpmax,Task,Isave,&
                      & Dsave)
      implicit none

      character*(*) Task
      integer Isave(2)
      real(wp) f , g , Stp , Ftol , Gtol , Xtol , Stpmin ,      &
                     & Stpmax
      real(wp) Dsave(13)
!     **********
!
!     Subroutine dcsrch
!
!     This subroutine finds a step that satisfies a sufficient
!     decrease condition and a curvature condition.
!
!     Each call of the subroutine updates an interval with
!     endpoints stx and sty. The interval is initially chosen
!     so that it contains a minimizer of the modified function
!
!           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
!
!     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
!     interval is chosen so that it contains a minimizer of f.
!
!     The algorithm is designed to find a step that satisfies
!     the sufficient decrease condition
!
!           f(stp) <= f(0) + ftol*stp*f'(0),
!
!     and the curvature condition
!
!           abs(f'(stp)) <= gtol*abs(f'(0)).
!
!     If ftol is less than gtol and if, for example, the function
!     is bounded below, then there is always a step which satisfies
!     both conditions.
!
!     If no step can be found that satisfies both conditions, then
!     the algorithm stops with a warning. In this case stp only
!     satisfies the sufficient decrease condition.
!
!     A typical invocation of dcsrch has the following outline:
!
!     task = 'START'
!  10 continue
!        call dcsrch( ... )
!        if (task == 'FG') then
!           Evaluate the function and the gradient at stp
!           goto 10
!           end if
!
!     NOTE: The user must no alter work arrays between calls.
!
!     The subroutine statement is
!
!        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
!                          task,isave,dsave)
!     where
!
!       f is a real(wp) variable.
!         On initial entry f is the value of the function at 0.
!            On subsequent entries f is the value of the
!            function at stp.
!         On exit f is the value of the function at stp.
!
!       g is a real(wp) variable.
!         On initial entry g is the derivative of the function at 0.
!            On subsequent entries g is the derivative of the
!            function at stp.
!         On exit g is the derivative of the function at stp.
!
!       stp is a real(wp) variable.
!         On entry stp is the current estimate of a satisfactory
!            step. On initial entry, a positive initial estimate
!            must be provided.
!         On exit stp is the current estimate of a satisfactory step
!            if task = 'FG'. If task = 'CONV' then stp satisfies
!            the sufficient decrease and curvature condition.
!
!       ftol is a real(wp) variable.
!         On entry ftol specifies a nonnegative tolerance for the
!            sufficient decrease condition.
!         On exit ftol is unchanged.
!
!       gtol is a real(wp) variable.
!         On entry gtol specifies a nonnegative tolerance for the
!            curvature condition.
!         On exit gtol is unchanged.
!
!       xtol is a real(wp) variable.
!         On entry xtol specifies a nonnegative relative tolerance
!            for an acceptable step. The subroutine exits with a
!            warning if the relative difference between sty and stx
!            is less than xtol.
!         On exit xtol is unchanged.
!
!       stpmin is a real(wp) variable.
!         On entry stpmin is a nonnegative lower bound for the step.
!         On exit stpmin is unchanged.
!
!       stpmax is a real(wp) variable.
!         On entry stpmax is a nonnegative upper bound for the step.
!         On exit stpmax is unchanged.
!
!       task is a character variable of length at least 60.
!         On initial entry task must be set to 'START'.
!         On exit task indicates the required action:
!
!            If task(1:2) = 'FG' then evaluate the function and
!            derivative at stp and call dcsrch again.
!
!            If task(1:4) = 'CONV' then the search is successful.
!
!            If task(1:4) = 'WARN' then the subroutine is not able
!            to satisfy the convergence conditions. The exit value of
!            stp contains the best point found during the search.
!
!            If task(1:5) = 'ERROR' then there is an error in the
!            input arguments.
!
!         On exit with convergence, a warning or an error, the
!            variable task contains additional information.
!
!       isave is an integer work array of dimension 2.
!
!       dsave is a real(wp) work array of dimension 13.
!
!     Subprograms called
!
!       MINPACK-2 ... dcstep
!
!     MINPACK-1 Project. June 1983.
!     Argonne National Laboratory.
!     Jorge J. More' and David J. Thuente.
!
!     MINPACK-2 Project. October 1993.
!     Argonne National Laboratory and University of Minnesota.
!     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
!
!     **********
      real(wp), parameter :: p5     = 0.5_wp
      real(wp), parameter :: p66    = 0.66_wp
      real(wp), parameter :: xtrapl = 1.1_wp
      real(wp), parameter :: xtrapu = 4.0_wp

      logical brackt
      integer stage
      real(wp) finit , ftest , fm , fx , fxm , fy , fym ,       &
                     & ginit , gtest , gm , gx , gxm , gy , gym , stx , &
                     & sty , stmin , stmax , width , width1

!     Initialization block.

      if ( Task(1:5)=='START' ) then

!        Check the input arguments for errors.

         if ( Stp<Stpmin ) Task = 'ERROR: STP < STPMIN'
         if ( Stp>Stpmax ) Task = 'ERROR: STP > STPMAX'
         if ( g>=zero ) Task = 'ERROR: INITIAL G >= ZERO'
         if ( Ftol<zero ) Task = 'ERROR: FTOL < ZERO'
         if ( Gtol<zero ) Task = 'ERROR: GTOL < ZERO'
         if ( Xtol<zero ) Task = 'ERROR: XTOL < ZERO'
         if ( Stpmin<zero ) Task = 'ERROR: STPMIN < ZERO'
         if ( Stpmax<Stpmin ) Task = 'ERROR: STPMAX < STPMIN'

!        Exit if there are errors on input.

         if ( Task(1:5)=='ERROR' ) return

!        Initialize local variables.

         brackt = .false.
         stage = 1
         finit = f
         ginit = g
         gtest = Ftol*ginit
         width = Stpmax - Stpmin
         width1 = width/p5

!        The variables stx, fx, gx contain the values of the step,
!        function, and derivative at the best step.
!        The variables sty, fy, gy contain the value of the step,
!        function, and derivative at sty.
!        The variables stp, f, g contain the values of the step,
!        function, and derivative at stp.

         stx = zero
         fx = finit
         gx = ginit
         sty = zero
         fy = finit
         gy = ginit
         stmin = zero
         stmax = Stp + xtrapu*Stp
         Task = 'FG'

         call save_locals()
         return

      else

!        Restore local variables.

         if ( Isave(1)==1 ) then
            brackt = .true.
         else
            brackt = .false.
         endif
         stage = Isave(2)
         ginit = Dsave(1)
         gtest = Dsave(2)
         gx = Dsave(3)
         gy = Dsave(4)
         finit = Dsave(5)
         fx = Dsave(6)
         fy = Dsave(7)
         stx = Dsave(8)
         sty = Dsave(9)
         stmin = Dsave(10)
         stmax = Dsave(11)
         width = Dsave(12)
         width1 = Dsave(13)

      endif

!     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
!     algorithm enters the second stage.

      ftest = finit + Stp*gtest
      if ( stage==1 .and. f<=ftest .and. g>=zero ) stage = 2

!     Test for warnings.

      if ( brackt .and. (Stp<=stmin .or. Stp>=stmax) )              &
          & Task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
      if ( brackt .and. stmax-stmin<=Xtol*stmax )                     &
          & Task = 'WARNING: XTOL TEST SATISFIED'
      if ( Stp==Stpmax .and. f<=ftest .and. g<=gtest )            &
          & Task = 'WARNING: STP = STPMAX'
      if ( Stp==Stpmin .and. (f>ftest .or. g>=gtest) )           &
          & Task = 'WARNING: STP = STPMIN'

!     Test for convergence.

      if ( f<=ftest .and. abs(g)<=Gtol*(-ginit) )                   &
          & Task = 'CONVERGENCE'

!     Test for termination.

      if ( Task(1:4)=='WARN' .or. Task(1:4)=='CONV' ) then
         call save_locals()
         return
      end if

!     A modified function is used to predict the step during the
!     first stage if a lower function value has been obtained but
!     the decrease is not sufficient.

      if ( stage==1 .and. f<=fx .and. f>ftest ) then

!        Define the modified function and derivative values.

         fm = f - Stp*gtest
         fxm = fx - stx*gtest
         fym = fy - sty*gtest
         gm = g - gtest
         gxm = gx - gtest
         gym = gy - gtest

!        Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx,fxm,gxm,sty,fym,gym,Stp,fm,gm,brackt,stmin,    &
                   & stmax)

!        Reset the function and derivative values for f.

         fx = fxm + stx*gtest
         fy = fym + sty*gtest
         gx = gxm + gtest
         gy = gym + gtest

      else

!       Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx,fx,gx,sty,fy,gy,Stp,f,g,brackt,stmin,stmax)

      endif

!     Decide if a bisection step is needed.

      if ( brackt ) then
         if ( abs(sty-stx)>=p66*width1 ) Stp = stx + p5*(sty-stx)
         width1 = width
         width = abs(sty-stx)
      endif

!     Set the minimum and maximum steps allowed for stp.

      if ( brackt ) then
         stmin = min(stx,sty)
         stmax = max(stx,sty)
      else
         stmin = Stp + xtrapl*(Stp-stx)
         stmax = Stp + xtrapu*(Stp-stx)
      endif

!     Force the step to be within the bounds stpmax and stpmin.

      Stp = max(Stp,Stpmin)
      Stp = min(Stp,Stpmax)

!     If further progress is not possible, let stp be the best
!     point obtained during the search.

      if ( brackt .and. (Stp<=stmin .or. Stp>=stmax) .or.           &
         & (brackt .and. stmax-stmin<=Xtol*stmax) ) Stp = stx

!     Obtain another function and derivative.

      Task = 'FG'

      call save_locals()

      contains

      subroutine save_locals()

      !! Save local variables.

         if ( brackt ) then
            Isave(1) = 1
         else
            Isave(1) = 0
         endif
         Isave(2) = stage
         Dsave(1) = ginit
         Dsave(2) = gtest
         Dsave(3) = gx
         Dsave(4) = gy
         Dsave(5) = finit
         Dsave(6) = fx
         Dsave(7) = fy
         Dsave(8) = stx
         Dsave(9) = sty
         Dsave(10) = stmin
         Dsave(11) = stmax
         Dsave(12) = width
         Dsave(13) = width1

      end subroutine save_locals

      end subroutine dcsrch

      subroutine dcstep(Stx,Fx,Dx,Sty,Fy,Dy,Stp,Fp,Dp,Brackt,Stpmin,    &
                      & Stpmax)
      implicit none

      logical Brackt
      real(wp) Stx , Fx , Dx , Sty , Fy , Dy , Stp , Fp , Dp ,  &
                     & Stpmin , Stpmax
!     **********
!
!     Subroutine dcstep
!
!     This subroutine computes a safeguarded step for a search
!     procedure and updates an interval that contains a step that
!     satisfies a sufficient decrease and a curvature condition.
!
!     The parameter stx contains the step with the least function
!     value. If brackt is set to .true. then a minimizer has
!     been bracketed in an interval with endpoints stx and sty.
!     The parameter stp contains the current step.
!     The subroutine assumes that if brackt is set to .true. then
!
!           min(stx,sty) < stp < max(stx,sty),
!
!     and that the derivative at stx is negative in the direction
!     of the step.
!
!     The subroutine statement is
!
!       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
!                         stpmin,stpmax)
!
!     where
!
!       stx is a real(wp) variable.
!         On entry stx is the best step obtained so far and is an
!            endpoint of the interval that contains the minimizer.
!         On exit stx is the updated best step.
!
!       fx is a real(wp) variable.
!         On entry fx is the function at stx.
!         On exit fx is the function at stx.
!
!       dx is a real(wp) variable.
!         On entry dx is the derivative of the function at
!            stx. The derivative must be negative in the direction of
!            the step, that is, dx and stp - stx must have opposite
!            signs.
!         On exit dx is the derivative of the function at stx.
!
!       sty is a real(wp) variable.
!         On entry sty is the second endpoint of the interval that
!            contains the minimizer.
!         On exit sty is the updated endpoint of the interval that
!            contains the minimizer.
!
!       fy is a real(wp) variable.
!         On entry fy is the function at sty.
!         On exit fy is the function at sty.
!
!       dy is a real(wp) variable.
!         On entry dy is the derivative of the function at sty.
!         On exit dy is the derivative of the function at the exit sty.
!
!       stp is a real(wp) variable.
!         On entry stp is the current step. If brackt is set to .true.
!            then on input stp must be between stx and sty.
!         On exit stp is a new trial step.
!
!       fp is a real(wp) variable.
!         On entry fp is the function at stp
!         On exit fp is unchanged.
!
!       dp is a real(wp) variable.
!         On entry dp is the the derivative of the function at stp.
!         On exit dp is unchanged.
!
!       brackt is an logical variable.
!         On entry brackt specifies if a minimizer has been bracketed.
!            Initially brackt must be set to .false.
!         On exit brackt specifies if a minimizer has been bracketed.
!            When a minimizer is bracketed brackt is set to .true.
!
!       stpmin is a real(wp) variable.
!         On entry stpmin is a lower bound for the step.
!         On exit stpmin is unchanged.
!
!       stpmax is a real(wp) variable.
!         On entry stpmax is an upper bound for the step.
!         On exit stpmax is unchanged.
!
!     MINPACK-1 Project. June 1983
!     Argonne National Laboratory.
!     Jorge J. More' and David J. Thuente.
!
!     MINPACK-2 Project. October 1993.
!     Argonne National Laboratory and University of Minnesota.
!     Brett M. Averick and Jorge J. More'.
!
!     **********

      real(wp),parameter :: p66 = 0.66_wp

      real(wp) gamma , p , q , r , s , sgnd , stpc , stpf ,     &
                     & stpq , theta

      sgnd = Dp*(Dx/abs(Dx))

!     First case: A higher function value. The minimum is bracketed.
!     If the cubic step is closer to stx than the quadratic step, the
!     cubic step is taken, otherwise the average of the cubic and
!     quadratic steps is taken.

      if ( Fp>Fx ) then
         theta = three*(Fx-Fp)/(Stp-Stx) + Dx + Dp
         s = max(abs(theta),abs(Dx),abs(Dp))
         gamma = s*sqrt((theta/s)**2-(Dx/s)*(Dp/s))
         if ( Stp<Stx ) gamma = -gamma
         p = (gamma-Dx) + theta
         q = ((gamma-Dx)+gamma) + Dp
         r = p/q
         stpc = Stx + r*(Stp-Stx)
         stpq = Stx + ((Dx/((Fx-Fp)/(Stp-Stx)+Dx))/two)*(Stp-Stx)
         if ( abs(stpc-Stx)<abs(stpq-Stx) ) then
            stpf = stpc
         else
            stpf = stpc + (stpq-stpc)/two
         endif
         Brackt = .true.

!     Second case: A lower function value and derivatives of opposite
!     sign. The minimum is bracketed. If the cubic step is farther from
!     stp than the secant step, the cubic step is taken, otherwise the
!     secant step is taken.

      elseif ( sgnd<zero ) then
         theta = three*(Fx-Fp)/(Stp-Stx) + Dx + Dp
         s = max(abs(theta),abs(Dx),abs(Dp))
         gamma = s*sqrt((theta/s)**2-(Dx/s)*(Dp/s))
         if ( Stp>Stx ) gamma = -gamma
         p = (gamma-Dp) + theta
         q = ((gamma-Dp)+gamma) + Dx
         r = p/q
         stpc = Stp + r*(Stx-Stp)
         stpq = Stp + (Dp/(Dp-Dx))*(Stx-Stp)
         if ( abs(stpc-Stp)>abs(stpq-Stp) ) then
            stpf = stpc
         else
            stpf = stpq
         endif
         Brackt = .true.

!     Third case: A lower function value, derivatives of the same sign,
!     and the magnitude of the derivative decreases.

      elseif ( abs(Dp)<abs(Dx) ) then

!        The cubic step is computed only if the cubic tends to infinity
!        in the direction of the step or if the minimum of the cubic
!        is beyond stp. Otherwise the cubic step is defined to be the
!        secant step.

         theta = three*(Fx-Fp)/(Stp-Stx) + Dx + Dp
         s = max(abs(theta),abs(Dx),abs(Dp))

!        The case gamma = 0 only arises if the cubic does not tend
!        to infinity in the direction of the step.

         gamma = s*sqrt(max(zero,(theta/s)**2-(Dx/s)*(Dp/s)))
         if ( Stp>Stx ) gamma = -gamma
         p = (gamma-Dp) + theta
         q = (gamma+(Dx-Dp)) + gamma
         r = p/q
         if ( r<zero .and. gamma/=zero ) then
            stpc = Stp + r*(Stx-Stp)
         elseif ( Stp>Stx ) then
            stpc = Stpmax
         else
            stpc = Stpmin
         endif
         stpq = Stp + (Dp/(Dp-Dx))*(Stx-Stp)

         if ( Brackt ) then

!           A minimizer has been bracketed. If the cubic step is
!           closer to stp than the secant step, the cubic step is
!           taken, otherwise the secant step is taken.

            if ( abs(stpc-Stp)<abs(stpq-Stp) ) then
               stpf = stpc
            else
               stpf = stpq
            endif
            if ( Stp>Stx ) then
               stpf = min(Stp+p66*(Sty-Stp),stpf)
            else
               stpf = max(Stp+p66*(Sty-Stp),stpf)
            endif
         else

!           A minimizer has not been bracketed. If the cubic step is
!           farther from stp than the secant step, the cubic step is
!           taken, otherwise the secant step is taken.

            if ( abs(stpc-Stp)>abs(stpq-Stp) ) then
               stpf = stpc
            else
               stpf = stpq
            endif
            stpf = min(Stpmax,stpf)
            stpf = max(Stpmin,stpf)
         endif

!     Fourth case: A lower function value, derivatives of the same sign,
!     and the magnitude of the derivative does not decrease. If the
!     minimum is not bracketed, the step is either stpmin or stpmax,
!     otherwise the cubic step is taken.

      else
         if ( Brackt ) then
            theta = three*(Fp-Fy)/(Sty-Stp) + Dy + Dp
            s = max(abs(theta),abs(Dy),abs(Dp))
            gamma = s*sqrt((theta/s)**2-(Dy/s)*(Dp/s))
            if ( Stp>Sty ) gamma = -gamma
            p = (gamma-Dp) + theta
            q = ((gamma-Dp)+gamma) + Dy
            r = p/q
            stpc = Stp + r*(Sty-Stp)
            stpf = stpc
         elseif ( Stp>Stx ) then
            stpf = Stpmax
         else
            stpf = Stpmin
         endif
      endif

!     Update the interval which contains a minimizer.

      if ( Fp>Fx ) then
         Sty = Stp
         Fy = Fp
         Dy = Dp
      else
         if ( sgnd<zero ) then
            Sty = Stx
            Fy = Fx
            Dy = Dx
         endif
         Stx = Stp
         Fx = Fp
         Dx = Dp
      endif

!     Compute the new step.

      Stp = stpf

      end subroutine dcstep

end module lbfgsb_module
