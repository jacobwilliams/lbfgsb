
module linpack
   use blas
contains
!
!  L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License"
!  or "3-clause license")
!  Please read attached file License.txt
!
subroutine dpofa(a,lda,n,info)
integer lda,n,info
double precision a(lda,*)
!
!     dpofa factors a double precision symmetric positive definite
!     matrix.
!
!     dpofa is usually called by dpoco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the symmetric matrix to be factored.  only the
!                diagonal and upper triangle are used.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix  r  so that  a = trans(r)*r
!                where  trans(r)  is the transpose.
!                the strict lower triangle is unaltered.
!                if  info .ne. 0 , the factorization is not complete.
!
!        info    integer
!                = 0  for normal return.
!                = k  signals an error condition.  the leading minor
!                     of order  k  is not positive definite.
!
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas ddot
!     fortran sqrt
!
!     internal variables
!
double precision t
double precision s
integer j,jm1,k
!     begin block with ...exits to 40
!
!
   do 30 j = 1, n
      info = j
      s = 0.0d0
      jm1 = j - 1
      if (jm1 .lt. 1) go to 20
      do 10 k = 1, jm1
         t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
         t = t/a(k,k)
         a(k,j) = t
         s = s + t*t
10       continue
20       continue
      s = a(j,j) - s
!     ......exit
      if (s .le. 0.0d0) go to 40
      a(j,j) = sqrt(s)
30    continue
   info = 0
40 continue
return
end

!====================== The end of dpofa ===============================

subroutine dtrsl(t,ldt,n,b,job,info)
integer ldt,n,job,info
double precision t(ldt,*),b(*)
!
!
!     dtrsl solves systems of the form
!
!                   t * x = b
!     or
!                   trans(t) * x = b
!
!     where t is a triangular matrix of order n. here trans(t)
!     denotes the transpose of the matrix t.
!
!     on entry
!
!         t         double precision(ldt,n)
!                   t contains the matrix of the system. the zero
!                   elements of the matrix are not referenced, and
!                   the corresponding elements of the array can be
!                   used to store other information.
!
!         ldt       integer
!                   ldt is the leading dimension of the array t.
!
!         n         integer
!                   n is the order of the system.
!
!         b         double precision(n).
!                   b contains the right hand side of the system.
!
!         job       integer
!                   job specifies what kind of system is to be solved.
!                   if job is
!
!                        00   solve t*x=b, t lower triangular,
!                        01   solve t*x=b, t upper triangular,
!                        10   solve trans(t)*x=b, t lower triangular,
!                        11   solve trans(t)*x=b, t upper triangular.
!
!     on return
!
!         b         b contains the solution, if info .eq. 0.
!                   otherwise b is unaltered.
!
!         info      integer
!                   info contains zero if the system is nonsingular.
!                   otherwise info contains the index of
!                   the first zero diagonal element of t.
!
!     linpack. this version dated 08/14/78 .
!     g. w. stewart, university of maryland, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!     fortran mod
!
!     internal variables
!
double precision temp
integer case,j,jj
!
!     begin block permitting ...exits to 150
!
!        check for zero diagonal elements.
!
   do 10 info = 1, n
!     ......exit
      if (t(info,info) .eq. 0.0d0) go to 150
10    continue
   info = 0
!
!        determine the task and go to it.
!
   case = 1
   if (mod(job,10) .ne. 0) case = 2
   if (mod(job,100)/10 .ne. 0) case = case + 2
   go to (20,50,80,110), case
!
!        solve t*x=b for t lower triangular
!
20    continue
      b(1) = b(1)/t(1,1)
      if (n .lt. 2) go to 40
      do 30 j = 2, n
         temp = -b(j-1)
         call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
         b(j) = b(j)/t(j,j)
30       continue
40       continue
   go to 140
!
!        solve t*x=b for t upper triangular.
!
50    continue
      b(n) = b(n)/t(n,n)
      if (n .lt. 2) go to 70
      do 60 jj = 2, n
         j = n - jj + 1
         temp = -b(j+1)
         call daxpy(j,temp,t(1,j+1),1,b(1),1)
         b(j) = b(j)/t(j,j)
60       continue
70       continue
   go to 140
!
!        solve trans(t)*x=b for t lower triangular.
!
80    continue
      b(n) = b(n)/t(n,n)
      if (n .lt. 2) go to 100
      do 90 jj = 2, n
         j = n - jj + 1
         b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
         b(j) = b(j)/t(j,j)
90       continue
100       continue
   go to 140
!
!        solve trans(t)*x=b for t upper triangular.
!
110    continue
      b(1) = b(1)/t(1,1)
      if (n .lt. 2) go to 130
      do 120 j = 2, n
         b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
         b(j) = b(j)/t(j,j)
120       continue
130       continue
140    continue
150 continue
return
end

!====================== The end of dtrsl ===============================


end module linpack