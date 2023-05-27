
module blas
  contains
!
!  L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License"
!  or "3-clause license")
!  Please read attached file License.txt
!

double precision function dnrm2(n,x,incx)
integer n,incx
double precision x(n)
!     **********
!
!     Function dnrm2
!
!     Given a vector x of length n, this function calculates the
!     Euclidean norm of x with stride incx.
!
!     The function statement is
!
!       double precision function dnrm2(n,x,incx)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
!       incx is a positive integer variable that specifies the
!         stride of the vector.
!
!     Subprograms called
!
!       FORTRAN-supplied ... abs, max, sqrt
!
!     MINPACK-2 Project. February 1991.
!     Argonne National Laboratory.
!     Brett M. Averick.
!
!     **********
integer i
double precision scale

dnrm2 = 0.0d0
scale = 0.0d0

do 10 i = 1, n, incx
   scale = max(scale, abs(x(i)))
10 continue

if (scale .eq. 0.0d0) return

do 20 i = 1, n, incx
   dnrm2 = dnrm2 + (x(i)/scale)**2
20 continue

dnrm2 = scale*sqrt(dnrm2)


return

end

!====================== The end of dnrm2 ===============================

subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
double precision dx(*),dy(*),da
integer i,incx,incy,ix,iy,m,mp1,n
!
if(n.le.0)return
if (da .eq. 0.0d0) return
if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
ix = 1
iy = 1
if(incx.lt.0)ix = (-n+1)*incx + 1
if(incy.lt.0)iy = (-n+1)*incy + 1
do 10 i = 1,n
  dy(iy) = dy(iy) + da*dx(ix)
  ix = ix + incx
  iy = iy + incy
10 continue
return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,4)
if( m .eq. 0 ) go to 40
do 30 i = 1,m
  dy(i) = dy(i) + da*dx(i)
30 continue
if( n .lt. 4 ) return
40 mp1 = m + 1
do 50 i = mp1,n,4
  dy(i) = dy(i) + da*dx(i)
  dy(i + 1) = dy(i + 1) + da*dx(i + 1)
  dy(i + 2) = dy(i + 2) + da*dx(i + 2)
  dy(i + 3) = dy(i + 3) + da*dx(i + 3)
50 continue
return
end

!====================== The end of daxpy ===============================

subroutine dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
double precision dx(*),dy(*)
integer i,incx,incy,ix,iy,m,mp1,n
!
if(n.le.0)return
if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
ix = 1
iy = 1
if(incx.lt.0)ix = (-n+1)*incx + 1
if(incy.lt.0)iy = (-n+1)*incy + 1
do 10 i = 1,n
  dy(iy) = dx(ix)
  ix = ix + incx
  iy = iy + incy
10 continue
return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,7)
if( m .eq. 0 ) go to 40
do 30 i = 1,m
  dy(i) = dx(i)
30 continue
if( n .lt. 7 ) return
40 mp1 = m + 1
do 50 i = mp1,n,7
  dy(i) = dx(i)
  dy(i + 1) = dx(i + 1)
  dy(i + 2) = dx(i + 2)
  dy(i + 3) = dx(i + 3)
  dy(i + 4) = dx(i + 4)
  dy(i + 5) = dx(i + 5)
  dy(i + 6) = dx(i + 6)
50 continue
return
end

!====================== The end of dcopy ===============================

double precision function ddot(n,dx,incx,dy,incy)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
double precision dx(*),dy(*),dtemp
integer i,incx,incy,ix,iy,m,mp1,n
!
ddot = 0.0d0
dtemp = 0.0d0
if(n.le.0)return
if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
ix = 1
iy = 1
if(incx.lt.0)ix = (-n+1)*incx + 1
if(incy.lt.0)iy = (-n+1)*incy + 1
do 10 i = 1,n
  dtemp = dtemp + dx(ix)*dy(iy)
  ix = ix + incx
  iy = iy + incy
10 continue
ddot = dtemp
return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,5)
if( m .eq. 0 ) go to 40
do 30 i = 1,m
  dtemp = dtemp + dx(i)*dy(i)
30 continue
if( n .lt. 5 ) go to 60
40 mp1 = m + 1
do 50 i = mp1,n,5
  dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
50 continue
60 ddot = dtemp
return
end

!====================== The end of ddot ================================

subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!
double precision da,dx(*)
integer i,incx,m,mp1,n,nincx
!
if( n.le.0 .or. incx.le.0 )return
if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
nincx = n*incx
do 10 i = 1,nincx,incx
  dx(i) = da*dx(i)
10 continue
return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,5)
if( m .eq. 0 ) go to 40
do 30 i = 1,m
  dx(i) = da*dx(i)
30 continue
if( n .lt. 5 ) return
40 mp1 = m + 1
do 50 i = mp1,n,5
  dx(i) = da*dx(i)
  dx(i + 1) = da*dx(i + 1)
  dx(i + 2) = da*dx(i + 2)
  dx(i + 3) = da*dx(i + 3)
  dx(i + 4) = da*dx(i + 4)
50 continue
return
end

!====================== The end of dscal ===============================

end module blas