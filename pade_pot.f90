module precs
  implicit none
  integer,    parameter, public :: prec=16
  real(prec), parameter, public :: pi = 3.1415926535897932384626433832795028841971693993751_prec

end module precs

module moments_module
!    This module contains a set of routines to calculate moments of Green function.
!
!   TODO
!   1. Moments from physics
!
use precs

implicit none

public def_moments_w

contains

subroutine def_moments_w( func, x, moment, lm1, lm2, lm3 )
!
!   This subroutine numerically defines moments of complex function
!   using high-energy expansion M1/x+M2/(x*x)+M3/(x*x*x) and least square method.
!   
  implicit none
  
  complex(prec) :: func(:)
  real(prec)    :: x(:)
  real(prec)    :: moment(3)
  logical       :: lm1, lm2, lm3
  
  integer    :: icase
  real(prec) :: a, b, c, d, e, f
  
  icase = 0
  if( lm1 ) icase = icase + 1
  if( lm2 ) icase = icase + 2
  if( lm3 ) icase = icase + 4
  
  if( icase == 0 ) return
  
  a = sum( 1/(x*x) )
  b = sum( 1/(x*x*x*x) )
  c = sum( 1/(x*x*x*x*x*x) )
  
  d = imag( sum( func/x ) )
  e = imag( sum( func/(x*x*x) ) )
  f = real( sum( func/(x*x) ) )
  
!    write(6,*) 'in moments module'
!    write(6,*) a, b, c, d, e, f

  select case (icase)
    case(7)
      moment(1) = ( b*e - c*d ) / ( a*c - b*b )
      moment(2) = - f / b
      moment(3) = ( a*e - b*d ) / ( a*c - b*b )
    case(6)
      moment(2) = - f / b
      moment(3) = ( moment(1)*b + e ) / c
    case(5)
      moment(1) = ( b*e - c*d ) / ( a*c - b*b )
      moment(3) = ( a*e - b*d ) / ( a*c - b*b )
    case(4)
      moment(3) = ( moment(1)*b + e ) / c
    case(3)
      moment(1) = ( moment(3)*b - d ) / a
      moment(2) = - f / b
    case(2)
      moment(2) = - f / b
    case(1)
      moment(1) = ( moment(3)*b - d ) / a
    case default
      stop
  end select

end subroutine def_moments_w

end module moments_module

program analytical_continuation
!
!    Makes and analytical continuation of the complex function
!    to arbitrary complex energy mesh
!
!   Written by A. Poteryaev <Alexander.Poteryaev _at_ cpht.polytechnique.fr>
!
use precs
use moments_module

implicit none

integer :: i, ios, n_in, n_pade, n1, n2, n3, n4, n_out_im, n_out_re
integer :: index4pade(8)
real(prec) :: e1, ref, imf, beta, emin, emax, eta, m0, simpson_integral
real(prec) :: moms(4)
real(prec), allocatable :: e_pade_re(:), sf_pade_re(:)
complex(prec) :: zz
complex(prec), allocatable :: z_in(:), f_in(:), z_pade_in(:), f_pade_in(:),                   &
                              z_pade_im(:), f_pade_im(:), f_pade_re(:)
!!!!!!!
! OEP
!!!!!!!
complex(prec) :: finf, fzero

logical :: print_usage = .false., hilb_trans

character(64) :: argname, filename_in, filename_out

!============================================================================================  
!   Input file

call getarg( 1, argname )
argname = trim(adjustl(argname))

!   Check the presence of argument
  
if( len_trim(argname) == 0 )  print_usage = .true.
  
filename_in = trim(adjustl(argname))
open( 10, file=filename_in, form='formatted', status='old', iostat=ios, position='rewind' )
if( ios /= 0 )  print_usage = .true.
  
!   Print help screen

if( print_usage )then
  write(6,"(A)")
  write(6,"(A)")' Analytical_continuation:  missing input file'
  write(6,"(A)")' Usage:'
  write(6,"(A)")'   analytical_continuation filename'
  write(6,"(A)")'    filename - (required) input file with a function on the Matsubara mesh'
  write(6,"(A)")
  stop
end if  

!   Read input file

i = 0
do
  i = i + 1
  read(10,*,iostat=ios) e1
  if( ios /= 0 ) exit
end do
n_in = i - 1
ios  = 0
rewind(10)

allocate( z_in(n_in), f_in(n_in) )

do i = 1, n_in 
  read(10,*,iostat=ios) e1, ref, imf
  z_in(i) = cmplx(0,e1,prec)
  f_in(i) = cmplx(ref,imf,prec)
end do
close(10)

if( ios /= 0 .or. n_in < 1 )  stop ' Problem with the input data format'

!    Some quantitative analysis

write(6,"(/,80('='),//,' Matsubara mesh has        ',i6,' points')") n_in
beta = 2*pi / imag(z_in(2)-z_in(1))
write(6,"(' Inverse temperature, beta = ',f10.6)") beta

write(6,"(' First Matsubara point is at ',f10.6,'    function value ',2f14.7)") imag(z_in(1)),    f_in(1)
write(6,"(' Last  Matsubara point is at ',f10.6,'    function value ',2f14.7)") imag(z_in(n_in)), f_in(n_in)

n1 = minloc(real(f_in,prec),dim=1)
n2 = maxloc(real(f_in,prec),dim=1)
write(6,"(' Minimum of real part at point ',i6,' value ',f10.6,'    function value ',2f14.7,/,      &
          ' Maximum of real part at point ',i6,' value ',f10.6,'    function value ',2f14.7)")      &
           n1, imag(z_in(n1)), f_in(n1), n2, imag(z_in(n2)), f_in(n2)
           
n1 = minloc(imag(f_in),dim=1)
n2 = maxloc(imag(f_in),dim=1)
write(6,"(' Minimum of imaginary part at point ',i6,' value ',f10.6,'    function value ',2f14.7,/,      &
          ' Maximum of imaginary part at point ',i6,' value ',f10.6,'    function value ',2f14.7)")      &
           n1, imag(z_in(n1)), f_in(n1), n2, imag(z_in(n2)), f_in(n2)
           
write(6,"(/,80('='))")

!==============================================================================  
!     Begining of interactive input

write(6,"(/,' Select points at low, middle and high energies for Pade approximation')")

index4pade(1) = 1
write(6,"(/,'   Low-energy input',/,' Enter starting point (default=1)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) index4pade(1)
end if
if( index4pade(1) < 1 ) index4pade(1) = 1

index4pade(2) = 10
write(6,"(' Enter end point at low energy (default=10)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) index4pade(2)
end if

index4pade(3) = 15
write(6,"(/,'   Middle-energy input (first part)')")
write(6,"('Enter starting point (default=15)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) index4pade(3)
end if

index4pade(4) = index4pade(3)
write(argname,*) index4pade(4)
write(6,*) 'Enter end point (default='//trim(adjustl(argname))//')'
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) index4pade(4)
end if

index4pade(5) = 30
write(6,"(/,'   Middle-energy input (second part)')")
write(6,"('Enter starting point (default=30)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) index4pade(5)
end if

index4pade(6) = index4pade(5)
write(argname,*) index4pade(6)
write(6,*) 'Enter end point (default='//trim(adjustl(argname))//')'
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) index4pade(6)
end if

index4pade(7) = n_in - 1
write(argname,*) index4pade(7)
write(6,"(/,'   High-energy input')")
write(6,*)'Enter starting point (default='//trim(adjustl(argname))//')'
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) index4pade(7)
end if

index4pade(8) = n_in
write(argname,*) index4pade(8)
write(6,*) 'Enter end point (default='//trim(adjustl(argname))//')'
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) index4pade(8)
end if

!   Prepare input for Pade

n_pade = sum( index4pade(2:8:2)-index4pade(1:7:2) ) + 4
  
write(6,"(' Number of points for Pade approximation ',i5)") n_pade
write(6,"(' Intervals to be used for Pade approximation', 4(2x,i0,'..',i0) )") index4pade
write(6,"(/,80('='),/)") 

allocate( z_pade_in(n_pade), f_pade_in(n_pade) )

n1 =      index4pade(2) - index4pade(1) + 1
n2 = n1 + index4pade(4) - index4pade(3) + 1 
n3 = n2 + index4pade(6) - index4pade(5) + 1
n4 = n3 + index4pade(8) - index4pade(7) + 1 

z_pade_in(1:n1)    = z_in(index4pade(1):index4pade(2))
z_pade_in(n1+1:n2) = z_in(index4pade(3):index4pade(4))
z_pade_in(n2+1:n3) = z_in(index4pade(5):index4pade(6))
z_pade_in(n3+1:n4) = z_in(index4pade(7):index4pade(8))

f_pade_in(1:n1)    = f_in(index4pade(1):index4pade(2))
f_pade_in(n1+1:n2) = f_in(index4pade(3):index4pade(4))
f_pade_in(n2+1:n3) = f_in(index4pade(5):index4pade(6))
f_pade_in(n3+1:n4) = f_in(index4pade(7):index4pade(8))

!==============================================================================  
!     Pade part

n_out_im = 4 * n_in + 10
allocate( z_pade_im(n_out_im), f_pade_im(n_out_im) )

if ( abs(imag(z_in(1))) < 1e-6/beta ) then
  write(6,"(' Bosonic frequencies are used ')")
!  bosonic_freq = .true.

  e1 = imag(z_in(2)) / real(5,prec)

  z_pade_im(6:10)        = (/ ( cmplx(0, e1*(i-1),prec),       i = 1,5 ) /)
  z_pade_im(11:n_out_im) = (/ ( cmplx(0, 0.5*pi*i/beta, prec), i = 0,n_out_im-9 ) /) + z_in(2)
else
!  bosonic_freq = .false.

  e1 = imag(z_in(1)) / real(11,prec)

  z_pade_im(1:10)        = (/ ( cmplx(0, e1*i,prec),       i = 1,10 ) /)
  z_pade_im(11:n_out_im) = (/ ( cmplx(0, 0.5*pi*i/beta, prec), i = 0,n_out_im-9 ) /) + z_in(1)
end if
  
!   Calculate Pade coefficients

call pade_coefficients( f_pade_in, z_pade_in, n_pade )

do i = 1, n_out_im
  call pade( f_pade_im(i), z_pade_im(i), z_pade_in, f_pade_in, n_pade )
end do

filename_out = trim(adjustl(filename_in))//'_im.dat'
open( 20, file=filename_out, form='formatted', position='rewind' )

do i = 1, n_out_im
  write(20,"(1x,f20.12,8(1x,2f20.12))") imag(z_pade_im(i)), f_pade_im(i)
end do

filename_out = trim(adjustl(filename_in))//'_re.dat'
open( 30, file=filename_out, form='formatted', position='rewind' )

!    Real axis mesh

write(6,"(10x,' Real axis mesh',/)")

emin = -20._prec
write(6,"(' Enter lowest energy (default emin=-20)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) emin
end if

emax = 20._prec
write(6,"(' Enter highest energy (default emax=20)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) emax
end if

n_out_re = 4001
write(6,"(' Enter number of points (default n_out_re=4001)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) n_out_re
end if

allocate( f_pade_re(n_out_re), e_pade_re(n_out_re), sf_pade_re(n_out_re) )

eta = 0.05_prec
write(6,"(' Enter offset to imaginary plan (default eta=0.05)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(adjustl(argname))
  read(argname,*) eta
end if

!!!!!!
! OEP
!!!!!!
hilb_trans = .false.
write(6,"(' Perform Hilbert transform? (y/n)')")
read(5,"(a)") argname
if( len_trim(argname) /= 0 )then
  argname = trim(argname)

  select case ( argname(1:1) )
    case ('y','Y')
      hilb_trans = .true.
    case default
  end select
end if

!    Analytical continuation to the real axis

e1 = ( emax - emin ) / real( n_out_re - 1 ,prec)
do i = 1, n_out_re
  e_pade_re(i) = emin + e1*(i-1)
  zz = cmplx(e_pade_re(i),eta,prec)

  call pade( f_pade_re(i), zz, z_pade_in, f_pade_in, n_pade )
  write(30,"(1x,f20.12,8(1x,2f20.12))") e_pade_re(i), f_pade_re(i)
end do

!    Calculation of moments 

write(6,"(/,80('='),//,10x,' Moments of the spectral function',/)")

!   On real axis
sf_pade_re = - imag(f_pade_re) / pi
if( any( sf_pade_re < 0 ) ) write(6,*) ' There are positive elements in spectral function'

moms(1) = simpson_integral( sf_pade_re, n_out_re, e1 )

sf_pade_re = sf_pade_re * e_pade_re
moms(2) = simpson_integral( sf_pade_re, n_out_re, e1 )

sf_pade_re = sf_pade_re * e_pade_re
moms(3) = simpson_integral( sf_pade_re, n_out_re, e1 )

sf_pade_re = sf_pade_re * e_pade_re
moms(4) = simpson_integral( sf_pade_re, n_out_re, e1 )

write(6,"(' On the real energy axis',/,' M0 = ',es17.8,/,                 &
          ' M1 = ',es17.8,/,' M2 = ',es17.8,/,' M3 = ',es17.8)") moms

!   On imaginary axis
call def_moments_w( f_in(n_in-100:n_in), imag(z_in(n_in-100:n_in)),      &
                    moms(1:3), .true., .true., .true. )

write(6,"(' On the imaginary energy axis',/,' M0 = ',es17.8,/,            &
          ' M1 = ',es17.8,/,' M2 = ',es17.8)") moms(1:3)

finf = cmplx(real(f_pade_in(1)), 0.d0, prec)

write(6,"(/,' Value at infinity: ', f20.15)") real(finf)

! Find zero
i = minloc(abs(e_pade_re), dim = 1)

ref = real(f_pade_re(i+1) - f_pade_re(i-1), prec) / ( 2.d0*e1 )

write(6,"(/,' d Re F/d w( w = 0 ): ', f20.15)") real(ref)

ref = imag( z_pade_im(2) - z_pade_im(1) )
if ( ref > 1e-10 ) then
  imf = imag(f_pade_im(2) - f_pade_im(1)) / ref 
  write(6,"(/,' d Im F/d iw_n( n = 1 ): ', f20.15)") real(imf)
end if

!
! Value at (0,0)
!
zz = cmplx(0.d0, 0d0, prec)
call pade( fzero, zz, z_pade_in, f_pade_in, n_pade )

write(6,"(/,' F( w = 0 ): ', f20.15, ' + i*', f20.15)") real(fzero), imag(fzero)

!!!!!!
! OEP
!!!!!!
if( hilb_trans) then
  filename_out = trim(adjustl(filename_in))//'_hilb.dat'
  open( 40, file=filename_out, form='formatted', position='rewind' )

! The imaginary part of the value at infinity is put to zero
! Makes sense for both the GF and self-energy

  sf_pade_re = imag(f_pade_re)

  write(6,"(/, ' Performing Hilbert transform...')")
  
! ******** DEBUG **********
! In fact it is not trivial to get a good estimate of the value
! at infinity in general
  finf = 0.0
! ************************

  call hilbert_transform( sf_pade_re, n_out_re, emin, e1, eta, z_in, n_in, finf, f_in )

    do i = 1, n_in
      write(40,"(1x,f20.12,8(1x,2f20.12))") imag(z_in(i)), f_in(i), f_pade_im(i)
    enddo 

    write(6,"(' Done')")

    close(40)
  endif

end program analytical_continuation


subroutine pade_coefficients( a, z, n )
!  
!             Computation of Pade coefficients.		    
!    Inputs:							    
!io    a  -  On input is a function for approximation
!            on output contains Pade coefficients		
!i     z  -  complex points in which function a is determined
!i     n  -  size of arrays.					
! 
!r   Remarks: 						
!         (J. of Low Temp. Phys., v29, n3/4, 1977)		       
!						    
  use precs
  
  implicit none
!--->  Passed variables
  integer,       intent(in   ) :: n
  complex(prec), intent(in   ) :: z(n)
  complex(prec), intent(inout) :: a(n)
!--->  Local variables 
  integer :: i, j
  complex(prec), allocatable :: g(:,:)
     
  allocate( g(n,n) )
  
  g(1,:) = a   

  do j = 2,n
    do i = 2,j
      g(i,j) = ( g(i-1,i-1)-g(i-1,j) ) / ( z(j)-z(i-1) ) / g(i-1,j)
    end do
  end do
  
  forall( i = 1:n ) a(i) = g(i,i)
    
  deallocate( g )
    
end subroutine pade_coefficients

subroutine pade( f, z1, z, p, n )
!
!	     Calculation of analytical function	
!    in the arbitrary complex point for a given Pade coefficients 
!								       
!    Inputs:							       
!i     p  -  Pade coefficients				       
!i     z  -  set of points from which analytical continue is performed
!i     n  -  size of arrays				       
!i     z1 -  complex point					       
!    Outputs:							       
!i     f  -  value of function
!
  use precs
  
  implicit none
!---> Passed variables
  integer,       intent(in   ) :: n
  complex(prec), intent(in   ) :: z1
  complex(prec), intent(in   ) :: z(n)
  complex(prec), intent(in   ) :: p(n)
  complex(prec), intent(  out) :: f
!---> Local variables
  integer       :: i, np
  real(prec)    :: cof1, cof2, d_sum
  complex(prec) :: a1, a2, b1, b2, anew, bnew

!  d_sum = 0.d0
!  do i = 1,n-1
!    d_sum = d_sum + abs( (z1-z(i))*p(i+1) )
!  end do  

!  if( d_sum >= 1.d0 )then
!    d_sum = nint(n*d_sum/(n-1)) 
!    d_sum = log10(d_sum)
!    np  = int(d_sum) + 1
!    cof1 = 1.d-250
!    cof2 = 1.d0
!    if( np > 275 ) cof2 = 1.d-275
!    else  
!    d_sum = nint((n-1)/(n*d_sum)) 
!    d_sum = log10(d_sum)
!    np  = -int(d_sum) - 1
!    cof1 = 1.d+250
!    cof2 = 1.d0
!    if( np < -275 ) cof2 = 1.d+275
!  end if

  cof1 = 1._prec
  cof2 = 1._prec

  a1 = cmplx(0,0,prec)
  a2 = cof1*p(1)
  b1 = cmplx(cof1,0,prec)
  b2 = cmplx(cof1,0,prec)
  
  do i = 1,n-1
    anew = a2 + ( z1 - z(i) ) * p(i+1) * a1
    bnew = b2 + ( z1 - z(i) ) * p(i+1) * b1
    a1   = a2
    b1   = b2
    a2   = anew
    b2   = bnew
  end do
    
  a2 = cof2 * a2
  b2 = cof2 * b2
  f  = a2 / b2

end subroutine pade

real(prec) function simpson_integral( f, n, de )
  use precs
  
  implicit none
  
  integer,    intent(in   ) :: n              !  Number of points on uniform mesh (odd preferred)
  real(prec), intent(in   ) :: f(n)           !  Function to integrate
  real(prec), intent(in   ) :: de             !  Energy step
  
  integer :: n1
  real(prec) :: a
  
  n1 = n
  if( mod( n,2 ) /= 1 ) n1 = n - 1
  
  a = f(1) + f(n1) + 2*sum(f(3:n1-2:2)) + 4*sum(f(2:n1-1:2))
  a = a * de / 3
  
  if( n1 /= n ) a = a + ( f(n1) + f(n) ) * de / 2
  
  simpson_integral = a
  
  if( n1 == n ) return
  
  a = f(2) + f(n) + 2*sum(f(4:n-2:2)) + 4*sum(f(3:n-1:2))
  a = a * de / 3
  
  if( n1 /= n ) a = a + ( f(1) + f(2) ) * de / 2
  
  simpson_integral = ( simpson_integral + a ) / 2

end function simpson_integral

complex(prec) function z_simpson_integral( f, n, de )
  use precs
  
  implicit none
  
  integer,    intent(in   ) :: n              !  Number of points on uniform mesh (odd preferred)
  complex(prec), intent(in   ) :: f(n)           !  Function to integrate
  real(prec), intent(in   ) :: de             !  Energy step
  
  integer :: n1
  complex(prec) :: a
  
  n1 = n
  if( mod( n,2 ) /= 1 ) n1 = n - 1
  
  a = f(1) + f(n1) + 2*sum(f(3:n1-2:2)) + 4*sum(f(2:n1-1:2))
  a = a * de / 3
  
  if( n1 /= n ) a = a + ( f(n1) + f(n) ) * de / 2
  
  z_simpson_integral = a
  
  if( n1 == n ) return
  
  a = f(2) + f(n) + 2*sum(f(4:n-2:2)) + 4*sum(f(3:n-1:2))
  a = a * de / 3
  
  if( n1 /= n ) a = a + ( f(1) + f(2) ) * de / 2
  
  z_simpson_integral = ( z_simpson_integral + a ) / 2

end function z_simpson_integral

subroutine hilbert_transform( f, ne, emin_h, de, eta, iom, nom, finf, fhilb )
!
!  Performs Hilbert transform from real energies to
!  imaginary frequencies
!
  use precs

  implicit none

  real(prec), parameter :: small = 1e-6

  integer,    intent(in   ) :: ne             !  Number of points on uniform mesh (odd preferred)
  real(prec), intent(in   ) :: f(ne)          !  Function to transform (e.g. Im Sigma)
  real(prec), intent(in   ) :: de             !  Energy step
  real(prec), intent(in   ) :: eta            !  Small displacement
  real(prec), intent(in   ) :: emin_h         !  Energy interval
  integer,    intent(in   ) :: nom            !  Number of imaginary frequencies
  complex(prec), intent(in  ):: iom(nom)      !  Imaginary frequencies
  complex(prec), intent(in  ):: finf          !  Value at infinity
  complex(prec), intent(out  ):: fhilb(nom)   !  Transformed function

  complex(prec), allocatable :: fint(:)
  real(prec) :: en
  complex(prec) :: zz, zom
  integer :: i, ie, izero

  complex(prec) :: z_simpson_integral

  allocate(fint(ne))

  do i = 1,nom

    zom = iom(i)

    if ( abs(zom) > small ) then
      do ie = 1,ne
        en = emin_h + ( ie - 1 ) * de
        zz = cmplx(en,  eta, prec)

        fint(ie) = f(ie) / ( zom - zz)
      enddo

      fhilb(i) = finf - 1.d0 / pi * z_simpson_integral(fint, ne, de)

    else
! Singularity
      izero = floor(-emin_h/de + 1.5)

      if ( izero < 1 ) then
        do ie = 1,ne
          en = emin_h + ( ie - 1 ) * de
          zz = cmplx(en, eta, prec)

          fint(ie) = f(ie) / ( -zz )
        enddo

        fhilb(i) = finf - 1.d0 / pi * z_simpson_integral(fint, ne, de)
      else

! Left form singularity
        do ie = 1, izero - 1
          en = emin_h + ( ie - 1 ) * de
          zz = cmplx(en, eta, prec)

          fint(ie) = f(ie) / ( -zz )
        enddo

        fhilb(i) = finf - 1.d0 / pi * z_simpson_integral(fint, izero - 1, de)

! Right form singularity
        do ie = izero + 1, ne
          en = emin_h + ( ie - 1 ) * de
          zz = cmplx(en, eta, prec)

          fint(ie - izero) = f(ie) / ( - zz)
        enddo

        fhilb(i) = fhilb(i) - 1.d0 / pi * z_simpson_integral(fint, ne - izero, de)

! Singularity contribution
        fhilb(i) = fhilb(i) + 1.d0/pi*( f(izero + 1) - f(izero - 1) )
      end if
    end if

  enddo
  
  deallocate(fint)

end subroutine hilbert_transform
