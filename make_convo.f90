program covar

implicit none
integer :: i,j,n,k,emin,emax
real*8 :: En, sig, pi, spec
real*8,allocatable :: E(:), mu(:) 

n=48 !55000
emin=10000
emax=40000
pi = 4.d0*datan(1.d0)
sig = 0.02d0

allocate(E(n),mu(n))
open(15,file='Absorption-imp.dat')

do k=1,n
read(15,*) E(k), mu(k)
enddo

do i=emin,emax
En = real(i)/10000.d0
!write(6,*) En
spec = 0.d0
do j=1,n
spec = spec + mu(j)/(sig*sqrt(2.d0*pi)) * exp(-(En - E(j))**2 / (2.d0*sig**2))
!write(6,*) (-(En - E(j))**2 / (2.d0*sig**2))
enddo
write(6,*) En, spec
enddo

end 

