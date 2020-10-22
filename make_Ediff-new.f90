program TDM_E

implicit none
integer :: i,j,n,nstates,nsys
real*8 :: dummy
real*8,allocatable :: Ene(:,:), TDM(:,:,:), sig(:,:), moy(:,:), moyTDM(:,:)

nsys = 2000
nstates= 49

allocate(Ene(nstates,nsys),TDM(nstates,nstates,nsys),sig(nstates,nstates),moy(nstates,nstates),&
moyTDM(nstates,nstates),source=0.d0)
open(10,file='Etransitions-he_ei-Voff.dat')
open(11,file='TransMat_ei-Voff.dat')
open(12,file='E-TDM-all-32-06-Voff.dat')
open(13,file='moy-sig-32-06-Voff.dat')

do n=1,nsys
read(10,*) dummy, (Ene(i,n), i=1,nstates)
enddo

do n=1,nsys
read(11,*) 
do j=1,nstates
read(11,*) (TDM(j,i,n), i=1,nstates)
enddo
read(11,*) 
enddo

!write(6,*) TDM(1,1,1), TDM(2,1,1)

do n=1,nsys
do i=1,nstates
do j=i,nstates
moy(i,j) = moy(i,j) + abs(Ene(i,n) - Ene(j,n))
!if ( i .eq. j ) then 
!moyTDM(1,j) = moyTDM(1,j) + TDM(1,i,n)**2
!write(6,*) moyTDM(i,j), TDM(i,1,n)**2 , i, j 
!elseif ( i .ne. 1 ) then
moyTDM(i,j) = moyTDM(i,j) + abs(TDM(1,i,n))
!moyTDM(i,j) = moyTDM(i,j) + abs(TDM(1,i,n))*abs(TDM(1,j,n))
!endif
!if ( ( abs(TDM(i,1,n))*abs(TDM(j,1,n) ) .ne. 0.0d0 ) .and. ( i .le. 25 ) .and. ( j .le. 25 ) ) then 
if ( ( abs(TDM(i,1,n))*abs(TDM(j,1,n) ) .ne. 0.0d0 ) .and. ( i .le. 25.) .and. ( j .le. 25 ) ) then 
write(12,*) abs(Ene(i,n)), abs(Ene(j,n)), abs(TDM(i,1,n))*abs(TDM(j,1,n)), i-1,j-1
endif
enddo
enddo
enddo

!do i=1,nstates
!write(6,*) i, moyTDM(i,1)
!enddo

do i=1,nstates
do j=i,nstates
do n=1,nsys
sig(i,j) = sig(i,j) + (abs(Ene(i,n) - Ene(j,n)) - (moy(i,j)/nsys))**2
enddo
enddo
enddo


do i=1,nstates
do j=i+1,nstates
if ( ( i .eq. 1 ) .and. (  moy(i,j)/nsys .le. 8.d0 ) ) then !.and. ( i .le. 25.) .and. ( j .le. 25 ) ) then
write(13,*) i-1,j-1, moy(i,j)/nsys, sqrt(sig(i,j)/nsys), (moyTDM(j,j)/nsys)**2
elseif ( ( i .ne. 1 ) .and. (  moy(i,j)/nsys .le. 8.d0 ) ) then !.and. ( i .le. 25.) .and. ( j .le. 25 ) ) then
!write(13,*) i-1,j-1, moy(i,j)/nsys, sqrt(sig(i,j)/nsys), moyTDM(i,j)/nsys
write(13,*) i-1,j-1, moy(i,j)/nsys, sqrt(sig(i,j)/nsys), (moyTDM(i,i)/nsys)*(moyTDM(j,j)/nsys)
endif
enddo
enddo


end 
