include 'specfun.f90'

program ModelQD_one

use Constants_au
use Variables_au
use Integrals
use Normal
use Make_Ham

implicit none

call getVariables

allocate(Eg(1000))
allocate(TDM(1000))
allocate(minEeA(10),minEhA(10),minEeB(10),minEhB(10),source=0._dp)
allocate(pow(0:ntime+1))
allocate(powtemp(0:ntime+1))
allocate(pow_gaus(0:ntime+1))
allocate(pulses(0:ntime+1))
allocate(pow_pol(npol,0:ntime+1))
allocate(pow_pol_gaus(npol,0:ntime+1))

open(newunit=Tmat_0_f      ,file='TransMat.dat')          
open(newunit=Tmat_ei_f     ,file='TransMat_ei.dat')       
open(newunit=H_0_f         ,file='Ham0.dat')          
open(newunit=Tmat_y_f      ,file='TransMat_ei_y.dat')     
open(newunit=Tmat_z_f      ,file='TransMat_ei_z.dat')     
open(newunit=Tmat_x_f      ,file='TransMat_ei_x.dat')     
open(newunit=H_dir_f       ,file='Ham_dir.dat')           
open(newunit=H_ex_f        ,file='Ham_ex.dat')            
open(newunit=H_JK_f        ,file='Ham_JK.dat')            
open(newunit=H_ei_f        ,file='Ham_ei.dat')            
open(newunit=Etr_0_f       ,file='Etransitions-he_0.dat') 
open(newunit=Etr_ei_f      ,file='Etransitions-he_ei.dat')
open(newunit=TDip_ei_f     ,file='TransDip_ei.dat')       
open(newunit=Abs_imp_f     ,file='Absorption-imp.dat')    
open(newunit=Liou_f        ,file='Liou.dat')
open(newunit=TransAbs      ,file='TransAbs.dat')
open(newunit=TransAbs_NR   ,file='TransAbs_NR.dat')
open(newunit=TransAbs_R    ,file='TransAbs_R.dat')
open(newunit=TransAbs_P    ,file='TransAbs_P.dat')
open(newunit=DipSpec       ,file='DipSpec.dat')
open(newunit=P_Match_f     ,file='Phase_Match.dat')
open(newunit=DipSpec_NR_f  ,file='DipSpec-NR.dat')
open(newunit=DipSpec_R_f   ,file='DipSpec-R.dat')
!open(newunit=TransAbs_7_f  ,file='TransAbs_7.dat')
!open(newunit=TransAbs_17_f ,file='TransAbs_17.dat') 
!open(newunit=TransAbs_33_f ,file='TransAbs_33.dat')
open(newunit=TransAbs_39_f ,file='TransAbs_39.dat')
open(newunit=TransAbs_41_f ,file='TransAbs_41.dat')
!open(newunit=TransAbs_43_f ,file='TransAbs_43.dat')
!open(newunit=TransAbs_44_f ,file='TransAbs_44.dat')
!open(newunit=DipSpec_7_f   ,file='DipSpec_7.dat')
!open(newunit=DipSpec_17_f  ,file='DipSpec_17.dat') 
!open(newunit=DipSpec_33_f  ,file='DipSpec_33.dat')
open(newunit=DipSpec_39_f  ,file='DipSpec_39.dat')
open(newunit=DipSpec_41_f  ,file='DipSpec_41.dat')
!open(newunit=DipSpec_43_f  ,file='DipSpec_43.dat')
!open(newunit=DipSpec_44_f  ,file='DipSpec_44.dat')
open(newunit=label_0       ,file='label_0.dat')

if ( inbox .eq. 'y' ) then
allocate(l1(npol),l2(npol),l3(npol),source=0._dp)
allocate(eTDM(3,3,3),source=0._dp)
allocate(rotmat(3,3))
allocate(Dcenter(3))
call get_phases
call make_eTDM
endif

matrices = (/ Tmat_0_f,Tmat_ei_f,Tmat_x_f,Tmat_y_f,Tmat_z_f,H_0_f,H_dir_f,H_ex_f,H_JK_f,H_ei_f /)

do n=1,nsys

if ( ( n .le. nQDA+nQDB ) .and. ( model .eq. "SB" ) ) then

if ( n .le. nQDA ) then
aRA      = r8_NORMAL_AB(aA,dispQD*aA)
call QDA
elseif ( ( n .ge. nQDA+1 ) .and. ( n .le. nQDA+nQDB ) ) then
aRA      = r8_NORMAL_AB(aB,dispQD*aB)
call QDA
endif

if ( inbox .eq. 'y' ) then
Dcenter(:) = vectorin((2._dp*t_au*2._dp*pi*cl/min(omega01,omega02,omega03))**3._dp)
endif

nstates=3 ; nstates2=nstates**2
!print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_singl
include 'Core.f90'

elseif ( ( n .le. nQDA+nQDB ) .and. ( model .eq. "FS" ) ) then

if ( n .le. nQDA ) then                                                                               
aRA      = r8_NORMAL_AB(aA,dispQD*aA)                                                               
call QDA
elseif ( ( n .ge. nQDA+1 ) .and. ( n .le. nQDA+nQDB ) ) then                                          
aRA      = r8_NORMAL_AB(aB,dispQD*aB)                                                               
call QDA
endif

if ( inbox .eq. 'y' ) then
Dcenter(:) = vectorin((2._dp*t_au*2._dp*pi*cl/min(omega01,omega02,omega03))**3._dp)
endif

nstates=55 ; nstates2=nstates**2
print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_FS_FO
include 'Core.f90'

elseif ( (n .gt. nQDA+nQDB) .and. ( model .eq. "FO" ) ) then

if ( ( n .ge. nQDA+nQDB+1 ) .and. ( n .le. nQDA+nQDB+nhomoA ) ) then                                
aRA      = r8_NORMAL_AB(aA,dispQD*aA)                                                               
aRB      = r8_NORMAL_AB(aA,dispQD*aA)                                                               
call QDA
call QDB
elseif ( ( n .ge. nQDA+nQDB+nhomoA+1 ) .and. ( n .le. nQDA+nQDB+nhomoA+nhomoB ) ) then  
aRA      = r8_NORMAL_AB(aB,dispQD*aB)                                                               
aRB      = r8_NORMAL_AB(aB,dispQD*aB)                                                               
call QDA
call QDB
elseif ( ( n .ge. nQDA+nQDB+nhomoA+nhomoB+1 ) .and. ( n .le. nQDA+nQDB+nhomoA+nhomoB+nhetero ) ) then 
aRA      = r8_NORMAL_AB(aA,dispQD*aA)                                                               
aRB      = r8_NORMAL_AB(aB,dispQD*aB)                                                               
call QDA
call QDB
endif

if ( inbox .eq. 'y' ) then
Dcenter(:) = vectorin((2._dp*t_au*2._dp*pi*cl/min(omega01,omega02,omega03))**3._dp)
endif

nstates=5 ; nstates2=nstates**2
print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_he_FO
include 'Core.f90'

elseif ( (n .gt. nQDA+nQDB) .and. ( model .eq. "FS" ) ) then

if ( ( n .ge. nQDA+nQDB+1 ) .and. ( n .le. nQDA+nQDB+nhomoA ) ) then
aRA      = r8_NORMAL_AB(aA,dispQD*aA)
aRB      = r8_NORMAL_AB(aA,dispQD*aA)
call QDA
call QDB
elseif ( ( n .ge. nQDA+nQDB+nhomoA+1 ) .and. ( n .le. nQDA+nQDB+nhomoA+nhomoB ) ) then
aRA      = r8_NORMAL_AB(aB,dispQD*aB)
aRB      = r8_NORMAL_AB(aB,dispQD*aB)
call QDA
call QDB
elseif ( ( n .ge. nQDA+nQDB+nhomoA+nhomoB+1 ) .and. ( n .le. nQDA+nQDB+nhomoA+nhomoB+nhetero ) ) then
aRA      = r8_NORMAL_AB(aA,dispQD*aA)
aRB      = r8_NORMAL_AB(aB,dispQD*aB)
call QDA
call QDB
endif

if ( inbox .eq. 'y' ) then
Dcenter(:) = vectorin((2._dp*t_au*2._dp*pi*cl/min(omega01,omega02,omega03))**3._dp)
endif

nstates=109 ; nstates2=nstates**2
print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_FS_FO
include 'Core.f90'

elseif ( (n .gt. nQDA+nQDB) .and. ( model .eq. "SB" ) ) then

if ( ( n .ge. nQDA+nQDB+1 ) .and. ( n .le. nQDA+nQDB+nhomoA ) ) then
aRA      = r8_NORMAL_AB(aA,dispQD*aA)
aRB      = r8_NORMAL_AB(aA,dispQD*aA)
call QDA
call QDB
elseif ( ( n .ge. nQDA+nQDB+nhomoA+1 ) .and. ( n .le. nQDA+nQDB+nhomoA+nhomoB ) ) then
aRA      = r8_NORMAL_AB(aB,dispQD*aB)
aRB      = r8_NORMAL_AB(aB,dispQD*aB)
call QDA
call QDB
elseif ( ( n .ge. nQDA+nQDB+nhomoA+nhomoB+1 ) .and. ( n .le. nQDA+nQDB+nhomoA+nhomoB+nhetero ) ) then
aRA      = r8_NORMAL_AB(aA,dispQD*aA)
aRB      = r8_NORMAL_AB(aB,dispQD*aB)
call QDA
call QDB
endif

if ( inbox .eq. 'y' ) then
Dcenter(:) = vectorin((2._dp*t_au*2._dp*pi*cl/min(omega01,omega02,omega03))**3._dp)
endif

nstates=9 ; nstates2=nstates**2
!print*, "Computing system number:    ", n, "which possesses", nstates, "states"
include 'allocate_core.f90'
call make_Ham_he
include 'Core.f90'

endif

enddo !end loop number of systems

!deallocate(minEeA,minEhA,minEeB,minEhB)
close(Tmat_0_f);close(Tmat_ei_f);close(H_0_f);close(Tmat_y_f);close(Tmat_z_f);close(Tmat_x_f);close(H_dir_f)
close(H_ex_f);close(H_JK_f);close(H_ei_f);close(Etr_0_f);close(TDip_ei_f);close(Liou_f);close(Etr_ei_f)

!if ( doAbs .eq. "y" ) then
!call Convolution
!endif

!Write DipSpec sum of all systems + multiplied by a gaussian 
timestep  = timestep  * t_au
totaltime = totaltime * t_au

do t=0,ntime

time = t*timestep

pow_gaus(t)=exp(-1._dp*((time-totaltime/2._dp)*timestep)**2._dp/(2._dp*(totaltime*timestep/15._dp)**2._dp))*(pow(t)/nsys)

if ( inbox .eq. 'y' ) then
pow_pol_gaus(39,t)=&
   exp(-1._dp*((time-totaltime/2._dp)*timestep)**2._dp/(2._dp*(totaltime*timestep/15._dp)**2._dp))*(pow_pol(39,t)/nsys)
pow_pol_gaus(41,t)=&
   exp(-1._dp*((time-totaltime/2._dp)*timestep)**2._dp/(2._dp*(totaltime*timestep/15._dp)**2._dp))*(pow_pol(41,t)/nsys)
!if ( nofiles .eq. 'n' ) then
write(DipSpec_39_f,*) time, dreal(pow_pol(39,t)), dimag(pow_pol(39,t))
write(DipSpec_41_f,*) time, dreal(pow_pol(41,t)), dimag(pow_pol(41,t))
!endif
endif

write(DipSpec,*) time, pow(t), pulses(t)

enddo 

if ( inbox .eq. 'y' ) then

do pol=1,npol
integPol     = dcmplx(0._dp,0._dp)
integPolconv = dcmplx(0._dp,0._dp)
do t=0,ntime
time = t*timestep
integPol = integPol + timestep*abs(pow_pol(pol,t) + pow_pol(pol,t+1))/2._dp
enddo
write(P_Match_f,*) pol, l1(pol), l2(pol), l3(pol), integPol
enddo

endif

FTscale = h/(elec*(2._dp**FTpow)*timestep)
!gives the scale of FT vectors
t=0
do while ( t*FTscale .le. 4._dp )
t=t+1
enddo
nFT=t

if ( doFT .eq. 'y' ) then

allocate(wft(0:nFT+1))
allocate(wft_pol(npol,0:nFT+1))
allocate(wftp(0:nFT+1))
allocate(wftf(0:nFT+1))
allocate(xpow_gaus(0:nint(2._dp**FTpow)))
allocate(xpulse(0:nint(2._dp**FTpow)))

do t=0,ntime
xpow_gaus(t) = dcmplx(pow_gaus(t)  ,0._dp)
xpulse(t)    = dcmplx(pulses(t),0._dp)
enddo

do t=ntime+1,nint(2.d0**FTpow)
xpow_gaus(t) = dcmplx(0._dp,0._dp)
xpulse(t)    = dcmplx(0._dp,0._dp)
enddo

call fft(xpulse)
call fft(xpow_gaus)

t=0
do while ( t*FTscale .le. 4._dp ) 
wftf(t)= -2._dp * dimag(sqrt(dreal(xpow_gaus(t))**2+dimag(xpow_gaus(t))**2) * dconjg(xpulse(t)))
!if ( nofiles .eq. 'n' ) then
write(TransAbs,*) t*FTscale, dreal(wftf(t)), dreal(xpow_gaus(t)), dimag(xpow_gaus(t)), abs(xpulse(t))
!endif
t = t + 1 
enddo

!endif

if ( inbox .eq. 'y' ) then

allocate(xpow_pol(44,0:nint(2._dp**FTpow)))
allocate(wftf_pol(44,0:nFT+1))

do t=0,ntime
xpow_pol(39,t)  = pow_pol_gaus(39,t)
xpow_pol(41,t)  = pow_pol_gaus(41,t)
enddo
do t=ntime+1,nint(2.d0**FTpow)
xpow_pol(39,t)  = dcmplx(0._dp,0._dp)
xpow_pol(41,t)  = dcmplx(0._dp,0._dp)
enddo

!do pol=1,npol
call fft(xpow_pol(39,:))
call fft(xpow_pol(41,:))
!enddo

t=0
do while ( t*FTscale .le. 4.d0 )
write(TransAbs_39_f,*) t*FTscale, dreal(xpow_pol(39,t)), dimag(xpow_pol(39,t)), abs(xpow_pol(39,t))
write(TransAbs_41_f,*) t*FTscale, dreal(xpow_pol(41,t)), dimag(xpow_pol(41,t)), abs(xpow_pol(41,t))
t = t + 1
enddo

deallocate(xpow_pol,wftf_pol)

endif

deallocate(pow,pow_gaus,xpow_gaus,xpulse,pulses,wft,wftp)

endif

!if ( doCovar .eq. 'y' ) then
!
!allocate(Scov(0:nFT+1,0:nFT+1))
!
!open(61,file='Allcov.dat')
!
!t=0
!do while ( t*FTscale .le. 4.d0 )
!t2=0
!do while ( t2*FTscale .le. 4.d0 )
!
!do k=1,nsys
!Scov(t,t2) = Scov(t,t2) + sum(dreal(wftf_s(k,t))*dreal(wftf_s(:,t2)))
!enddo
!
!t2 = t2 + 5
!enddo
!t = t + 5
!enddo
!
!t=0
!do while ( t*FTscale .le. 4.d0 )
!t2=0
!do while ( t2*FTscale .le. 4.d0 )
!write(61,'(2f10.6,ES15.6E3)') t*h/(elec*5.24288d-12), t2*h/(elec*5.24288d-12), Scov(t,t2)
!t2 = t2 + 5
!enddo
!write(61,*)
!t = t + 5
!enddo
!
!endif

!call system('sleep 1')
!call system('find . -size 0 -delete')

!!$OMP END DO

!$OMP END PARALLEL DO

end
