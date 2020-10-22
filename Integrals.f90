module Integrals

      use Constants_au
      use Variables_au
      use Vectors
      use omp_lib

contains

!TDM homodimer h1e from fit
real(dp) function TransDip_Fit_h1e_ho(a,b)

      implicit none
      real(dp) :: d1,d2,d3,d4,a,b

d1 = -1.53548e-07_dp
d2 = 52.6179_dp
d3 = 4.36573_dp
d4 = 9.86824_dp

TransDip_Fit_h1e_ho =  ( d1 + d2 / ((a*1e9_dp)**d3) * exp(-1.0_dp*d4*b*1e9_dp) )*1.0e-33_dp

end function TransDip_Fit_h1e_ho

!TDM homodimer h2e from fit
real(dp) function TransDip_Fit_h2e_ho(a,b)

      implicit none
      real(dp) :: d1,d2,d3,d4,a,b

d1   = 2.20868e-06_dp      
d2   = 39.8905_dp          
d3   = 9.61212_dp          
d4   = 8.88849_dp

!link=1.0nm
!d1              = 1.09021e-05    
!d2              = 0.0028037      
!d3              = 1.82566        
!d4              = 2.70059

TransDip_Fit_h2e_ho = ( d1 + d2 / ((a*1e9_dp)**d3) * exp(-1.0_dp*d4*b*1e9_dp) )*1.0e-33_dp

end function TransDip_Fit_h2e_ho

!TDM homodimer h1e from fit
real(dp) function TransDip_Fit_h1e_he(a,b)
      implicit none
      real(dp) :: d1,d2,d3,d4,a,b
if ( idlink .eq. 20 ) then
d1              = 0.0087115_dp
d2              = 0.418797_dp  
d3              = 1.03237_dp   
d4              = 1.23661_dp   
else if ( idlink .eq. 55 ) then
d1              = 0.00119598_dp 
d2              = 0.0314746_dp  
d3              = 1.0704_dp     
d4              = 2.18365_dp    
endif
TransDip_Fit_h1e_he =  d1 + d2 / ((a*1e9_dp)**d3*(b*1e9_dp)**d4)
end function TransDip_Fit_h1e_he

!TDM homodimer h1e from fit
real(dp) function TransDip_Fit_h2e_he(a,b)
      implicit none
      real(dp) :: d1,d2,d3,d4,a,b
if ( idlink .eq. 20 ) then
d1              = 0.00536146_dp  
d2              = 0.834794_dp    
d3              = 1.03382_dp     
d4              = 1.11718_dp     
else if ( idlink .eq. 55 ) then
d1              = 0.0163095_dp  
d2              = 0.752361_dp   
d3              = 1.17482_dp    
d4              = 1.87855_dp    
endif
TransDip_Fit_h2e_he =  d1 + d2 / ((a*1e9_dp)**d3*(b*1e9_dp)**d4)
end function TransDip_Fit_h2e_he

!TDM homodimer h1e from fit
real(dp) function TransDip_Fit_h1h1_he(a,b)
      implicit none
      real(dp) :: d1,d2,d3,d4,a,b
if ( idlink .eq. 20 ) then
d1              = 0.000465077_dp  
d2              = 0.0461328_dp    
d3              = 1.11291_dp      
d4              = 1.11376_dp      
else if ( idlink .eq. 55 ) then
d1              = 1.27936e-05_dp   
d2              = 0.00108558_dp  
d3              = 1.15109_dp     
d4              = 1.15094_dp     
endif
TransDip_Fit_h1h1_he =  d1 + d2 / ((a*1e9_dp)**d3*(b*1e9_dp)**d4)
end function TransDip_Fit_h1h1_he

!TDM homodimer h1e from fit
real(dp) function TransDip_Fit_h2h2_he(a,b)
      implicit none
      real(dp) :: d1,d2,d3,d4,a,b
if ( idlink .eq. 20 ) then
d1                  = 0.00359745_dp 
d2                  = 0.242846_dp   
d3                  = 1.27172_dp    
d4                  = 1.26927_dp    
else if ( idlink .eq. 55 ) then
d1                   = 0.000166763_dp 
d2                   = 0.00893376_dp  
d3                   = 1.56922_dp     
d4                   = 1.57073_dp     
endif
TransDip_Fit_h2h2_he =  d1 + d2 / ((a*1e9_dp)**d3*(b*1e9_dp)**d4)
end function TransDip_Fit_h2h2_he

!TDM homodimer h1e from fit
real(dp) function TransDip_Fit_ee_he(a,b)
      implicit none
      real(dp) :: d1,d2,d3,d4,a,b
if ( idlink .eq. 20 ) then
d1                = 0.0119711_dp 
d2                = 2.18843_dp   
d3                = 0.904435_dp  
d4                = 0.906561_dp  
else if ( idlink .eq. 55 ) then
d1               = 0.00336278_dp 
d2               = 0.278173_dp   
d3               = 1.30876_dp    
d4               = 1.30738_dp    
endif
TransDip_Fit_ee_he =  d1 + d2 / ((a*1e9_dp)**d3*(b*1e9_dp)**d4)
end function TransDip_Fit_ee_he

!TDM homodimer h1e from fit
real(dp) function TransDip_Fit_h1h2_he(a,b)
      implicit none
      real(dp) :: d1,d2,d3,d4,a,b
if ( idlink .eq. 20 ) then
d1              = 0.00184203_dp   
d2              = 0.107537_dp     
d3              = 1.12685_dp      
d4              = 1.31525_dp      
else if ( idlink .eq. 55 ) then
d1                = 8.62101e-05_dp  
d2                = 0.0032382_dp    
d3                = 1.19545_dp      
d4                = 1.69668_dp      
endif
TransDip_Fit_h1h2_he =  d1 + d2 / ((a*1e9_dp)**d3*(b*1e9_dp)**d4)
end function TransDip_Fit_h1h2_he

!Transition dipole inside dot volume, spherical symmetry
real(dp) function TransDip_Num(A1,B1,kin1,kout1,A2,B2,kin2,kout2,r)

      implicit none
      integer :: i, m
      real(dp) :: x, fin, fout, interval
      real(dp), intent(in) :: A1,A2,B1,B2,kin1,kin2,kout1,kout2,r

      i=0
      m=1000

      interval= r/m

      do i=0,m-1
      x=0.5*interval+i*interval
      fin = fin + elec*x*A1*sin(kin1*x)*A2*sin(kin2*x)
      enddo

      do i=m,2*m
      x=0.5*interval+i*interval
      fout = fout + elec*x*B1*exp(-1.0*kout1*x)*B2*exp(-1.0*kout2*x)
      enddo

      TransDip_Num=abs(fin+fout)*interval

end function TransDip_Num

!Transition dipole one dot, analytic, spherical symmetry
real(dp) function TransDip_Ana(A1,A2,B1,B2,kin1,kin2,kout1,kout2,r)

      implicit none
      real(dp) :: A1,A2,B1,B2,kin1,kin2,kout1,kout2,r,rmin,rmax

      rmin=0.0_dp
      rmax=2._dp*r

      TransDip_Ana=(-1.0_dp*elec*(((A1*A2*0.5_dp*((r*sin(r*(kin1-kin2))/(kin1-kin2))-(r*sin(r*(kin1+kin2))/(kin1+kin2))+&
                       (cos(r*(kin1-kin2))/(kin1-kin2)**2)-(cos(r*(kin1+kin2))/(kin1+kin2)**2))) - &
                       A1*A2*0.5_dp*((cos(rmin*(kin1-kin2))/(kin1-kin2)**2)-(cos(rmin*(kin1+kin2))/(kin1+kin2)**2))) + &
                       (B1*B2*((-1.0_dp*(exp(-1.0_dp*rmax*(kout1+kout2))*(kout1*rmax+kout2*rmax+1_dp)/(kout1+kout2)**2))+&
                       ((exp(-1.0_dp*r*(kout1+kout2))*(kout1*r+kout2*r+1._dp)/(kout1+kout2)**2))))))/Cm_to_D

end function TransDip_Ana

real(dp) function pulse(t)
real(dp) :: t

pulse = pulse1 * Ed01 * cos(omega01*(t-t01)+phase01) * exp(-1._dp*((t-t01)**2._dp/(2._dp*width01**2._dp))) + &
        pulse2 * Ed02 * cos(omega02*(t-t02)+phase02) * exp(-1._dp*((t-t02)**2._dp/(2._dp*width02**2._dp))) + &
        pulse3 * Ed03 * cos(omega03*(t-t03)+phase03) * exp(-1._dp*((t-t03)**2._dp/(2._dp*width03**2._dp)))

!write(6,*) Ed01*E_au,Ed02*E_au,Ed03*E_au
!write(6,*) omega01/(t_au*(2._dp*pi)),omega02/(t_au*(2._dp*pi)),omega03/(t_au*(2._dp*pi))
!write(6,*) phase01,phase02,phase03
!write(6,*) width01/(1.e-15_dp/t_au),width02/(1.e-15_dp/t_au),width03/(1.e-15_dp/t_au)

end function pulse

real(dp) function pulsex(t)
real(dp) :: t

pulsex= &
pulse1*Pe1(1)*Ed01*cos(-dot_product(k_1(:),Dcenter(:))+omega01*(t-t01)+phase01) &
                 * exp(-1._dp*((t-t01)**2._dp/(2._dp*width01**2._dp)))+&
pulse2*Pe2(1)*Ed02*cos(-dot_product(k_2(:),Dcenter(:))+omega02*(t-t02)+phase02) &
                 * exp(-1._dp*((t-t02)**2._dp/(2._dp*width02**2._dp)))+&
pulse3*Pe3(1)*Ed03*cos(-dot_product(k_3(:),Dcenter(:))+omega03*(t-t03)+phase03) &
                 * exp(-1._dp*((t-t03)**2._dp/(2._dp*width03**2._dp)))

end function pulsex

real(dp) function pulsey(t)
real(dp) :: t

pulsey=&
pulse1*Pe1(2)*Ed01* cos(-dot_product(k_1(:),Dcenter(:))+omega01*(t-t01)+phase01) & 
                  * exp(-1._dp*((t-t01)**2._dp/(2._dp*width01**2._dp)))+&
pulse2*Pe2(2)*Ed02* cos(-dot_product(k_2(:),Dcenter(:))+omega02*(t-t02)+phase02) & 
                  * exp(-1._dp*((t-t02)**2._dp/(2._dp*width02**2._dp)))+&
pulse3*Pe3(2)*Ed03* cos(-dot_product(k_3(:),Dcenter(:))+omega03*(t-t03)+phase03) & 
                  * exp(-1._dp*((t-t03)**2._dp/(2._dp*width03**2._dp)))

end function pulsey

real(dp) function pulsez(t)
real(dp) :: t

pulsez=&
pulse1*Pe1(3)*Ed01* cos(-dot_product(k_1(:),Dcenter(:))+omega01*(t-t01)+phase01) & 
                  * exp(-1._dp*((t-t01)**2._dp/(2._dp*width01**2._dp)))+&
pulse2*Pe2(3)*Ed02* cos(-dot_product(k_2(:),Dcenter(:))+omega02*(t-t02)+phase02) & 
                  * exp(-1._dp*((t-t02)**2._dp/(2._dp*width02**2._dp)))+&
pulse3*Pe3(3)*Ed03* cos(-dot_product(k_3(:),Dcenter(:))+omega03*(t-t03)+phase03) & 
                  * exp(-1._dp*((t-t03)**2._dp/(2._dp*width03**2._dp)))

end function pulsez

subroutine RK_0_ei

!!!Opens output files
if ( ( nofiles .eq. 'n' ) .and. ( Dyn_0 .eq. 'y' ) ) then
write(popc,'(a5,i5.5,a4)') 'Popc-', n, '.dat'
write(Re_c,'(a5,i5.5,a4)') 'Re_c-', n, '.dat'
write(Im_c,'(a5,i5.5,a4)') 'Im_c-', n, '.dat'
open(newunit=popc_0_f   ,file=popc)
open(newunit=Re_c_0_f   ,file=Re_c)
open(newunit=Im_c_0_f   ,file=Im_c)
endif
if ( ( nofiles .eq. 'n' ) .and. ( Dyn_ei .eq. 'y' ) ) then
write(popc_ei,'(a8,i5.5,a4)') 'Popc_ei-', n, '.dat'
write(Re_c_ei,'(a8,i5.5,a4)') 'Re_c_ei-', n, '.dat'
write(Im_c_ei,'(a8,i5.5,a4)') 'Im_c_ei-', n, '.dat'
open(newunit=popc_ei_f  ,file=popc_ei)
open(newunit=Re_c_ei_f  ,file=Re_c_ei)
open(newunit=Im_c_ei_f  ,file=Im_c_ei)
endif
if ( ( nofiles .eq. 'n' ) .and. ( Dyn_L .eq. 'y' ) ) then
write(Re_c_L,'(a7,i5.5,a4)') 'Re_c_L-', n, '.dat'
write(Im_c_L,'(a7,i5.5,a4)') 'Im_c_L-', n, '.dat'
open(newunit=Re_c_L_f   ,file=Re_c_L)
open(newunit=Im_c_L_f   ,file=Im_c_L)
endif

!!!!!INITIAL POPULATIONS
c0         = 0._dp
c0(0)      = 1._dp
xc0        = dcmplx(c0,0._dp)
xc         = dcmplx(0._dp,0._dp)
xc_ei      = dcmplx(0._dp,0._dp)
xc_L       = dcmplx(0._dp,0._dp)
xc(:,0)    = xc0(:)
xc_ei(:,0) = xc0(:)
xc_L(0,0)  = dcmplx(1._dp,0._dp)
powtemp    = 0._dp
f_ana = 0

do t=0,ntime

time = t*timestep

if ( n .eq. 1 ) then
pulses(t) = pulse1 * Ed01 * cos(omega01*(time-t01)+phase01) * exp(-(time-t01)**2/(2.0_dp*(width01**2._dp))) + &
            pulse2 * Ed02 * cos(omega02*(time-t02)+phase02) * exp(-(time-t02)**2/(2.0_dp*(width02**2._dp))) + &
            pulse3 * Ed03 * cos(omega03*(time-t03)+phase03) * exp(-(time-t03)**2/(2.0_dp*(width03**2._dp)))

!write(6,*) t, pulses(t)
endif

if ( Dyn_0 .eq. 'y' ) then

k1 = 0._dp ; k2 = 0._dp ; k3 = 0._dp ; k4 = 0._dp ; k5 = 0._dp ; k6 = 0._dp ; k7 = 0._dp ; k8 = 0._dp

if ( rdm_ori .eq. 'n' ) then 

tp1 = pulse(time)
tp2 = pulse(time+(timestep/9._dp))
tp3 = pulse(time+(timestep/6._dp))
tp4 = pulse(time+(timestep/3._dp))
tp5 = pulse(time+0.5_dp*timestep)
tp6 = pulse(time+(2._dp*timestep/3._dp))
tp7 = pulse(time+(5._dp*timestep/6._dp))
tp8 = pulse(time+timestep)

do i=0,nstates-1             
k1(i) = k1(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam(i,:)*tp1)* &
xc(:,t)) 
enddo

do i=0,nstates-1
k2(i) = k2(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam(i,:)*tp2)* &
(xc(:,t) + k1(:)/9._dp))
enddo

do i=0,nstates-1
k3(i) = k3(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam(i,:)*tp3)* &
(xc(:,t)+(k1(:)+3._dp*k2(:))/24._dp))
enddo

do i=0,nstates-1
k4(i) = k4(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam(i,:)*tp4)* &
(xc(:,t)+(k1(:)-3._dp*k2(:)+4._dp*k3(:))/6._dp))
enddo

do i=0,nstates-1
k5(i) = k5(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam(i,:)*tp5)* &
(xc(:,t)+(-5._dp*k1(:)+27._dp*k2(:)-24._dp*k3(:)+6._dp*k4(:))*0.125_dp))
enddo

do i=0,nstates-1
k6(i) = k6(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam(i,:)*tp6)* &
(xc(:,t)+(221._dp*k1(:)-981._dp*k2(:)+867._dp*k3(:)-102._dp*k4(:)+k5(:))/9._dp))
enddo

do i=0,nstates-1
k7(i) = k7(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam(i,:)*tp7)* &
(xc(:,t)+(-183._dp*k1(:)+678._dp*k2(:)-472._dp*k3(:)-66._dp*k4(:)+80._dp*k5(:)+3._dp*k6(:))/48._dp))
enddo

do i=0,nstates-1
k8(i) = k8(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam(i,:)*tp8)* &
(xc(:,t)+&
(716._dp*k1(:)-2079._dp*k2(:)+1002._dp*k3(:)+834._dp*k4(:)-454._dp*k5(:)-9._dp*k6(:)+72._dp*k7(:))/82._dp))
enddo

elseif ( rdm_ori .eq. 'y' ) then

tpx1 = pulsex(time)
tpx2 = pulsex(time+(timestep/9._dp))
tpx3 = pulsex(time+(timestep/6._dp))
tpx4 = pulsex(time+(timestep/3._dp))
tpx5 = pulsex(time+0.5_dp*timestep)
tpx6 = pulsex(time+(2._dp*timestep/3._dp))
tpx7 = pulsex(time+(5._dp*timestep/6._dp))
tpx8 = pulsex(time+timestep)

tpy1 = pulsey(time)
tpy2 = pulsey(time+(timestep/9._dp))
tpy3 = pulsey(time+(timestep/6._dp))
tpy4 = pulsey(time+(timestep/3._dp))
tpy5 = pulsey(time+0.5_dp*timestep)
tpy6 = pulsey(time+(2._dp*timestep/3._dp))
tpy7 = pulsey(time+(5._dp*timestep/6._dp))
tpy8 = pulsey(time+timestep)

tpz1 = pulsez(time)
tpz2 = pulsez(time+(timestep/9._dp))
tpz3 = pulsez(time+(timestep/6._dp))
tpz4 = pulsez(time+(timestep/3._dp))
tpz5 = pulsez(time+0.5_dp*timestep)
tpz6 = pulsez(time+(2._dp*timestep/3._dp))
tpz7 = pulsez(time+(5._dp*timestep/6._dp))
tpz8 = pulsez(time+timestep)

do i=0,nstates-1
k1(i) = k1(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam_l(i,:,1)*tpx1 - &
TransHam_l(i,:,2)*tpy1 - TransHam_l(i,:,3)*tpz1) * &
xc(:,t)) 
enddo

do i=0,nstates-1
k2(i) = k2(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam_l(i,:,1)*tpx2 - &
TransHam_l(i,:,2)*tpy2 - TransHam_l(i,:,3)*tpz2) * &
(xc(:,t) + k1(:)/9._dp))
enddo

do i=0,nstates-1
k3(i) = k3(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam_l(i,:,1)*tpx3 - &
TransHam_l(i,:,2)*tpy3 - TransHam_l(i,:,3)*tpz3) * &
(xc(:,t)+(k1(:)+3._dp*k2(:))/24._dp))
enddo

do i=0,nstates-1
k4(i) = k4(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam_l(i,:,1)*tpx4 - &
TransHam_l(i,:,2)*tpy4 - TransHam_l(i,:,3)*tpz4) * &
(xc(:,t)+(k1(:)-3._dp*k2(:)+4._dp*k3(:))/6._dp))
enddo

do i=0,nstates-1
k5(i) = k5(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam_l(i,:,1)*tpx5 - &
TransHam_l(i,:,2)*tpy5 - TransHam_l(i,:,3)*tpz5) * &
(xc(:,t)+(-5._dp*k1(:)+27._dp*k2(:)-24._dp*k3(:)+6._dp*k4(:))*0.125_dp))
enddo

do i=0,nstates-1
k6(i) = k6(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam_l(i,:,1)*tpx6 - &
TransHam_l(i,:,2)*tpy6 - TransHam_l(i,:,3)*tpz6) * &
(xc(:,t)+(221._dp*k1(:)-981._dp*k2(:)+867._dp*k3(:)-102._dp*k4(:)+k5(:))/9._dp))
enddo

do i=0,nstates-1
k7(i) = k7(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam_l(i,:,1)*tpx7 - &
TransHam_l(i,:,2)*tpy7 - TransHam_l(i,:,3)*tpz7) * &
(xc(:,t)+(-183._dp*k1(:)+678._dp*k2(:)-472._dp*k3(:)-66._dp*k4(:)+80._dp*k5(:)+3._dp*k6(:))/48._dp))
enddo

do i=0,nstates-1
k8(i) = k8(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham(i,:) - TransHam_l(i,:,1)*tpx8 - &
TransHam_l(i,:,2)*tpy8 - TransHam_l(i,:,3)*tpz8) * &
(xc(:,t)+&
(716._dp*k1(:)-2079._dp*k2(:)+1002._dp*k3(:)+834._dp*k4(:)-454._dp*k5(:)-9._dp*k6(:)+72._dp*k7(:))/82._dp))
enddo

endif

xc(:,t+1)=xc(:,t)+(41._dp*(k1(:)+k8(:))+216._dp*(k3(:)+k7(:))+27._dp*(k4(:)+k6(:))+272._dp*k5(:))/840._dp

if ( ( MOD(t,10) .eq. 0 ) .and. ( nofiles .eq. 'n' ) ) then
cnorm2 = sum(dreal(xc(:,t))**2 + dimag(xc(:,t))**2)
write(popc_0_f,form_pop) time*t_au, (dreal(xc(i,t))**2+dimag(xc(i,t))**2, i=0,nstates-1), cnorm2
write(Re_c_0_f,form_com) time*t_au, (dreal(xc(i,t)), i=0,nstates-1)
write(Im_c_0_f,form_com) time*t_au, (dimag(xc(i,t)), i=0,nstates-1)
endif

endif

if ( ( Dyn_ei .eq. 'y' ) .and. ( f_ana .eq. 0 ) ) then

k1 = 0._dp ; k2 = 0._dp ; k3 = 0._dp ; k4 = 0._dp ; k5 = 0._dp ; k6 = 0._dp ; k7 = 0._dp ; k8 = 0._dp

if ( rdm_ori .eq. 'n' ) then 

tp1 = pulse(time)
tp2 = pulse(time+(timestep/9._dp))
tp3 = pulse(time+(timestep/6._dp))
tp4 = pulse(time+(timestep/3._dp))
tp5 = pulse(time+0.5_dp*timestep)
tp6 = pulse(time+(2._dp*timestep/3._dp))
tp7 = pulse(time+(5._dp*timestep/6._dp))
tp8 = pulse(time+timestep)

do i=0,nstates-1             
k1(i) = k1(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei(i,:)*tp1)* &
xc_ei(:,t)) 
enddo

do i=0,nstates-1
k2(i) = k2(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei(i,:)*tp2)* &
(xc_ei(:,t) + k1(:)/9._dp))
enddo

do i=0,nstates-1
k3(i) = k3(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei(i,:)*tp3)* &
(xc_ei(:,t)+(k1(:)+3._dp*k2(:))/24._dp))
enddo

do i=0,nstates-1
k4(i) = k4(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei(i,:)*tp4)* &
(xc_ei(:,t)+(k1(:)-3._dp*k2(:)+4._dp*k3(:))/6._dp))
enddo

do i=0,nstates-1
k5(i) = k5(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei(i,:)*tp5)* &
(xc_ei(:,t)+(-5._dp*k1(:)+27._dp*k2(:)-24._dp*k3(:)+6._dp*k4(:))*0.125_dp))
enddo

do i=0,nstates-1
k6(i) = k6(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei(i,:)*tp6)* &
(xc_ei(:,t)+(221._dp*k1(:)-981._dp*k2(:)+867._dp*k3(:)-102._dp*k4(:)+k5(:))/9._dp))
enddo

do i=0,nstates-1
k7(i) = k7(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei(i,:)*tp7)* &
(xc_ei(:,t)+(-183._dp*k1(:)+678._dp*k2(:)-472._dp*k3(:)-66._dp*k4(:)+80._dp*k5(:)+3._dp*k6(:))/48._dp))
enddo

do i=0,nstates-1
k8(i) = k8(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei(i,:)*tp8)* &
(xc_ei(:,t)+&
(716._dp*k1(:)-2079._dp*k2(:)+1002._dp*k3(:)+834._dp*k4(:)-454._dp*k5(:)-9._dp*k6(:)+72._dp*k7(:))/82._dp))
enddo

elseif ( rdm_ori .eq. 'y' ) then

tpx1 = pulsex(time)
tpx2 = pulsex(time+(timestep/9._dp))
tpx3 = pulsex(time+(timestep/6._dp))
tpx4 = pulsex(time+(timestep/3._dp))
tpx5 = pulsex(time+0.5_dp*timestep)
tpx6 = pulsex(time+(2._dp*timestep/3._dp))
tpx7 = pulsex(time+(5._dp*timestep/6._dp))
tpx8 = pulsex(time+timestep)

tpy1 = pulsey(time)
tpy2 = pulsey(time+(timestep/9._dp))
tpy3 = pulsey(time+(timestep/6._dp))
tpy4 = pulsey(time+(timestep/3._dp))
tpy5 = pulsey(time+0.5_dp*timestep)
tpy6 = pulsey(time+(2._dp*timestep/3._dp))
tpy7 = pulsey(time+(5._dp*timestep/6._dp))
tpy8 = pulsey(time+timestep)

tpz1 = pulsez(time)
tpz2 = pulsez(time+(timestep/9._dp))
tpz3 = pulsez(time+(timestep/6._dp))
tpz4 = pulsez(time+(timestep/3._dp))
tpz5 = pulsez(time+0.5_dp*timestep)
tpz6 = pulsez(time+(2._dp*timestep/3._dp))
tpz7 = pulsez(time+(5._dp*timestep/6._dp))
tpz8 = pulsez(time+timestep)

do i=0,nstates-1
k1(i) = k1(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei_l(i,:,1)*tpx1 - &
TransHam_ei_l(i,:,2)*tpy1 - TransHam_ei_l(i,:,3)*tpz1) * &
xc_ei(:,t)) 
enddo

do i=0,nstates-1
k2(i) = k2(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei_l(i,:,1)*tpx2 - &
TransHam_ei_l(i,:,2)*tpy2 - TransHam_ei_l(i,:,3)*tpz2) * &
(xc_ei(:,t) + k1(:)/9._dp))
enddo

do i=0,nstates-1
k3(i) = k3(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei_l(i,:,1)*tpx3 - &
TransHam_ei_l(i,:,2)*tpy3 - TransHam_ei_l(i,:,3)*tpz3) * &
(xc_ei(:,t)+(k1(:)+3._dp*k2(:))/24._dp))
enddo

do i=0,nstates-1
k4(i) = k4(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei_l(i,:,1)*tpx4 - &
TransHam_ei_l(i,:,2)*tpy4 - TransHam_ei_l(i,:,3)*tpz4) * &
(xc_ei(:,t)+(k1(:)-3._dp*k2(:)+4._dp*k3(:))/6._dp))
enddo

do i=0,nstates-1
k5(i) = k5(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei_l(i,:,1)*tpx5 - &
TransHam_ei_l(i,:,2)*tpy5 - TransHam_ei_l(i,:,3)*tpz5) * &
(xc_ei(:,t)+(-5._dp*k1(:)+27._dp*k2(:)-24._dp*k3(:)+6._dp*k4(:))*0.125_dp))
enddo

do i=0,nstates-1
k6(i) = k6(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei_l(i,:,1)*tpx6 - &
TransHam_ei_l(i,:,2)*tpy6 - TransHam_ei_l(i,:,3)*tpz6) * &
(xc_ei(:,t)+(221._dp*k1(:)-981._dp*k2(:)+867._dp*k3(:)-102._dp*k4(:)+k5(:))/9._dp))
enddo

do i=0,nstates-1
k7(i) = k7(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei_l(i,:,1)*tpx7 - &
TransHam_ei_l(i,:,2)*tpy7 - TransHam_ei_l(i,:,3)*tpz7) * &
(xc_ei(:,t)+(-183._dp*k1(:)+678._dp*k2(:)-472._dp*k3(:)-66._dp*k4(:)+80._dp*k5(:)+3._dp*k6(:))/48._dp))
enddo

do i=0,nstates-1
k8(i) = k8(i) - dcmplx(timestep,0._dp) * sum(dcmplx(0._dp,Ham_l(i,:) - TransHam_ei_l(i,:,1)*tpx8 - &
TransHam_ei_l(i,:,2)*tpy8 - TransHam_ei_l(i,:,3)*tpz8) * &
(xc_ei(:,t)+&
(716._dp*k1(:)-2079._dp*k2(:)+1002._dp*k3(:)+834._dp*k4(:)-454._dp*k5(:)-9._dp*k6(:)+72._dp*k7(:))/82._dp))
enddo

endif

xc_ei(:,t+1)=xc_ei(:,t)+(41._dp*(k1(:)+k8(:))+216._dp*(k3(:)+k7(:))+27._dp*(k4(:)+k6(:))+272._dp*k5(:))/840._dp

if ( ( xc_ei(0,t+1) .eq. xc_ei(0,t) ) .and. ( time .gt. max(pulse1*t01,pulse2*t02,pulse3*t03) ) .and. ( f_ana .eq. 0 ) ) then
f_ana = 1
t_ana = t
time_ana = time
endif

if ( ( MOD(t,10) .eq. 0 ) .and. ( nofiles .eq. 'n' ) ) then
cnorm2 = sum(dreal(xc_ei(:,t))**2 + dimag(xc_ei(:,t))**2)
write(popc_ei_f,form_pop) time*t_au, (dreal(xc_ei(i,t))**2+dimag(xc_ei(i,t))**2, i=0,nstates-1), cnorm2
write(Re_c_ei_f,form_com) time*t_au, (dreal(xc_ei(i,t)), i=0,nstates-1)
write(Im_c_ei_f,form_com) time*t_au, (dimag(xc_ei(i,t)), i=0,nstates-1)
endif

elseif ( ( Dyn_ei .eq. 'y' ) .and. ( f_ana .eq. 1) ) then

xc_ei(:,t) = xc_ei(:,t_ana)*exp(-im*lambda(:)*(time-time_ana))

if ( ( MOD(t,10) .eq. 0 ) .and. ( nofiles .eq. 'n' ) ) then
cnorm2 = sum(dreal(xc_ei(:,t))**2 + dimag(xc_ei(:,t))**2)
write(popc_ei_f,form_pop) time*t_au, (dreal(xc_ei(i,t))**2+dimag(xc_ei(i,t))**2, i=0,nstates-1), cnorm2
write(Re_c_ei_f,form_com) time*t_au, (dreal(xc_ei(i,t)), i=0,nstates-1)
write(Im_c_ei_f,form_com) time*t_au, (dimag(xc_ei(i,t)), i=0,nstates-1)
endif

endif

if ( ( Dyn_L .eq. 'y' ) ) then 

tp1 = pulse(time)
tp2 = pulse(time+(timestep/9._dp))
tp3 = pulse(time+(timestep/6._dp))
tp4 = pulse(time+(timestep/3._dp))
tp5 = pulse(time+0.5_dp*timestep)
tp6 = pulse(time+(2._dp*timestep/3._dp))
tp7 = pulse(time+(5._dp*timestep/6._dp))
tp8 = pulse(time+timestep)

k1_L = 0._dp ; k2_L = 0._dp ; k3_L = 0._dp ; k4_L = 0._dp ; k5_L = 0._dp ; k6_L = 0._dp ; k7_L = 0._dp ; k8_L=0._dp 

do i=0,nstates2-1
k1_L(i) = k1_L(i) + dcmplx(timestep,0._dp) * sum((dcmplx(0._dp,-1._dp*lfield(i,:)*(merge_diag(i,:)-merge_odiag(i,:)*tp1)) &
- ( time - merge(t01,t02, ( time .le. t02 ) ) )*sigD(i,:) &
!- (time-t01)*sigD(i,:) &
) * xc_L(:,t))
enddo

do i=0,nstates2-1
k2_L(i) = k2_L(i) + dcmplx(timestep,0._dp) * sum((dcmplx(0._dp,-1._dp*lfield(i,:)*(merge_diag(i,:)-merge_odiag(i,:)*tp2)) &
!- (time-t01)*sigD(i,:) &
- ( time - merge(t01,t02, ( time .le. t02 ) ) )*sigD(i,:) &
) * (xc_L(:,t)+(1._dp/9._dp)*k1_L(:)))
enddo

do i=0,nstates2-1
k3_L(i) = k3_L(i) + dcmplx(timestep,0._dp) * sum((dcmplx(0._dp,-1._dp*lfield(i,:)*(merge_diag(i,:)-merge_odiag(i,:)*tp3)) &
!- (time-t01)*sigD(i,:) &
- ( time - merge(t01,t02, ( time .le. t02 ) ) )*sigD(i,:) &
) * (xc_L(:,t)+(1._dp/24._dp)*(k1_L(:)+3._dp*k2_L(:))))
enddo

do i=0,nstates2-1
k4_L(i) = k4_L(i) + dcmplx(timestep,0._dp) * sum((dcmplx(0._dp,-1._dp*lfield(i,:)*(merge_diag(i,:)-merge_odiag(i,:)*tp4)) &
!- (time-t01)*sigD(i,:) &
- ( time - merge(t01,t02, ( time .le. t02 ) ) )*sigD(i,:) &
) * (xc_L(:,t)+(1._dp/6._dp)*(k1_L(:)-3._dp*k2_L(:)+4._dp*k3_L(:))))
enddo

do i=0,nstates2-1
k5_L(i) = k5_L(i) + dcmplx(timestep,0._dp) * sum((dcmplx(0._dp,-1._dp*lfield(i,:)*(merge_diag(i,:)-merge_odiag(i,:)*tp5)) &
!- (time-t01)*sigD(i,:) &
- ( time - merge(t01,t02, ( time .le. t02 ) ) )*sigD(i,:) &
) * (xc_L(:,t)+(1._dp/8._dp)*(-5._dp*k1_L(:)+27._dp*k2_L(:)-24._dp*k3_L(:)+6._dp*k4_L(:))))
enddo

do i=0,nstates2-1
k6_L(i) = k6_L(i) + dcmplx(timestep,0._dp) * sum((dcmplx(0._dp,-1._dp*lfield(i,:)*(merge_diag(i,:)-merge_odiag(i,:)*tp6)) &
!- (time-t01)*sigD(i,:) &
- ( time - merge(t01,t02, ( time .le. t02 ) ) )*sigD(i,:) &
) * (xc_L(:,t)+(1._dp/9._dp)*(221._dp*k1_L(:)-981._dp*k2_L(:)+867._dp*k3_L(:)-102._dp*k4_L(:)+k5_L(:))))
enddo

do i=0,nstates2-1
k7_L(i) = k7_L(i) + dcmplx(timestep,0._dp) * sum((dcmplx(0._dp,-1._dp*lfield(i,:)*(merge_diag(i,:)-merge_odiag(i,:)*tp7)) &
!- (time-t01)*sigD(i,:) &
- ( time - merge(t01,t02, ( time .le. t02 ) ) )*sigD(i,:) &
) * (xc_L(:,t)+(1._dp/48._dp)*(-183._dp*k1_L(:)+678._dp*k2_L(:)-472._dp*k3_L(:)-66._dp*k4_L(:)+80._dp*k5_L(:)+3._dp*k6_L(:))))
enddo

do i=0,nstates2-1
k8_L(i) = k8_L(i) + dcmplx(timestep,0._dp) * sum((dcmplx(0._dp,-1._dp*lfield(i,:)*(merge_diag(i,:)-merge_odiag(i,:)*tp8)) &
- ( time - merge(t01,t02, ( time .le. t02 ) ) )*sigD(i,:) &
!- (time-t01)*sigD(i,:) &
) * (xc_L(:,t)+(1._dp/82._dp) * (716._dp*k1_L(:)-2079._dp*k2_L(:)+1002._dp*k3_L(:)+834._dp*k4_L(:)-454._dp*k5_L(:)-&
                9._dp*k6_L(:)+72._dp*k7_L(:))))
enddo

xc_L(:,t+1)=xc_L(:,t)+(41._dp*(k1_L(:)+k8_L(:))+216._dp*(k3_L(:)+k7_L(:))+27._dp*(k4_L(:)+k6_L(:))+272._dp*k5_L(:))&
              /840._dp

if ( ( MOD(t,10) .eq. 0 ) .and. ( nofiles .eq. 'n' ) ) then
cnorm2 = sum(dreal(xc_L(:,t))**2 + dimag(xc_L(:,t))**2)
write(Re_c_L_f,form_com_L) time*t_au, (dreal(xc_L(i,t)), i=0,nstates2-1)
write(Im_c_L_f,form_com_L) time*t_au, (dimag(xc_L(i,t)), i=0,nstates2-1)
endif

endif

do i=0,nstates-1
powtemp(t) = powtemp(t) + 2._dp * sum(TransHam_ei(i,:) * dreal(dconjg(xc_ei(i,t))*xc_ei(:,t)))
enddo

pow(t) = pow(t) + powtemp(t)

if ( inbox .eq. 'y' ) then
!do pol=1,npol
!pow_pol(pol,t) = pow_pol(pol,t) + dcmplx(pow_s(n,t),0.d0)*&
!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(pol)*Pe1(:)+l2(pol)*Pe2(:)+l3(pol)*Pe3(:),Dcenter(n,:)),0.e0_dp))
!enddo

!NORMAL POLARIZATION
!pow_pol(39,t) = pow_pol(39,t) + dcmplx(powtemp(t),0._dp)*&
!   exp(-im*dcmplx(dot_product(l1(39)*k_1(:)+l2(39)*k_2(:)+l3(39)*k_3(:),Dcenter(:)),0._dp))
!
!pow_pol(41,t) = pow_pol(41,t) + dcmplx(powtemp(t),0._dp)*&
!   exp(-im*dcmplx(dot_product(l1(41)*k_1(:)+l2(41)*k_2(:)+l3(41)*k_3(:),Dcenter(:)),0._dp))

!!FAR PROBE POLARIZATION 
pow_pol(39,t) = pow_pol(39,t) + dcmplx(powtemp(t),0._dp)*&
   exp(-im*dcmplx(dot_product(l1(39)*k_1(:)+l2(39)*k_2(:)+l3(39)*k_3(:),Dcenter(:)),0._dp))

pow_pol(41,t) = pow_pol(41,t) + dcmplx(powtemp(t),0._dp)*&
   exp(-im*dcmplx(dot_product(l1(41)*k_1(:)+l2(41)*k_2(:)+l3(41)*k_3(:),Dcenter(:)),0._dp))


!pow_pol(43,t) = pow_pol(43,t) + dcmplx(pow_s(n,t),0._dp)*&
!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(43)*Pe1(:)+l2(43)*Pe2(:)+l3(43)*Pe3(:),Dcenter(n,:)),0._dp))
!
!pow_pol(44,t) = pow_pol(44,t) + dcmplx(pow_s(n,t),0._dp)*&
!   exp(-im*(1._dp/545.e-9_dp)*dcmplx(dot_product(l1(44)*Pe1(:)+l2(44)*Pe2(:)+l3(44)*Pe3(:),Dcenter(n,:)),0._dp))

endif

enddo

if ( ( nofiles .eq. 'n' ) .and. ( Dyn_0 .eq. 'y' ) ) then
close(popc_0_f) ; close(Re_c_0_f) ; close(Im_c_0_f)
endif
if ( ( nofiles .eq. 'n' ) .and. ( Dyn_ei .eq. 'y' ) ) then
close(Re_c_ei_f) ; close(Im_c_ei_f) ; close(popc_ei_f)
endif
if ( ( nofiles .eq. 'n' ) .and. ( Dyn_L .eq. 'y' ) ) then
close(Re_c_L_f) ; close(Im_c_L_f)
endif

end subroutine RK_0_ei

subroutine Convolution

rewind Abs_imp_f
open(newunit=abso,file="Absorption.dat")

k=0

do
read(Abs_imp_f,*,iostat=io)
if (io .ne. 0) exit
k = k + 1
enddo

allocate(dipole(2,k))

rewind Abs_imp_f

do i=1,k
read(Abs_imp_f,*) dipole(1,i), dipole(2,i)
!print*, dipole(1,i), dipole(2,i)
enddo

sigma_conv=0.03d0
Emin=1.d0
Emax=7.d0
Estep=1000.d0

EminID = nint(Emin*Estep)
EmaxID = nint(Emax*Estep)

allocate(spec(EmaxID))

do j=EminID,EmaxID
do i=1,k
spec(j) = spec(j) + (dipole(2,i))**2 * exp(-1.d0*((j/1000.d0)-dipole(1,i))**2/(2*sigma_conv**2))
enddo
write(abso,*) j/1000.d0, spec(j), j
enddo

deallocate(spec,dipole)

end

recursive subroutine fft(x)
  complex(kind=dp), dimension(:), intent(inout)  :: x
  complex(kind=dp)                               :: t
  integer                                        :: N
  integer                                        :: i
  complex(kind=dp), dimension(:), allocatable    :: even, odd

  N=size(x)

  if(N .le. 1) return

  allocate(odd((N+1)/2))
  allocate(even(N/2))

  ! divide
  odd =x(1:N:2)
  even=x(2:N:2)

  ! conquer
  call fft(odd)
  call fft(even)

  ! combine
  do i=1,N/2
     t=exp(cmplx(0._dp,-2._dp*pi*real(i-1,dp)/real(N,dp),kind=dp))*even(i)
     x(i)     = odd(i) + t
     x(i+N/2) = odd(i) - t
  end do

  deallocate(odd)
  deallocate(even)

end subroutine fft

!real(dp) function E_particle(aR,mel,V0,n)
function E_particle(aR,mel,V0) result(minE)

integer :: je, jh, nsteps, i, n
real(dp) :: delta, Ef, aR, V0, mel
real(dp),dimension(10) :: minE
real(dp),allocatable :: diff(:), E(:)

delta=  0.00001e-18_dp
Ef=     1.28174e-18_dp
nsteps= int(Ef/delta)
!V0=V0/elec

allocate(diff(nsteps),E(nsteps),source=0._dp)

je = 1
diffe = 0._dp
do i=1,nsteps
E(i)=delta*real(i,dp)
diff(i) = abs(sqrt(2._dp*mel*E(i))/hbar * aR * 1._dp/tan(sqrt(2._dp*mel*E(i))/hbar*aR)- 1._dp + (mel/m0) + (mel*aR)/(hbar)&
           * sqrt((2._dp/m0)*(V0-E(i))))

if (diff(i) .le. diff(i-1)) then
        minE(je) = E(i)
elseif ( (diff(i) .ge. diff(i-1)) .and. (E(i-1) .eq. minE(je)) ) then
        je=je+1
endif
enddo

!return minE
!E_particle = minE(n)

!deallocate(diff,E,minE)

end function

subroutine get_kAB_A

!!wave vectors in and out
kineA=sqrt(2._dp*me*minEeA(1))/hbar
kouteA=sqrt(2._dp*m0*(V0eA-minEeA(1)))/hbar

kinh1A=sqrt(2._dp*mh*minEhA(1))/hbar
kouth1A=sqrt(2._dp*m0*(V0hA-minEhA(1)))/hbar

kinh2A=sqrt(2._dp*mh*minEhA(2))/hbar
kouth2A=sqrt(2._dp*m0*(V0hA-minEhA(2)))/hbar

!!normalization factors
AeA=1._dp/sqrt(aRA/2._dp-sin(2._dp*kineA*aRA)/(4._dp*kineA)+sin(kineA*aRA)**2/(2._dp*kouteA))
BeA=AeA*sin(kineA*aRA)*exp(kouteA*aRA)

Ah1A=1._dp/sqrt(aRA/2._dp-sin(2._dp*kinh1A*aRA)/(4._dp*kinh1A)+sin(kinh1A*aRA)**2/(2._dp*kouth1A))
Bh1A=Ah1A*sin(kinh1A*aRA)*exp(kouth1A*aRA)

Ah2A=1._dp/sqrt(aRA/2._dp-sin(2._dp*kinh2A*aRA)/(4._dp*kinh2A)+sin(kinh2A*aRA)**2/(2._dp*kouth2A))
Bh2A=Ah2A*sin(kinh2A*aRA)*exp(kouth2A*aRA)

end subroutine

subroutine get_kAB_B

!!wave vectors in and out
kineB=sqrt(2._dp*me*minEeB(1))/hbar
kouteB=sqrt(2._dp*m0*(V0eB-minEeB(1)))/hbar

kinh1B=sqrt(2._dp*mh*minEhB(1))/hbar
kouth1B=sqrt(2._dp*m0*(V0hB-minEhB(1)))/hbar

kinh2B=sqrt(2._dp*mh*minEhB(2))/hbar
kouth2B=sqrt(2._dp*m0*(V0hB-minEhB(2)))/hbar

!!normalization factors
AeB=1._dp/sqrt(aRB/2._dp-sin(2._dp*kineB*aRB)/(4._dp*kineB)+sin(kineB*aRB)**2/(2._dp*kouteB))
BeB=AeB*sin(kineB*aRB)*exp(kouteB*aRB)

Ah1B=1._dp/sqrt(aRB/2._dp-sin(2._dp*kinh1B*aRB)/(4._dp*kinh1B)+sin(kinh1B*aRB)**2/(2._dp*kouth1B))
Bh1B=Ah1B*sin(kinh1B*aRB)*exp(kouth1B*aRB)

Ah2B=1._dp/sqrt(aRB/2._dp-sin(2._dp*kinh2B*aRB)/(4._dp*kinh2B)+sin(kinh2B*aRB)**2/(2._dp*kouth2B))
Bh2B=Ah2B*sin(kinh2B*aRB)*exp(kouth2B*aRB)

end subroutine

subroutine QDA

V0eA   = -elec*(-3.49_dp+2.47_dp*(2.e9_dp*aRA)**(-1.32_dp))
V0hA   = -elec*(-5.23_dp-0.74_dp*(2.e9_dp*aRA)**(-0.95_dp))

minEeA = E_particle(aRA,me,V0eA)
minEhA = E_particle(aRA,mh,V0hA)

call get_kAB_A

TransDip_Ana_h1eA  = abs(TransDip_Ana(AeA,Ah1A,BeA,Bh1A,kineA,kinh1A,kouteA,kouth1A,aRA))
TransDip_Ana_h1h2A = abs(TransDip_Ana(Ah1A,Ah2A,Bh1A,Bh2A,kinh1A,kinh2A,kouth1A,kouth2A,aRA))
TransDip_Ana_h2eA  = abs(TransDip_Ana(AeA,Ah2A,BeA,Bh2A,kineA,kinh2A,kouteA,kouth2A,aRA))

end subroutine

subroutine QDB

V0eB   = -elec*(-3.49_dp+2.47_dp*(2.e9_dp*aRB)**(-1.32_dp))
V0hB   = -elec*(-5.23_dp-0.74_dp*(2.e9_dp*aRB)**(-0.95_dp))

minEeB = E_particle(aRB,me,V0eB)
minEhB = E_particle(aRB,mh,V0hB)

call get_kAB_B

TransDip_Ana_h1eB  = abs(TransDip_Ana(AeB,Ah1B,BeB,Bh1B,kineB,kinh1B,kouteB,kouth1B,aRB))
TransDip_Ana_h1h2B = abs(TransDip_Ana(Ah1B,Ah2B,Bh1B,Bh2B,kinh1B,kinh2B,kouth1B,kouth2B,aRB))
TransDip_Ana_h2eB  = abs(TransDip_Ana(AeB,Ah2B,BeB,Bh2B,kineB,kinh2B,kouteB,kouth2B,aRB))

end subroutine

end module Integrals
