module Make_Ham

use Constants_au
use Variables_au
use Integrals
use Normal

real(dp), parameter :: a11_1d_ho = 0.0191697_dp  
real(dp), parameter :: a11_2d_ho = 0.242775_dp  
real(dp), parameter :: a11_3d_ho = 1.28322_dp 
real(dp), parameter :: a11_1e_ho = 0.0131599_dp       
real(dp), parameter :: a11_2e_ho = 0.199444_dp        
real(dp), parameter :: a11_3e_ho = 1.10133_dp
real(dp), parameter :: a22_1d_ho = 0.0178298_dp  
real(dp), parameter :: a22_2d_ho = 0.239631_dp   
real(dp), parameter :: a22_3d_ho = 1.25828_dp
real(dp), parameter :: a22_1e_ho = 0.00477726_dp      
real(dp), parameter :: a22_2e_ho = 0.0379395_dp       
real(dp), parameter :: a22_3e_ho = 1.66235_dp
real(dp), parameter :: a12_1d_ho = -0.00288509_dp
real(dp), parameter :: a12_2d_ho = 0.0260066_dp
real(dp), parameter :: a12_3d_ho = 0.57151_dp
real(dp), parameter :: a12_1e_ho = -0.00651569_dp  
real(dp), parameter :: a12_2e_ho = 0.0530662_dp     
real(dp), parameter :: a12_3e_ho = 1.71364_dp

contains

subroutine make_Ham_he

if ( idlink .eq. 20 ) then
include 'Parameters-dir-02.f90'
include 'Parameters-ex-02.f90'
elseif ( idlink .eq. 55 ) then
include 'Parameters-dir-55.f90'
include 'Parameters-ex-55.f90'
endif

Ham_0(1)     = minEeA(1) + minEhA(1)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aRA*1e9_dp)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aRA*1e9_dp)**a11_3e_ho))

Ham_0(2)     = minEeA(1) + minEhA(2) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aRA*1e9_dp)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aRA*1e9_dp)**a22_3e_ho))

Ham_0(3)     = minEeB(1) + minEhB(1) + V0 
Ham_dir(3,3) = elec*(a11_1d_ho + a11_2d_ho / ((aRB*1e9_dp)**a11_3d_ho))
Ham_ex(3,3)  = elec*(a11_1e_ho + a11_2e_ho / ((aRB*1e9_dp)**a11_3e_ho))

Ham_0(4)     = minEeB(1) + minEhB(2) + V0 
Ham_dir(4,4) = elec*(a22_1d_ho + a22_2d_ho / ((aRB*1e9_dp)**a22_3d_ho))
Ham_ex(4,4)  = elec*(a22_1e_ho + a22_2e_ho / ((aRB*1e9_dp)**a22_3e_ho))

Ham_0(5)     = minEeB(1) + minEhA(1) + V0 
Ham_dir(5,5) = elec*(a55_1d_he + a55_2d_he / ((aRA*1e9_dp)**a55_3d_he * (aRB*1d9)**a55_4d_he)) 
Ham_ex(5,5)  = elec*(a55_1e_he + a55_2e_he / ((aRA*1e9_dp)**a55_3e_he * (aRB*1d9)**a55_4e_he))
                                  
Ham_0(6)     = minEeB(1) + minEhA(2) + V0 
Ham_dir(6,6) = elec*(a66_1d_he + a66_2d_he / ((aRA*1e9_dp)**a66_3d_he * (aRB*1d9)**a66_4d_he)) 
Ham_ex(6,6)  = elec*(a66_1e_he + a66_2e_he / ((aRA*1e9_dp)**a66_3e_he * (aRB*1d9)**a66_4e_he)) 
                                  
Ham_0(7)     = minEeA(1) + minEhB(1) + V0 
Ham_dir(7,7) = elec*(a55_1d_he + a55_2d_he / ((aRB*1e9_dp)**a55_3d_he * (aRA*1d9)**a55_4d_he)) 
Ham_ex(7,7)  = elec*(a55_1e_he + a55_2e_he / ((aRB*1e9_dp)**a55_3e_he * (aRA*1d9)**a55_4e_he))

Ham_0(8)     = minEeA(1) + minEhB(2) + V0 
Ham_dir(8,8) = elec*(a66_1d_he + a66_2d_he / ((aRB*1e9_dp)**a66_3d_he * (aRA*1d9)**a66_4d_he)) 
Ham_ex(8,8)  = elec*(a66_1e_he + a66_2e_he / ((aRB*1e9_dp)**a66_3e_he * (aRA*1d9)**a66_4e_he))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aRA*1e9_dp)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aRA*1e9_dp)**a12_3e_ho))

Ham_dir(1,3) = elec*(a13_1d_he + a13_2d_he / ((aRA*1e9_dp)**a13_3d_he * (aRB*1d9)**a13_4d_he))
Ham_ex(1,3)  = elec*(a13_1e_he + a13_2e_he / ((aRA*1e9_dp)**a13_3e_he * (aRB*1d9)**a13_4e_he))

Ham_dir(1,4) = elec*(a14_1d_he + a14_2d_he / ((aRA*1e9_dp)**a14_3d_he * (aRB*1d9)**a14_4d_he))
Ham_ex(1,4)  = elec*(a14_1e_he + a14_2e_he / ((aRA*1e9_dp)**a14_3e_he * (aRB*1d9)**a14_4e_he))

Ham_dir(1,5) = elec*(a15_1d_he + a15_2d_he / ((aRA*1e9_dp)**a15_3d_he * (aRB*1d9)**a15_4d_he))
Ham_ex(1,5)  = elec*(a15_1e_he + a15_2e_he / ((aRA*1e9_dp)**a15_3e_he * (aRB*1d9)**a15_4e_he))

Ham_dir(1,6) = elec*(a16_1d_he + a16_2d_he / ((aRA*1e9_dp)**a16_3d_he * (aRB*1d9)**a16_4d_he))
Ham_ex(1,6)  = elec*(a16_1e_he + a16_2e_he / ((aRA*1e9_dp)**a16_3e_he * (aRB*1d9)**a16_4e_he))

Ham_dir(1,7) = elec*(a17_1d_he + a17_2d_he / ((aRA*1e9_dp)**a17_3d_he * (aRB*1d9)**a17_4d_he))
Ham_ex(1,7)  = elec*(a17_1e_he + a17_2e_he / ((aRA*1e9_dp)**a17_3e_he * (aRB*1d9)**a17_4e_he))

Ham_dir(1,8) = elec*(a18_1d_he + a18_2d_he / ((aRA*1e9_dp)**a18_3d_he * (aRB*1d9)**a18_4d_he))
Ham_ex(1,8)  = elec*(a18_1e_he + a18_2e_he / ((aRA*1e9_dp)**a18_3e_he * (aRB*1d9)**a18_4e_he))

Ham_dir(2,3) = elec*(a14_1d_he + a14_2d_he / ((aRB*1e9_dp)**a14_3d_he * (aRA*1d9)**a14_4d_he))
Ham_ex(2,3)  = elec*(a14_1e_he + a14_2e_he / ((aRB*1e9_dp)**a14_3e_he * (aRA*1d9)**a14_4e_he))

Ham_dir(2,4) = elec*(a24_1d_he + a24_2d_he / ((aRA*1e9_dp)**a24_3d_he * (aRB*1d9)**a24_4d_he))
Ham_ex(2,4)  = elec*(a24_1e_he + a24_2e_he / ((aRA*1e9_dp)**a24_3e_he * (aRB*1d9)**a24_4e_he))

Ham_dir(2,5) = elec*(a25_1d_he + a25_2d_he / ((aRA*1e9_dp)**a25_3d_he * (aRB*1d9)**a25_4d_he))
Ham_ex(2,5)  = elec*(a25_1e_he + a25_2e_he / ((aRA*1e9_dp)**a25_3e_he * (aRB*1d9)**a25_4e_he))

Ham_dir(2,6) = elec*(a26_1d_he + a26_2d_he / ((aRA*1e9_dp)**a26_3d_he * (aRB*1d9)**a26_4d_he))
Ham_ex(2,6)  = elec*(a26_1e_he + a26_2e_he / ((aRA*1e9_dp)**a26_3e_he * (aRB*1d9)**a26_4e_he))

Ham_dir(2,7) = elec*(a27_1d_he + a27_2d_he / ((aRA*1e9_dp)**a27_3d_he * (aRB*1d9)**a27_4d_he))
Ham_ex(2,7)  = elec*(a27_1e_he + a27_2e_he / ((aRA*1e9_dp)**a27_3e_he * (aRB*1d9)**a27_4e_he))

Ham_dir(2,8) = elec*(a28_1d_he + a28_2d_he / ((aRA*1e9_dp)**a28_3d_he * (aRB*1d9)**a28_4d_he))
Ham_ex(2,8)  = elec*(a28_1e_he + a28_2e_he / ((aRA*1e9_dp)**a28_3e_he * (aRB*1d9)**a28_4e_he))

Ham_dir(3,4) = elec*(a12_1d_ho + a12_2d_ho / ((aRB*1e9_dp)**a12_3d_ho)) 
Ham_ex(3,4)  = elec*(a12_1e_ho + a12_2e_ho / ((aRB*1e9_dp)**a12_3e_ho))

Ham_dir(3,5) = elec*(a17_1d_he + a17_2d_he / ((aRB*1e9_dp)**a17_3d_he * (aRA*1d9)**a17_4d_he)) 
Ham_ex(3,5)  = elec*(a17_1e_he + a17_2e_he / ((aRB*1e9_dp)**a17_3e_he * (aRA*1d9)**a17_4e_he))

Ham_dir(3,6) = elec*(a18_1d_he + a18_2d_he / ((aRB*1e9_dp)**a18_3d_he * (aRA*1d9)**a18_4d_he)) 
Ham_ex(3,6)  = elec*(a18_1e_he + a18_2e_he / ((aRB*1e9_dp)**a18_3e_he * (aRA*1d9)**a18_4e_he))

Ham_dir(3,7) = elec*(a15_1d_he + a15_2d_he / ((aRB*1e9_dp)**a15_3d_he * (aRA*1d9)**a15_4d_he)) 
Ham_ex(3,7)  = elec*(a15_1e_he + a15_2e_he / ((aRB*1e9_dp)**a15_3e_he * (aRA*1d9)**a15_4e_he))

Ham_dir(3,8) = elec*(a16_1d_he + a16_2d_he / ((aRB*1e9_dp)**a16_3d_he * (aRA*1d9)**a16_4d_he)) 
Ham_ex(3,8)  = elec*(a16_1e_he + a16_2e_he / ((aRB*1e9_dp)**a16_3e_he * (aRA*1d9)**a16_4e_he))

Ham_dir(4,5) = elec*(a27_1d_he + a27_2d_he / ((aRB*1e9_dp)**a27_3d_he * (aRA*1d9)**a27_4d_he)) 
Ham_ex(4,5)  = elec*(a27_1e_he + a27_2e_he / ((aRB*1e9_dp)**a27_3e_he * (aRA*1d9)**a27_4e_he))

Ham_dir(4,6) = elec*(a28_1d_he + a28_2d_he / ((aRB*1e9_dp)**a28_3d_he * (aRA*1d9)**a28_4d_he)) 
Ham_ex(4,6)  = elec*(a28_1e_he + a28_2e_he / ((aRB*1e9_dp)**a28_3e_he * (aRA*1d9)**a28_4e_he))

Ham_dir(4,7) = elec*(a25_1d_he + a25_2d_he / ((aRB*1e9_dp)**a25_3d_he * (aRA*1d9)**a25_4d_he)) 
Ham_ex(4,7)  = elec*(a25_1e_he + a25_2e_he / ((aRB*1e9_dp)**a25_3e_he * (aRA*1d9)**a25_4e_he))

Ham_dir(4,8) = elec*(a26_1d_he + a26_2d_he / ((aRB*1e9_dp)**a26_3d_he * (aRA*1d9)**a26_4d_he)) 
Ham_ex(4,8)  = elec*(a26_1e_he + a26_2e_he / ((aRB*1e9_dp)**a26_3e_he * (aRA*1d9)**a26_4e_he))

Ham_dir(5,6) = elec*(a56_1d_he + a56_2d_he / ((aRA*1e9_dp)**a56_3d_he * (aRB*1d9)**a56_4d_he)) 
Ham_ex(5,6)  = elec*(a56_1e_he + a56_2e_he / ((aRA*1e9_dp)**a56_3e_he * (aRB*1d9)**a56_4e_he))

Ham_dir(5,7) = elec*(a57_1d_he + a57_2d_he / ((aRA*1e9_dp)**a57_3d_he * (aRB*1d9)**a57_4d_he))  
Ham_ex(5,7)  = elec*(a57_1e_he + a57_2e_he / ((aRA*1e9_dp)**a57_3e_he * (aRB*1d9)**a57_4e_he))   

Ham_dir(5,8) = elec*(a58_1d_he + a58_2d_he / ((aRA*1e9_dp)**a58_3d_he * (aRB*1d9)**a58_4d_he))  
Ham_ex(5,8)  = elec*(a58_1e_he + a58_2e_he / ((aRA*1e9_dp)**a58_3e_he * (aRB*1d9)**a58_4e_he))   

Ham_dir(6,7) = elec*(a58_1d_he + a58_2d_he / ((aRB*1e9_dp)**a58_3d_he * (aRA*1d9)**a58_4d_he))  
Ham_ex(6,7)  = elec*(a58_1e_he + a58_2e_he / ((aRB*1e9_dp)**a58_3e_he * (aRA*1d9)**a58_4e_he))   

Ham_dir(6,8) = elec*(a68_1d_he + a68_2d_he / ((aRA*1e9_dp)**a68_3d_he * (aRB*1d9)**a68_4d_he)) 
Ham_ex(6,8)  = elec*(a68_1e_he + a68_2e_he / ((aRA*1e9_dp)**a68_3e_he * (aRB*1d9)**a68_4e_he))   

Ham_dir(7,8) = elec*(a56_1d_he + a56_2d_he / ((aRB*1e9_dp)**a56_3d_he * (aRA*1d9)**a56_4d_he)) 
Ham_ex(7,8)  = elec*(a78_1e_he + a78_2e_he / ((aRB*1e9_dp)**a78_3e_he * (aRA*1d9)**a78_4e_he))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1._dp * Ham_dir(i,j) + Ham_ex(i,j)
    endif
  enddo
enddo

do i=1,nstates-1
  do j=i+1,nstates-1
    Ham(j,i) = Ham(i,j)
    Ham_dir(j,i) = Ham_dir(i,j)
    Ham_ex(j,i) = Ham_ex(i,j)
  enddo
enddo

Ham = Ham/Energ_au
Ham_dir = Ham_dir/Energ_au
Ham_ex = Ham_ex/Energ_au

if ( inbox .eq. "y" ) then

if ( rdm_ori .eq. "y" ) then

call random_number(psirot)
call random_number(phirot)
call random_number(thetarot)

psirot   = psirot   * 2._dp*pi
phirot   = phirot   * 2._dp*pi
thetarot = thetarot * 2._dp*pi

rotmat(1,1)= cos(phirot)*cos(thetarot)*cos(psirot) - sin(phirot)*sin(psirot)
rotmat(1,2)=-cos(phirot)*cos(thetarot)*sin(psirot) - sin(phirot)*cos(psirot)
rotmat(1,3)= cos(phirot)*sin(thetarot)
rotmat(2,1)= sin(phirot)*cos(thetarot)*cos(psirot) + cos(phirot)*sin(psirot)
rotmat(2,2)=-sin(phirot)*cos(thetarot)*sin(psirot) + cos(phirot)*cos(psirot)
rotmat(2,3)= sin(phirot)*sin(thetarot)
rotmat(3,1)=-sin(thetarot)*cos(psirot)
rotmat(3,2)= sin(thetarot)*sin(psirot)
rotmat(3,3)= cos(thetarot)

TransHam_l(0,1,:) = matmul(rotmat,eTDM(0,1,:))*(TransDip_Ana_h1eA)
TransHam_l(0,2,:) = matmul(rotmat,eTDM(0,2,:))*(TransDip_Ana_h2eA)
TransHam_l(0,3,:) = matmul(rotmat,eTDM(0,3,:))*(TransDip_Ana_h1eB)
TransHam_l(0,4,:) = matmul(rotmat,eTDM(0,4,:))*(TransDip_Ana_h2eB)
TransHam_l(0,5,:) = matmul(rotmat,eTDM(0,5,:))*(TransDip_Fit_h1e_he(aRB,aRA))
TransHam_l(0,6,:) = matmul(rotmat,eTDM(0,6,:))*(TransDip_Fit_h2e_he(aRB,aRA))
TransHam_l(0,7,:) = matmul(rotmat,eTDM(0,7,:))*(TransDip_Fit_h1e_he(aRA,aRB))
TransHam_l(0,8,:) = matmul(rotmat,eTDM(0,8,:))*(TransDip_Fit_h2e_he(aRA,aRB))
TransHam_l(1,2,:) = matmul(rotmat,eTDM(1,2,:))*(TransDip_Ana_h1h2A)
TransHam_l(1,5,:) = matmul(rotmat,eTDM(1,5,:))*(TransDip_Fit_ee_he(aRB,aRA))
TransHam_l(1,7,:) = matmul(rotmat,eTDM(1,7,:))*(TransDip_Fit_h1h1_he(aRB,aRA))
TransHam_l(1,8,:) = matmul(rotmat,eTDM(1,8,:))*(TransDip_Fit_h1h2_he(aRA,aRB))
TransHam_l(2,6,:) = matmul(rotmat,eTDM(2,6,:))*(TransDip_Fit_ee_he(aRB,aRA))
TransHam_l(2,7,:) = matmul(rotmat,eTDM(2,7,:))*(TransDip_Fit_h1h2_he(aRB,aRA))
TransHam_l(2,8,:) = matmul(rotmat,eTDM(2,8,:))*(TransDip_Fit_h2h2_he(aRB,aRA))
TransHam_l(3,4,:) = matmul(rotmat,eTDM(3,4,:))*(TransDip_Ana_h1h2B)
TransHam_l(3,5,:) = matmul(rotmat,eTDM(3,5,:))*(TransDip_Fit_h1h1_he(aRB,aRA))
TransHam_l(3,6,:) = matmul(rotmat,eTDM(3,6,:))*(TransDip_Fit_h1h2_he(aRB,aRA))
TransHam_l(3,7,:) = matmul(rotmat,eTDM(3,7,:))*(TransDip_Fit_ee_he(aRB,aRA))
TransHam_l(4,5,:) = matmul(rotmat,eTDM(4,5,:))*(TransDip_Fit_h1h2_he(aRA,aRB))
TransHam_l(4,6,:) = matmul(rotmat,eTDM(4,6,:))*(TransDip_Fit_h2h2_he(aRB,aRA))
TransHam_l(4,8,:) = matmul(rotmat,eTDM(4,8,:))*(TransDip_Fit_ee_he(aRB,aRA))
TransHam_l(5,6,:) = matmul(rotmat,eTDM(5,6,:))*(TransDip_Ana_h1h2A)
TransHam_l(7,8,:) = matmul(rotmat,eTDM(7,8,:))*(TransDip_Ana_h1h2B)

elseif ( rdm_ori .eq. "n" ) then

TransHam_l(0,1,:) = eTDM(0,1,:)*(TransDip_Ana_h1eA)
TransHam_l(0,2,:) = eTDM(0,2,:)*(TransDip_Ana_h2eA)
TransHam_l(0,3,:) = eTDM(0,3,:)*(TransDip_Ana_h1eB)
TransHam_l(0,4,:) = eTDM(0,4,:)*(TransDip_Ana_h2eB)
TransHam_l(0,5,:) = eTDM(0,5,:)*(TransDip_Fit_h1e_he(aRB,aRA))
TransHam_l(0,6,:) = eTDM(0,6,:)*(TransDip_Fit_h2e_he(aRB,aRA))
TransHam_l(0,7,:) = eTDM(0,7,:)*(TransDip_Fit_h1e_he(aRA,aRB))
TransHam_l(0,8,:) = eTDM(0,8,:)*(TransDip_Fit_h2e_he(aRA,aRB))
TransHam_l(1,2,:) = eTDM(1,2,:)*(TransDip_Ana_h1h2A)
TransHam_l(1,5,:) = eTDM(1,5,:)*(TransDip_Fit_ee_he(aRB,aRA))
TransHam_l(1,7,:) = eTDM(1,7,:)*(TransDip_Fit_h1h1_he(aRB,aRA))
TransHam_l(1,8,:) = eTDM(1,8,:)*(TransDip_Fit_h1h2_he(aRA,aRB))
TransHam_l(2,6,:) = eTDM(2,6,:)*(TransDip_Fit_ee_he(aRB,aRA))
TransHam_l(2,7,:) = eTDM(2,7,:)*(TransDip_Fit_h1h2_he(aRB,aRA))
TransHam_l(2,8,:) = eTDM(2,8,:)*(TransDip_Fit_h2h2_he(aRB,aRA))
TransHam_l(3,4,:) = eTDM(3,4,:)*(TransDip_Ana_h1h2B)
TransHam_l(3,5,:) = eTDM(3,5,:)*(TransDip_Fit_h1h1_he(aRB,aRA))
TransHam_l(3,6,:) = eTDM(3,6,:)*(TransDip_Fit_h1h2_he(aRB,aRA))
TransHam_l(3,7,:) = eTDM(3,7,:)*(TransDip_Fit_ee_he(aRB,aRA))
TransHam_l(4,5,:) = eTDM(4,5,:)*(TransDip_Fit_h1h2_he(aRA,aRB))
TransHam_l(4,6,:) = eTDM(4,6,:)*(TransDip_Fit_h2h2_he(aRB,aRA))
TransHam_l(4,8,:) = eTDM(4,8,:)*(TransDip_Fit_ee_he(aRB,aRA))
TransHam_l(5,6,:) = eTDM(5,6,:)*(TransDip_Ana_h1h2A)
TransHam_l(7,8,:) = eTDM(7,8,:)*(TransDip_Ana_h1h2B)

endif

do i=0,nstates-1
do j=i+1,nstates-1
TransHam_l(j,i,:) = TransHam_l(i,j,:)
enddo
enddo

TransHam_l = TransHam_l/D_to_au

elseif ( inbox .eq. "n" ) then

TransHam(0,1) = TransDip_Ana_h1eA
TransHam(0,2) = TransDip_Ana_h2eA
TransHam(0,3) = TransDip_Ana_h1eB
TransHam(0,4) = TransDip_Ana_h2eB
TransHam(0,5) = TransDip_Fit_h1e_he(aRB,aRA)
TransHam(0,6) = TransDip_Fit_h2e_he(aRB,aRA)
TransHam(0,7) = TransDip_Fit_h1e_he(aRA,aRB)
TransHam(0,8) = TransDip_Fit_h2e_he(aRA,aRB)
TransHam(1,2) = TransDip_Ana_h1h2A
TransHam(1,5) = TransDip_Fit_ee_he(aRB,aRA)
TransHam(1,7) = TransDip_Fit_h1h1_he(aRB,aRA)
TransHam(1,8) = TransDip_Fit_h1h2_he(aRA,aRB)
TransHam(2,6) = TransDip_Fit_ee_he(aRB,aRA)
TransHam(2,7) = TransDip_Fit_h1h2_he(aRB,aRA)
TransHam(2,8) = TransDip_Fit_h2h2_he(aRB,aRA)
TransHam(3,4) = TransDip_Ana_h1h2B
TransHam(3,5) = TransDip_Fit_h1h1_he(aRB,aRA)
TransHam(3,6) = TransDip_Fit_h1h2_he(aRB,aRA)
TransHam(3,7) = TransDip_Fit_ee_he(aRB,aRA)
TransHam(4,5) = TransDip_Fit_h1h2_he(aRA,aRB)
TransHam(4,6) = TransDip_Fit_h2h2_he(aRB,aRA)
TransHam(4,8) = TransDip_Fit_ee_he(aRB,aRA)
TransHam(5,6) = TransDip_Ana_h1h2A
TransHam(7,8) = TransDip_Ana_h1h2B

do i=0,nstates-1
do j=i+1,nstates-1
TransHam(j,i) = TransHam(i,j)
enddo
enddo

TransHam = TransHam/D_to_au

endif

end subroutine make_Ham_he

subroutine make_Ham_he_FO

if ( idlink .eq. 20 ) then
include 'Parameters-dir-02.f90'
include 'Parameters-ex-02.f90'
elseif ( idlink .eq. 55 ) then
include 'Parameters-dir-55.f90'
include 'Parameters-ex-55.f90'
endif

Ham_0(1)     = minEeA(1) + minEhA(1)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aRA*1e9_dp)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aRA*1e9_dp)**a11_3e_ho))

Ham_0(2)     = minEeA(1) + minEhA(2) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aRA*1e9_dp)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aRA*1e9_dp)**a22_3e_ho))

Ham_0(3)     = minEeB(1) + minEhB(1) + V0 
Ham_dir(3,3) = elec*(a11_1d_ho + a11_2d_ho / ((aRB*1e9_dp)**a11_3d_ho))
Ham_ex(3,3)  = elec*(a11_1e_ho + a11_2e_ho / ((aRB*1e9_dp)**a11_3e_ho))

Ham_0(4)     = minEeB(1) + minEhB(2) + V0 
Ham_dir(4,4) = elec*(a22_1d_ho + a22_2d_ho / ((aRB*1e9_dp)**a22_3d_ho))
Ham_ex(4,4)  = elec*(a22_1e_ho + a22_2e_ho / ((aRB*1e9_dp)**a22_3e_ho))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aRA*1e9_dp)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aRA*1e9_dp)**a12_3e_ho))

Ham_dir(1,3) = elec*(a13_1d_he + a13_2d_he / ((aRA*1e9_dp)**a13_3d_he * (aRB*1e9_dp)**a13_4d_he))
Ham_ex(1,3)  = elec*(a13_1e_he + a13_2e_he / ((aRA*1e9_dp)**a13_3e_he * (aRB*1e9_dp)**a13_4e_he))

Ham_dir(1,4) = elec*(a14_1d_he + a14_2d_he / ((aRA*1e9_dp)**a14_3d_he * (aRB*1e9_dp)**a14_4d_he))
Ham_ex(1,4)  = elec*(a14_1e_he + a14_2e_he / ((aRA*1e9_dp)**a14_3e_he * (aRB*1e9_dp)**a14_4e_he))

Ham_dir(2,3) = elec*(a14_1d_he + a14_2d_he / ((aRB*1e9_dp)**a14_3d_he * (aRA*1e9_dp)**a14_4d_he))
Ham_ex(2,3)  = elec*(a14_1e_he + a14_2e_he / ((aRB*1e9_dp)**a14_3e_he * (aRA*1e9_dp)**a14_4e_he))

Ham_dir(2,4) = elec*(a24_1d_he + a24_2d_he / ((aRA*1e9_dp)**a24_3d_he * (aRB*1e9_dp)**a24_4d_he))
Ham_ex(2,4)  = elec*(a24_1e_he + a24_2e_he / ((aRA*1e9_dp)**a24_3e_he * (aRB*1e9_dp)**a24_4e_he))

Ham_dir(3,4) = elec*(a12_1d_ho + a12_2d_ho / ((aRB*1e9_dp)**a12_3d_ho)) 
Ham_ex(3,4)  = elec*(a12_1e_ho + a12_2e_ho / ((aRB*1e9_dp)**a12_3e_ho))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1._dp * Ham_dir(i,j) + Ham_ex(i,j)
    endif
  enddo
enddo

do i=1,nstates-1
  do j=i+1,nstates-1
    Ham(j,i) = Ham(i,j)
    Ham_dir(j,i) = Ham_dir(i,j)
    Ham_ex(j,i) = Ham_ex(i,j)
  enddo
enddo

Ham = Ham/Energ_au
Ham_dir = Ham_dir/Energ_au
Ham_ex = Ham_ex/Energ_au

if ( inbox .eq. "y" ) then

if ( rdm_ori .eq. "y" ) then

call random_number(psirot)
call random_number(phirot)
call random_number(thetarot)

psirot   = psirot   * 2._dp*pi
phirot   = phirot   * 2._dp*pi
thetarot = thetarot * 2._dp*pi

rotmat(1,1)= cos(phirot)*cos(thetarot)*cos(psirot) - sin(phirot)*sin(psirot)
rotmat(1,2)=-cos(phirot)*cos(thetarot)*sin(psirot) - sin(phirot)*cos(psirot)
rotmat(1,3)= cos(phirot)*sin(thetarot)
rotmat(2,1)= sin(phirot)*cos(thetarot)*cos(psirot) + cos(phirot)*sin(psirot)
rotmat(2,2)=-sin(phirot)*cos(thetarot)*sin(psirot) + cos(phirot)*cos(psirot)
rotmat(2,3)= sin(phirot)*sin(thetarot)
rotmat(3,1)=-sin(thetarot)*cos(psirot)
rotmat(3,2)= sin(thetarot)*sin(psirot)
rotmat(3,3)= cos(thetarot)

TransHam_l(0,1,:) = matmul(rotmat,eTDM(0,1,:))*(TransDip_Ana_h1eA)
TransHam_l(0,2,:) = matmul(rotmat,eTDM(0,2,:))*(TransDip_Ana_h2eA)
TransHam_l(0,3,:) = matmul(rotmat,eTDM(0,3,:))*(TransDip_Ana_h1eB)
TransHam_l(0,4,:) = matmul(rotmat,eTDM(0,4,:))*(TransDip_Ana_h2eB)
TransHam_l(1,2,:) = matmul(rotmat,eTDM(1,2,:))*(TransDip_Ana_h1h2A)
TransHam_l(3,4,:) = matmul(rotmat,eTDM(3,4,:))*(TransDip_Ana_h1h2B)

elseif ( rdm_ori .eq. "n" ) then

TransHam_l(0,1,:) = eTDM(0,1,:)*(TransDip_Ana_h1eA)
TransHam_l(0,2,:) = eTDM(0,2,:)*(TransDip_Ana_h2eA)
TransHam_l(0,3,:) = eTDM(0,3,:)*(TransDip_Ana_h1eB)
TransHam_l(0,4,:) = eTDM(0,4,:)*(TransDip_Ana_h2eB)
TransHam_l(1,2,:) = eTDM(1,2,:)*(TransDip_Ana_h1h2A)
TransHam_l(3,4,:) = eTDM(3,4,:)*(TransDip_Ana_h1h2B)

endif

do i=0,nstates-1
do j=i+1,nstates-1
TransHam_l(j,i,:) = TransHam_l(i,j,:)
enddo
enddo

TransHam_l = TransHam_l/D_to_au

elseif ( inbox .eq. "n" ) then

TransHam(0,1) = TransDip_Ana_h1eA
TransHam(0,2) = TransDip_Ana_h2eA
TransHam(0,3) = TransDip_Ana_h1eB
TransHam(0,4) = TransDip_Ana_h2eB
TransHam(1,2) = TransDip_Ana_h1h2A
TransHam(3,4) = TransDip_Ana_h1h2B

do i=0,nstates-1
do j=i+1,nstates-1
TransHam(j,i) = TransHam(i,j)
enddo
enddo

TransHam = TransHam/D_to_au

endif

end subroutine make_Ham_he_FO

subroutine make_Ham_singl

Ham_0(1)   = minEeA(1) + minEhA(1)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aRA*1.e9_dp)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aRA*1.e9_dp)**a11_3e_ho))

Ham_0(2)   = minEeA(1) + minEhA(2) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aRA*1.e9_dp)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aRA*1.e9_dp)**a22_3e_ho))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aRA*1.e9_dp)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aRA*1.e9_dp)**a12_3e_ho))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1._dp * Ham_dir(i,j) + Ham_ex(i,j)
    endif
  enddo
enddo

do i=1,nstates-1
  do j=i+1,nstates-1
    Ham(j,i) = Ham(i,j)
    Ham_dir(j,i) = Ham_dir(i,j)
    Ham_ex(j,i) = Ham_ex(i,j)
  enddo
enddo

Ham = Ham/Energ_au
Ham_dir = Ham_dir/Energ_au
Ham_ex = Ham_ex/Energ_au

if ( inbox .eq. "y" ) then

if ( rdm_ori .eq. "y" ) then

call random_number(psirot)
call random_number(phirot)
call random_number(thetarot)

psirot   = psirot   * 2._dp*pi
phirot   = phirot   * 2._dp*pi
thetarot = thetarot * 2._dp*pi

rotmat(1,1)= cos(phirot)*cos(thetarot)*cos(psirot) - sin(phirot)*sin(psirot)
rotmat(1,2)=-cos(phirot)*cos(thetarot)*sin(psirot) - sin(phirot)*cos(psirot)
rotmat(1,3)= cos(phirot)*sin(thetarot)
rotmat(2,1)= sin(phirot)*cos(thetarot)*cos(psirot) + cos(phirot)*sin(psirot)
rotmat(2,2)=-sin(phirot)*cos(thetarot)*sin(psirot) + cos(phirot)*cos(psirot)
rotmat(2,3)= sin(phirot)*sin(thetarot)
rotmat(3,1)=-sin(thetarot)*cos(psirot)
rotmat(3,2)= sin(thetarot)*sin(psirot)
rotmat(3,3)= cos(thetarot)

TransHam_l(0,1,:) = matmul(rotmat,eTDM(0,1,:))*(TransDip_Ana_h1eA)
TransHam_l(0,2,:) = matmul(rotmat,eTDM(0,2,:))*(TransDip_Ana_h2eA)
TransHam_l(1,2,:) = matmul(rotmat,eTDM(1,2,:))*(TransDip_Ana_h1h2A)

elseif ( rdm_ori .eq. "n" ) then

TransHam_l(0,1,:) = eTDM(0,1,:)*(TransDip_Ana_h1eA)
TransHam_l(0,2,:) = eTDM(0,2,:)*(TransDip_Ana_h2eA)
TransHam_l(1,2,:) = eTDM(1,2,:)*(TransDip_Ana_h1h2A)

endif

do i=0,nstates-1
do j=i+1,nstates-1
TransHam_l(j,i,:) = TransHam_l(i,j,:)
enddo
enddo

TransHam_l = TransHam_l/D_to_au

elseif  ( inbox .eq. "n" ) then

TransHam(0,1) = TransDip_Ana_h1eA
TransHam(0,2) = TransDip_Ana_h2eA
TransHam(1,2) = TransDip_Ana_h1h2A

do i=0,nstates-1
do j=i+1,nstates-1
TransHam(j,i) = TransHam(i,j)
enddo
enddo

TransHam = TransHam/D_to_au

endif

end subroutine make_Ham_singl

!subroutine make_Ham_FS_FO
!
!Eg(1) = minEeA(1) + minEhA(1)  + V0 - Eeh1A*elec
!Eg(2) = Eg(1) + elec*(-0.523865_dp + Eg(1)/elec * 0.293665_dp)
!Eg(3) = minEeB(1) + minEhB(1)  + V0 - Eeh1A*elec
!Eg(4) = Eg(3) + elec*(-0.523865_dp + Eg(3)/elec * 0.293665_dp)
!
!TDM(1) = TransDip_Ana_h1eA
!TDM(2) = TransDip_Ana_h2eA
!TDM(3) = TransDip_Ana_h1eB
!TDM(4) = TransDip_Ana_h2eB
!
!j = 0
!k = 1
!i = 0
!
!if ( n .le. nQDA+nQDB ) then
!nbands = 2
!elseif ( n .gt. nQDA+nQDB ) then
!nbands = 4 
!endif
!
!do while ( k .le. nbands ) 
!
!Ham(j+1 ,j+1 ) = Eg(k)  - elec* (Dso1/3._dp + 3._dp*Kas)
!Ham(j+2 ,j+2 ) = Eg(k)  - elec* (Dso1/3._dp + 3._dp*Kcs)
!Ham(j+3 ,j+3 ) = Eg(k)  - elec* Kas
!Ham(j+4 ,j+4 ) = Eg(k)  - elec* 3._dp*Kas
!Ham(j+5 ,j+5 ) = Eg(k)  - elec* (3._dp*Kbs - Dxf)
!Ham(j+6 ,j+6 ) = Eg(k)  - elec* (3._dp*Kbs - Dxf)
!Ham(j+7 ,j+7 ) = Eg(k)  - elec* Kcs
!Ham(j+8 ,j+8 ) = Eg(k)  - elec* 3._dp*Kcs
!Ham(j+9 ,j+9 ) = Eg(k)  - elec* (Dso1/(-3._dp) + 3._dp*Kas)
!Ham(j+10,j+10) = Eg(k)  - elec* (Kbs - Dxf)
!Ham(j+11,j+11) = Eg(k)  - elec* (3._dp*Kbs - Dxf)
!Ham(j+12,j+12) = Eg(k)  - elec* (Dso1/(-3._dp) + 3._dp*Kcs)
!
!Ham(j+3 ,j+4 ) = elec*(-1._dp*Dso1/3._dp)
!Ham(j+4 ,j+3 ) = Ham(j+3,j+4)
!Ham(j+3 ,j+5 ) = Ham(j+3,j+4)
!Ham(j+5 ,j+3 ) = Ham(j+3,j+4)
!Ham(j+6 ,j+7 ) = Ham(j+3,j+4)
!Ham(j+7 ,j+6 ) = Ham(j+3,j+4)
!Ham(j+6 ,j+8 ) = Ham(j+3,j+4)
!Ham(j+8 ,j+6 ) = Ham(j+3,j+4)
!Ham(j+9 ,j+10) = Ham(j+3,j+4)
!Ham(j+10,j+9 ) = Ham(j+3,j+4)
!Ham(j+9 ,j+11) = Ham(j+3,j+4)
!Ham(j+11,j+9 ) = Ham(j+3,j+4)
!Ham(j+10,j+12) = Ham(j+3,j+4)
!Ham(j+12,j+10) = Ham(j+3,j+4)
!
!Ham(j+4 ,j+5 ) = elec*(Dso1/3._dp)
!Ham(j+5 ,j+4 ) = Ham(j+4,j+5)
!Ham(j+7 ,j+8 ) = Ham(j+4,j+5)
!Ham(j+8 ,j+7 ) = Ham(j+4,j+5)
!Ham(j+11,j+12) = Ham(j+4,j+5)
!Ham(j+12,j+11) = Ham(j+4,j+5)
!
!TransHam(0 ,j+3 ) = abs(TDM(k))/sqrt(3._dp) 
!TransHam(0 ,j+7 ) = abs(TDM(k))/sqrt(3._dp)
!TransHam(0 ,j+10) = abs(TDM(k))/sqrt(3._dp)
!
!k = k + 1
!j = j + 12
!
!enddo
!
!k = 1
!i = 0
!
!do while ( k .le. nbands ) 
!
!Ham(j+1 ,j+1 ) = 2*Eg(k)  - elec*(2._dp*Kbs + 2._dp*Kcs)
!Ham(j+2 ,j+2 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2._dp*Kcs - Dxf + Dso1/3._dp + Kpp)
!Ham(j+3 ,j+3 ) = 2*Eg(k)  - elec*(2._dp*Kas + Kbs + Kcs - Dxf + Dso1/3._dp + Kpp)
!Ham(j+4 ,j+4 ) = 2*Eg(k)  - elec*(2._dp*Kas + 2._dp*Kbs)
!Ham(j+5 ,j+5 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2._dp*Kcs - Dxf - Kpp)
!Ham(j+6 ,j+6 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2._dp*Kcs - Dxf + Kpp)
!Ham(j+7 ,j+7 ) = 2*Eg(k)  - elec*(Kas + 2._dp*Kbs + Kcs + Kpp)
!Ham(j+8 ,j+8 ) = 2*Eg(k)  - elec*(Kas + 2._dp*Kbs + Kcs + Kpp)
!Ham(j+9 ,j+9 ) = 2*Eg(k)  - elec*(2._dp*Kas + Kbs + Kcs - Dxf - Kpp)
!Ham(j+10,j+10) = 2*Eg(k)  - elec*(2._dp*Kas + Kbs + Kcs - Dxf + Kpp)
!Ham(j+11,j+11) = 2*Eg(k)  - elec*(Kas + Kbs + 2._dp*Kcs - Dxf - Dso1/3._dp + Kpp)
!Ham(j+12,j+12) = 2*Eg(k)  - elec*(2._dp*Kas + 2._dp*Kcs - 2._dp*Dxf)
!Ham(j+13,j+13) = 2*Eg(k)  - elec*(Kas + 2._dp*Kbs + Kcs - Kpp)
!Ham(j+14,j+14) = 2*Eg(k)  - elec*(Kas + 2._dp*Kbs + Kcs + Kpp)
!Ham(j+15,j+15) = 2*Eg(k)  - elec*(2._dp*Kas + Kbs + Kcs - Dxf - Dso1/3._dp + Kpp)
!
!Ham(j+1 ,j+4 ) = elec*Kpp
!Ham(j+4 ,j+1 ) = Ham(j+1,j+4)
!
!Ham(j+1 ,j+12 ) = elec*Kpp
!Ham(j+12 ,j+1 ) = Ham(j+1,j+12)
!
!Ham(j+4 ,j+12 ) = elec*Kpp
!Ham(j+12 ,j+4 ) = Ham(j+4,j+12)
!
!Ham(j+1 ,j+2 ) = elec*sqrt(2._dp)*Dso1/3._dp
!Ham(j+2 ,j+1 ) = Ham(j+1,j+2)
!Ham(j+3 ,j+4 ) = Ham(j+1,j+2)
!Ham(j+4 ,j+3 ) = Ham(j+1,j+2)
!Ham(j+11,j+12) = Ham(j+1,j+2)
!Ham(j+12,j+11) = Ham(j+1,j+2)
!Ham(j+12,j+15) = Ham(j+1,j+2)
!Ham(j+15,j+12) = Ham(j+1,j+2)
!
!Ham(j+5 ,j+6 ) = elec*Dso1/3._dp
!Ham(j+6 ,j+5 ) = Ham(j+5,j+6)
!Ham(j+6 ,j+7 ) = Ham(j+5,j+6)
!Ham(j+7 ,j+6 ) = Ham(j+5,j+6)
!Ham(j+8 ,j+10) = Ham(j+5,j+6)
!Ham(j+10,j+8 ) = Ham(j+5,j+6)
!Ham(j+9 ,j+10) = Ham(j+5,j+6)
!Ham(j+10,j+9 ) = Ham(j+5,j+6)
!Ham(j+11,j+13) = Ham(j+5,j+6)
!Ham(j+13,j+11) = Ham(j+5,j+6)
!Ham(j+11,j+14) = Ham(j+5,j+6)
!Ham(j+14,j+11) = Ham(j+5,j+6)
!Ham(j+15,j+13) = Ham(j+5,j+6)
!Ham(j+13,j+15) = Ham(j+5,j+6)
!Ham(j+14,j+15) = Ham(j+5,j+6)
!Ham(j+15,j+14) = Ham(j+5,j+6)
!
!Ham(j+5 ,j+7 ) = elec*Dso1/(-3._dp)
!Ham(j+7 ,j+5 ) = Ham(j+5 ,j+7 )
!Ham(j+8 ,j+9 ) = Ham(j+5 ,j+7 )
!Ham(j+9 ,j+8 ) = Ham(j+5 ,j+7 )
!
!Ham(j+13,j+14) = elec*2._dp*Dso1/3._dp
!Ham(j+14,j+13) = Ham(j+13,j+14)
!
!TransHam(0 ,j+1 ) = abs(TDM(k))/sqrt(6._dp)
!TransHam(0 ,j+4 ) = abs(TDM(k))/sqrt(6._dp)
!TransHam(0 ,j+5 ) = abs(TDM(k))/sqrt(6._dp)
!TransHam(0 ,j+9 ) = abs(TDM(k))/sqrt(6._dp)
!TransHam(0 ,j+12) = abs(TDM(k))/sqrt(6._dp)
!TransHam(0 ,j+13) = abs(TDM(k))/sqrt(6._dp)
!
!k = k + 1
!j = j + 15
!
!enddo
!
!Ham = Ham/Energ_au
!TransHam = TransHam/D_to_au
!
!do i=0,nstates-1
!do j=i+1,nstates-1
!TransHam(j,i) = TransHam(i,j)
!enddo
!enddo
!
!end subroutine make_Ham_FS_FO

subroutine make_Ham_FS_FO_P

Eg(1) = minEeA(1) + minEhA(1)  + V0 - Eeh1A*elec
Eg(2) = Eg(1) + elec*(-0.523865d0 + Eg(1)/elec * 0.293665d0)
Eg(3) = minEeB(1) + minEhB(1)  + V0 - Eeh1A*elec
Eg(4) = Eg(3) + elec*(-0.523865d0 + Eg(3)/elec * 0.293665d0)

TDM(1) = TransDip_Ana_h1eA
TDM(2) = TransDip_Ana_h2eA
TDM(3) = TransDip_Ana_h1eB
TDM(4) = TransDip_Ana_h2eB

j = 0
k = 1
i = 0

Ham(j+1 ,j+1 ) = Eg(k)  - elec* (2._dp*Dsop/3._dp + 3._dp*Kp)
Ham(j+2 ,j+2 ) = Eg(k)  - elec* (2._dp*Dsop/3._dp + 3._dp*Kp)
Ham(j+3 ,j+3 ) = Eg(k)  - elec* Kpp
Ham(j+4 ,j+4 ) = Eg(k)  - elec* 3._dp*Kpp
Ham(j+5 ,j+5 ) = Eg(k)  - elec* (3._dp*Kpp - Dsop/3._dp)
Ham(j+6 ,j+6 ) = Eg(k)  - elec* (3._dp*Kpp - Dsop/3._dp)
Ham(j+7 ,j+7 ) = Eg(k)  - elec* (3._dp*Kpp - Dsop/3._dp)
Ham(j+8 ,j+8 ) = Eg(k)  - elec* (3._dp*Kpp - Dsop/3._dp)
Ham(j+9 ,j+9 ) = Eg(k)  - elec* Kpp
Ham(j+10,j+10) = Eg(k)  - elec* 3._dp*Kpp
Ham(j+11,j+11) = Eg(k)  - elec* (3._dp*Kpp - 2._dp*Dsop/3._dp)
Ham(j+12,j+12) = Eg(k)  - elec* Kpp
Ham(j+13,j+13) = Eg(k)  - elec* 3._dp*Kpp
Ham(j+14,j+14) = Eg(k)  - elec* Kpp
Ham(j+15,j+15) = Eg(k)  - elec* 3._dp*Kpp
Ham(j+16,j+16) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+17,j+17) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+18,j+18) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+19,j+19) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+20,j+20) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+21,j+21) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+22,j+22) = Eg(k)  - elec* Kcs
Ham(j+23,j+23) = Eg(k)  - elec* 3._dp*Kpp
Ham(j+24,j+24) = Eg(k)  - elec* Kpp
Ham(j+25,j+25) = Eg(k)  - elec* 3._dp*Kpp
Ham(j+26,j+26) = Eg(k)  - elec* (3._dp*Kpp - 2*Dsop/3._dp)
Ham(j+27,j+27) = Eg(k)  - elec* (3._dp*Kpp - Dsop/3._dp)
Ham(j+28,j+28) = Eg(k)  - elec* (3._dp*Kpp - Dsop/3._dp)
Ham(j+29,j+29) = Eg(k)  - elec* (Kpp - Dxf)
Ham(j+30,j+30) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+31,j+31) = Eg(k)  - elec* (Kpp - Dxf)
Ham(j+32,j+32) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+33,j+33) = Eg(k)  - elec* (Kpp - Dxf)
Ham(j+34,j+34) = Eg(k)  - elec* (3._dp*Kpp - Dxf)
Ham(j+35,j+35) = Eg(k)  - elec* (3._dp*Kpp - Dsop/3._dp)
Ham(j+36,j+36) = Eg(k)  - elec* (3._dp*Kpp - Dsop/3._dp)

Ham(j+3 ,j+4 ) = elec*(-1._dp*Dso1/3._dp)
Ham(j+4 ,j+3 ) = Ham(j+3,j+4)
Ham(j+3 ,j+5 ) = Ham(j+3,j+4)
Ham(j+5 ,j+3 ) = Ham(j+3,j+4)
Ham(j+6 ,j+7 ) = Ham(j+3,j+4)
Ham(j+7 ,j+6 ) = Ham(j+3,j+4)
Ham(j+6 ,j+8 ) = Ham(j+3,j+4)
Ham(j+8 ,j+6 ) = Ham(j+3,j+4)
Ham(j+9 ,j+10) = Ham(j+3,j+4)
Ham(j+10,j+9 ) = Ham(j+3,j+4)
Ham(j+9 ,j+11) = Ham(j+3,j+4)
Ham(j+11,j+9 ) = Ham(j+3,j+4)
Ham(j+10,j+12) = Ham(j+3,j+4)
Ham(j+12,j+10) = Ham(j+3,j+4)

Ham(j+4 ,j+5 ) = elec*(Dso1/3._dp)
Ham(j+5 ,j+4 ) = Ham(j+4,j+5)
Ham(j+7 ,j+8 ) = Ham(j+4,j+5)
Ham(j+8 ,j+7 ) = Ham(j+4,j+5)
Ham(j+11,j+12) = Ham(j+4,j+5)
Ham(j+12,j+11) = Ham(j+4,j+5)

TransHam(0 ,j+3 ) = abs(TDM(k))/sqrt(3._dp) 
TransHam(0 ,j+7 ) = abs(TDM(k))/sqrt(3._dp)
TransHam(0 ,j+10) = abs(TDM(k))/sqrt(3._dp)

Ham = Ham/Energ_au

do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

end subroutine make_Ham_FS_FO_P

subroutine make_eTDM
eTDM(0,1,:) = vector(1._dp) 
eTDM(1,0,:) = eTDM(0,1,:) 
eTDM(0,2,:) = vector(1._dp)
eTDM(2,0,:) = eTDM(0,2,:) 
eTDM(1,2,:) = vector(1._dp)
eTDM(1,2,:) = eTDM(1,2,:) 
!eTDM(0,3,:) = vector(1._dp)
!eTDM(0,4,:) = vector(1._dp)
!eTDM(0,5,:) = vector(1._dp)
!eTDM(0,6,:) = vector(1._dp)
!eTDM(0,7,:) = vector(1._dp)
!eTDM(0,8,:) = vector(1._dp)
!eTDM(1,2,:) = vector(1._dp)
!eTDM(1,5,:) = vector(1._dp)
!eTDM(1,7,:) = vector(1._dp)
!eTDM(1,8,:) = vector(1._dp)
!eTDM(2,6,:) = vector(1._dp)
!eTDM(2,7,:) = vector(1._dp)
!eTDM(2,8,:) = vector(1._dp)
!eTDM(3,4,:) = vector(1._dp)
!eTDM(3,5,:) = vector(1._dp)
!eTDM(3,6,:) = vector(1._dp)
!eTDM(3,7,:) = vector(1._dp)
!eTDM(4,5,:) = vector(1._dp)
!eTDM(4,6,:) = vector(1._dp)
!eTDM(4,8,:) = vector(1._dp)
!eTDM(5,6,:) = vector(1._dp)
!eTDM(7,8,:) = vector(1._dp)

end subroutine

subroutine get_phases
l1( 1 )=   1._dp ; l2( 1 )=  0._dp ; l3( 1 )=  0._dp 
l1( 2 )=  -1._dp ; l2( 2 )=  0._dp ; l3( 2 )=  0._dp 
l1( 3 )=   0._dp ; l2( 3 )=  1._dp ; l3( 3 )=  0._dp 
l1( 4 )=   0._dp ; l2( 4 )= -1._dp ; l3( 4 )=  0._dp 
l1( 5 )=   0._dp ; l2( 5 )=  0._dp ; l3( 5 )=  1._dp 
l1( 6 )=   0._dp ; l2( 6 )=  0._dp ; l3( 6 )= -1._dp 
l1( 7 )=   3._dp ; l2( 7 )=  0._dp ; l3( 7 )=  0._dp 
l1( 8 )=  -3._dp ; l2( 8 )=  0._dp ; l3( 8 )=  0._dp 
l1( 9 )=   0._dp ; l2( 9 )=  3._dp ; l3( 9 )=  0._dp 
l1(10 )=   0._dp ; l2(10 )= -3._dp ; l3(10 )=  0._dp 
l1(11 )=   0._dp ; l2(11 )=  0._dp ; l3(11 )=  3._dp 
l1(12 )=   0._dp ; l2(12 )=  0._dp ; l3(12 )= -3._dp 
l1(13 )=   2._dp ; l2(13 )=  1._dp ; l3(13 )=  0._dp 
l1(14 )=  -2._dp ; l2(14 )= -1._dp ; l3(14 )=  0._dp 
l1(15 )=   2._dp ; l2(15 )= -1._dp ; l3(15 )=  0._dp 
l1(16 )=  -2._dp ; l2(16 )=  1._dp ; l3(16 )=  0._dp 
l1(17 )=   2._dp ; l2(17 )=  0._dp ; l3(17 )=  1._dp 
l1(18 )=  -2._dp ; l2(18 )=  0._dp ; l3(18 )= -1._dp 
l1(19 )=   2._dp ; l2(19 )=  0._dp ; l3(19 )= -1._dp 
l1(20 )=  -2._dp ; l2(20 )=  0._dp ; l3(20 )=  1._dp 
l1(21 )=   1._dp ; l2(21 )=  2._dp ; l3(21 )=  0._dp 
l1(22 )=  -1._dp ; l2(22 )= -2._dp ; l3(22 )=  0._dp 
l1(23 )=   1._dp ; l2(23 )= -2._dp ; l3(23 )=  0._dp 
l1(24 )=  -1._dp ; l2(24 )=  2._dp ; l3(24 )=  0._dp 
l1(25 )=   0._dp ; l2(25 )=  2._dp ; l3(25 )=  1._dp 
l1(26 )=   0._dp ; l2(26 )= -2._dp ; l3(26 )= -1._dp 
l1(27 )=   0._dp ; l2(27 )=  2._dp ; l3(27 )= -1._dp 
l1(28 )=   0._dp ; l2(28 )= -2._dp ; l3(28 )=  1._dp 
l1(29 )=   0._dp ; l2(29 )=  1._dp ; l3(29 )=  2._dp 
l1(30 )=   0._dp ; l2(30 )= -1._dp ; l3(30 )= -2._dp 
l1(31 )=   0._dp ; l2(31 )=  1._dp ; l3(31 )= -2._dp 
l1(32 )=   0._dp ; l2(32 )= -1._dp ; l3(32 )=  2._dp 
l1(33 )=   1._dp ; l2(33 )=  0._dp ; l3(33 )=  2._dp 
l1(34 )=  -1._dp ; l2(34 )=  0._dp ; l3(34 )= -2._dp 
l1(35 )=   1._dp ; l2(35 )=  0._dp ; l3(35 )= -2._dp 
l1(36 )=  -1._dp ; l2(36 )=  0._dp ; l3(36 )=  2._dp 
l1(37 )=   1._dp ; l2(37 )=  1._dp ; l3(37 )= -1._dp 
l1(38 )=  -1._dp ; l2(38 )= -1._dp ; l3(38 )=  1._dp 
l1(39 )=   1._dp ; l2(39 )= -1._dp ; l3(39 )=  1._dp 
l1(40 )=  -1._dp ; l2(40 )=  1._dp ; l3(40 )= -1._dp 
l1(41 )=  -1._dp ; l2(41 )=  1._dp ; l3(41 )=  1._dp 
l1(42 )=   1._dp ; l2(42 )= -1._dp ; l3(42 )= -1._dp 
l1(43 )=   1._dp ; l2(43 )=  1._dp ; l3(43 )=  1._dp 
l1(44 )=  -1._dp ; l2(44 )= -1._dp ; l3(44 )= -1._dp 
l1(45 )=   5._dp ; l2(45 )=  0._dp ; l3(45 )=  0._dp 
l1(46 )=  -5._dp ; l2(46 )=  0._dp ; l3(46 )=  0._dp
l1(47 )=   0._dp ; l2(47 )=  5._dp ; l3(47 )=  0._dp
l1(48 )=   0._dp ; l2(48 )= -5._dp ; l3(48 )=  0._dp
l1(49 )=   0._dp ; l2(49 )=  0._dp ; l3(49 )=  5._dp
l1(50 )=   0._dp ; l2(50 )=  0._dp ; l3(50 )= -5._dp
l1(51 )=   3._dp ; l2(51 )=  2._dp ; l3(51 )=  0._dp
l1(52 )=  -3._dp ; l2(52 )= -2._dp ; l3(52 )=  0._dp
l1(53 )=   3._dp ; l2(53 )= -2._dp ; l3(53 )=  0._dp
l1(54 )=  -3._dp ; l2(54 )=  2._dp ; l3(54 )=  0._dp
l1(55 )=   2._dp ; l2(55 )=  3._dp ; l3(55 )=  0._dp
l1(56 )=  -2._dp ; l2(56 )= -3._dp ; l3(56 )=  0._dp
l1(57 )=   2._dp ; l2(57 )= -3._dp ; l3(57 )=  0._dp
l1(58 )=  -2._dp ; l2(58 )=  3._dp ; l3(58 )=  0._dp
l1(59 )=   0._dp ; l2(59 )=  3._dp ; l3(59 )=  2._dp
l1(60 )=   0._dp ; l2(60 )= -3._dp ; l3(60 )= -2._dp
l1(61 )=   0._dp ; l2(61 )=  3._dp ; l3(61 )= -2._dp
l1(62 )=   0._dp ; l2(62 )= -3._dp ; l3(62 )=  2._dp
l1(63 )=   0._dp ; l2(63 )=  2._dp ; l3(63 )=  3._dp
l1(64 )=   0._dp ; l2(64 )= -2._dp ; l3(64 )= -3._dp
l1(65 )=   0._dp ; l2(65 )=  2._dp ; l3(65 )= -3._dp
l1(66 )=   0._dp ; l2(66 )= -2._dp ; l3(66 )=  3._dp
l1(67 )=   3._dp ; l2(67 )=  0._dp ; l3(67 )=  2._dp
l1(68 )=  -3._dp ; l2(68 )=  0._dp ; l3(68 )= -2._dp
l1(69 )=   3._dp ; l2(69 )=  0._dp ; l3(69 )= -2._dp
l1(70 )=  -3._dp ; l2(70 )=  0._dp ; l3(70 )=  2._dp
l1(71 )=   2._dp ; l2(71 )=  0._dp ; l3(71 )=  3._dp
l1(72 )=  -2._dp ; l2(72 )=  0._dp ; l3(72 )= -3._dp
l1(73 )=   2._dp ; l2(73 )=  0._dp ; l3(73 )= -3._dp
l1(74 )=  -2._dp ; l2(74 )=  0._dp ; l3(74 )=  3._dp
l1(75 )=   3._dp ; l2(75 )=  1._dp ; l3(75 )=  1._dp
l1(76 )=  -3._dp ; l2(76 )= -1._dp ; l3(76 )= -1._dp
l1(77 )=   3._dp ; l2(77 )=  1._dp ; l3(77 )= -1._dp
l1(78 )=  -3._dp ; l2(78 )= -1._dp ; l3(78 )=  1._dp
l1(79 )=   3._dp ; l2(79 )= -1._dp ; l3(79 )=  1._dp 
l1(80 )=  -3._dp ; l2(80 )=  1._dp ; l3(80 )= -1._dp
l1(81 )=   3._dp ; l2(81 )= -1._dp ; l3(81 )= -1._dp
l1(82 )=  -3._dp ; l2(82 )=  1._dp ; l3(82 )=  1._dp
l1(83 )=   1._dp ; l2(83 )=  3._dp ; l3(83 )=  1._dp
l1(84 )=  -1._dp ; l2(84 )= -3._dp ; l3(84 )= -1._dp
l1(85 )=   1._dp ; l2(85 )=  3._dp ; l3(85 )= -1._dp
l1(86 )=  -1._dp ; l2(86 )= -3._dp ; l3(86 )=  1._dp
l1(87 )=  -1._dp ; l2(87 )=  3._dp ; l3(87 )=  1._dp
l1(88 )=   1._dp ; l2(88 )= -3._dp ; l3(88 )= -1._dp
l1(89 )=   1._dp ; l2(89 )= -3._dp ; l3(89 )=  1._dp
l1(90 )=  -1._dp ; l2(90 )=  3._dp ; l3(90 )= -1._dp
l1(91 )=   1._dp ; l2(91 )=  1._dp ; l3(91 )=  3._dp
l1(92 )=  -1._dp ; l2(92 )= -1._dp ; l3(92 )= -3._dp
l1(93 )=   1._dp ; l2(93 )= -1._dp ; l3(93 )=  3._dp
l1(94 )=  -1._dp ; l2(94 )=  1._dp ; l3(94 )= -3._dp
l1(95 )=  -1._dp ; l2(95 )=  1._dp ; l3(95 )=  3._dp
l1(96 )=   1._dp ; l2(96 )= -1._dp ; l3(96 )= -3._dp
l1(97 )=   1._dp ; l2(97 )=  1._dp ; l3(97 )= -3._dp
l1(98 )=  -1._dp ; l2(98 )= -1._dp ; l3(98 )=  3._dp
l1(99 )=   2._dp ; l2(99 )=  2._dp ; l3(99 )=  1._dp
l1(100)=  -2._dp ; l2(100)= -2._dp ; l3(100)= -1._dp
l1(101)=  -2._dp ; l2(101)=  2._dp ; l3(101)=  1._dp
l1(102)=   2._dp ; l2(102)= -2._dp ; l3(102)= -1._dp
l1(103)=   2._dp ; l2(103)= -2._dp ; l3(103)=  1._dp
l1(104)=  -2._dp ; l2(104)=  2._dp ; l3(104)= -1._dp
l1(105)=   2._dp ; l2(105)=  2._dp ; l3(105)= -1._dp
l1(106)=  -2._dp ; l2(106)= -2._dp ; l3(106)=  1._dp
l1(107)=   2._dp ; l2(107)=  1._dp ; l3(107)=  2._dp
l1(108)=  -2._dp ; l2(108)= -1._dp ; l3(108)= -2._dp
l1(109)=  -2._dp ; l2(109)=  1._dp ; l3(109)=  2._dp
l1(110)=   2._dp ; l2(110)= -1._dp ; l3(110)= -2._dp
l1(111)=   2._dp ; l2(111)=  1._dp ; l3(111)= -2._dp
l1(112)=  -2._dp ; l2(112)= -1._dp ; l3(112)=  2._dp
l1(113)=   2._dp ; l2(113)= -1._dp ; l3(113)=  2._dp
l1(114)=  -2._dp ; l2(114)=  1._dp ; l3(114)= -2._dp
l1(115)=   1._dp ; l2(115)=  2._dp ; l3(115)=  2._dp
l1(116)=  -1._dp ; l2(116)= -2._dp ; l3(116)= -2._dp
l1(117)=   1._dp ; l2(117)= -2._dp ; l3(117)=  2._dp
l1(118)=  -1._dp ; l2(118)=  2._dp ; l3(118)= -2._dp
l1(119)=   1._dp ; l2(119)=  2._dp ; l3(119)= -2._dp
l1(120)=  -1._dp ; l2(120)= -2._dp ; l3(120)=  2._dp
l1(121)=  -1._dp ; l2(121)=  2._dp ; l3(121)=  2._dp
l1(122)=   1._dp ; l2(122)= -2._dp ; l3(122)= -2._dp
l1(123)=   4._dp ; l2(123)=  1._dp ; l3(123)=  0._dp
l1(124)=  -4._dp ; l2(124)= -1._dp ; l3(124)=  0._dp 
l1(125)=   4._dp ; l2(125)= -1._dp ; l3(125)=  0._dp
l1(126)=  -4._dp ; l2(126)=  1._dp ; l3(126)=  0._dp
l1(127)=   1._dp ; l2(127)=  4._dp ; l3(127)=  0._dp
l1(128)=  -1._dp ; l2(128)= -4._dp ; l3(128)=  0._dp
l1(129)=  -1._dp ; l2(129)=  4._dp ; l3(129)=  0._dp
l1(130)=   1._dp ; l2(130)= -4._dp ; l3(130)=  0._dp
l1(131)=   1._dp ; l2(131)=  0._dp ; l3(131)=  4._dp
l1(132)=  -1._dp ; l2(132)=  0._dp ; l3(132)= -4._dp
l1(133)=  -1._dp ; l2(133)=  0._dp ; l3(133)=  4._dp
l1(134)=   1._dp ; l2(134)=  0._dp ; l3(134)= -4._dp
l1(135)=   4._dp ; l2(135)=  0._dp ; l3(135)=  1._dp
l1(136)=  -4._dp ; l2(136)=  0._dp ; l3(136)= -1._dp
l1(137)=   4._dp ; l2(137)=  0._dp ; l3(137)= -1._dp
l1(138)=  -4._dp ; l2(138)=  0._dp ; l3(138)=  1._dp
l1(139)=   0._dp ; l2(139)=  4._dp ; l3(139)=  1._dp
l1(140)=   0._dp ; l2(140)= -4._dp ; l3(140)= -1._dp
l1(141)=   0._dp ; l2(141)=  4._dp ; l3(141)= -1._dp
l1(142)=   0._dp ; l2(142)= -4._dp ; l3(142)=  1._dp
l1(143)=   0._dp ; l2(143)=  1._dp ; l3(143)=  4._dp
l1(144)=   0._dp ; l2(144)= -1._dp ; l3(144)= -4._dp
l1(145)=   0._dp ; l2(145)= -1._dp ; l3(145)=  4._dp
l1(146)=   0._dp ; l2(146)=  1._dp ; l3(146)= -4._dp
end subroutine

subroutine make_Ham_FS_FO

if ( idlink .eq. 20 ) then
include 'Parameters-dir-02.f90'
include 'Parameters-ex-02.f90'
elseif ( idlink .eq. 55 ) then
include 'Parameters-dir-55.f90'
include 'Parameters-ex-55.f90'
endif

Ham_0(1)     = minEeA(1) + minEhA(1)  + V0
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aRA*1e9_dp)**a11_3d_ho))
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aRA*1e9_dp)**a11_3e_ho))

Ham_0(2)     = minEeA(1) + minEhA(2) + V0
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aRA*1e9_dp)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aRA*1e9_dp)**a22_3e_ho))

if ( n .gt. nQDA+nQDB ) then
Ham_0(3)     = minEeB(1) + minEhB(1)  + V0
Ham_dir(3,3) = elec*(a11_1d_ho + a11_2d_ho / ((aRB*1e9_dp)**a11_3d_ho))
Ham_ex(3,3)  = elec*(a11_1e_ho + a11_2e_ho / ((aRB*1e9_dp)**a11_3e_ho))

Ham_0(4)     = minEeB(1) + minEhB(2) + V0
Ham_dir(4,4) = elec*(a22_1d_ho + a22_2d_ho / ((aRB*1e9_dp)**a22_3d_ho))
Ham_ex(4,4)  = elec*(a22_1e_ho + a22_2e_ho / ((aRB*1e9_dp)**a22_3e_ho))
endif

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aRA*1e9_dp)**a12_3d_ho))
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aRA*1e9_dp)**a12_3e_ho))

if ( n .gt. nQDA+nQDB ) then
Ham_dir(3,4) = elec*(a12_1d_ho + a12_2d_ho / ((aRB*1e9_dp)**a12_3d_ho))
Ham_ex(3,4)  = elec*(a12_1e_ho + a12_2e_ho / ((aRB*1e9_dp)**a12_3e_ho))

Ham_dir(1,3) = elec*(a13_1d_he + a13_2d_he / ((aRA*1e9_dp)**a13_3d_he * (aRB*1e9_dp)**a13_4d_he))
Ham_ex(1,3)  = elec*(a13_1e_he + a13_2e_he / ((aRA*1e9_dp)**a13_3e_he * (aRB*1e9_dp)**a13_4e_he))

Ham_dir(1,4) = elec*(a14_1d_he + a14_2d_he / ((aRA*1e9_dp)**a14_3d_he * (aRB*1e9_dp)**a14_4d_he))
Ham_ex(1,4)  = elec*(a14_1e_he + a14_2e_he / ((aRA*1e9_dp)**a14_3e_he * (aRB*1e9_dp)**a14_4e_he))

Ham_dir(2,3) = elec*(a14_1d_he + a14_2d_he / ((aRB*1e9_dp)**a14_3d_he * (aRA*1e9_dp)**a14_4d_he))
Ham_ex(2,3)  = elec*(a14_1e_he + a14_2e_he / ((aRB*1e9_dp)**a14_3e_he * (aRA*1e9_dp)**a14_4e_he))

Ham_dir(2,4) = elec*(a24_1d_he + a24_2d_he / ((aRA*1e9_dp)**a24_3d_he * (aRB*1e9_dp)**a24_4d_he))
Ham_ex(2,4)  = elec*(a24_1e_he + a24_2e_he / ((aRA*1e9_dp)**a24_3e_he * (aRB*1e9_dp)**a24_4e_he))
endif

Eg(1) = Ham_0(1) 
Eg(2) = Ham_0(2) 
if ( n .gt. nQDA+nQDB ) then
Eg(3) = Ham_0(3)
Eg(4) = Ham_0(4)
endif

TDM(1) = TransDip_Ana_h1eA
TDM(2) = TransDip_Ana_h2eA
if ( n .gt. nQDA+nQDB ) then
TDM(3) = TransDip_Ana_h1eB
TDM(4) = TransDip_Ana_h2eB
endif

if ( n .le. nQDA+nQDB ) then
nbands = 2
elseif ( n .gt. nQDA+nQDB ) then
nbands = 4 
endif

j=0
k=1

do while ( k .le. nbands ) 

Ham(j+1 ,j+1 ) = Eg(k) + Ham_ex(k,k)  - elec* Dso1/3.d0 
Ham(j+2 ,j+2 ) = Eg(k) + Ham_ex(k,k)  - elec* Dso1/3.d0 
Ham(j+3 ,j+3 ) = Eg(k) - Ham_dir(k,k) + Ham_ex(k,k)
Ham(j+4 ,j+4 ) = Eg(k) + Ham_ex(k,k) 
Ham(j+5 ,j+5 ) = Eg(k) + Ham_ex(k,k)  + elec* Dxf
Ham(j+6 ,j+6 ) = Eg(k) + Ham_ex(k,k)  + elec* Dxf
Ham(j+7 ,j+7 ) = Eg(k) - Ham_dir(k,k) + Ham_ex(k,k) 
Ham(j+8 ,j+8 ) = Eg(k) + Ham_ex(k,k) 
Ham(j+9 ,j+9 ) = Eg(k) + Ham_ex(k,k)  + elec* Dso1/3.d0
Ham(j+10,j+10) = Eg(k) - Ham_dir(k,k) + Ham_ex(k,k) + elec* Dxf
Ham(j+11,j+11) = Eg(k) + Ham_ex(k,k)  + elec* Dxf
Ham(j+12,j+12) = Eg(k) + Ham_ex(k,k)  + elec* Dso1/3.d0

Ham(j+3 ,j+4 ) = elec*(-1.d0*Dso1/3.d0) 
Ham(j+3 ,j+5 ) = elec*(-1.d0*Dso1/3.d0) 
Ham(j+6 ,j+7 ) = elec*(-1.d0*Dso1/3.d0) 
Ham(j+6 ,j+8 ) = elec*(-1.d0*Dso1/3.d0) 
Ham(j+9 ,j+10) = elec*(-1.d0*Dso1/3.d0) 
Ham(j+9 ,j+11) = elec*(-1.d0*Dso1/3.d0) 
Ham(j+10,j+12) = elec*(-1.d0*Dso1/3.d0)

Ham(j+4 ,j+5 ) = elec*(Dso1/3.d0) 
Ham(j+7 ,j+8 ) = elec*(Dso1/3.d0) 
Ham(j+11,j+12) = elec*(Dso1/3.d0) 

TransHam(0 ,j+3 ) = abs(TDM(k)) 
TransHam(0 ,j+7 ) = abs(TDM(k))
TransHam(0 ,j+10) = abs(TDM(k))

j = j + 12
k = k + 1

enddo

if ( Biex .eq. 'y' ) then

do k = 1,nbands
   do l = k,nbands

!Ham(j+1 ,j+1 ) = 2*Eg(k)  - elec*(2._dp*Kbs + 2._dp*Kcs)
!Ham(j+2 ,j+2 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2._dp*Kcs - Dxf + Dso1/3._dp + Kpp)
!Ham(j+3 ,j+3 ) = 2*Eg(k)  - elec*(2._dp*Kas + Kbs + Kcs - Dxf + Dso1/3._dp + Kpp)
!Ham(j+4 ,j+4 ) = 2*Eg(k)  - elec*(2._dp*Kas + 2._dp*Kbs)
!Ham(j+5 ,j+5 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2._dp*Kcs - Dxf - Kpp)
!Ham(j+6 ,j+6 ) = 2*Eg(k)  - elec*(Kas + Kbs + 2._dp*Kcs - Dxf + Kpp)
!Ham(j+7 ,j+7 ) = 2*Eg(k)  - elec*(Kas + 2._dp*Kbs + Kcs + Kpp)
!Ham(j+8 ,j+8 ) = 2*Eg(k)  - elec*(Kas + 2._dp*Kbs + Kcs + Kpp)
!Ham(j+9 ,j+9 ) = 2*Eg(k)  - elec*(2._dp*Kas + Kbs + Kcs - Dxf - Kpp)
!Ham(j+10,j+10) = 2*Eg(k)  - elec*(2._dp*Kas + Kbs + Kcs - Dxf + Kpp)
!Ham(j+11,j+11) = 2*Eg(k)  - elec*(Kas + Kbs + 2._dp*Kcs - Dxf - Dso1/3._dp + Kpp)
!Ham(j+12,j+12) = 2*Eg(k)  - elec*(2._dp*Kas + 2._dp*Kcs - 2._dp*Dxf)
!Ham(j+13,j+13) = 2*Eg(k)  - elec*(Kas + 2._dp*Kbs + Kcs - Kpp)
!Ham(j+14,j+14) = 2*Eg(k)  - elec*(Kas + 2._dp*Kbs + Kcs + Kpp)
!Ham(j+15,j+15) = 2*Eg(k)  - elec*(2._dp*Kas + Kbs + Kcs - Dxf - Dso1/3._dp + Kpp)
Ham(j+1 ,j+1 ) = Eg(k) + Eg(l) + 4.d0 * ( - Ham_dir(k,l) + Ham_ex(k,l) )  
Ham(j+2 ,j+2 ) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec*(- Dxf + Dso1/3._dp + Kpp)
Ham(j+3 ,j+3 ) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec*(- Dxf + Dso1/3._dp + Kpp)
Ham(j+4 ,j+4 ) = Eg(k) + Eg(l) + 4.d0 * ( - Ham_dir(k,l) + Ham_ex(k,l) ) 
Ham(j+5 ,j+5 ) = Eg(k) + Eg(l) + 4.d0 * ( - Ham_dir(k,l) + Ham_ex(k,l) ) - elec*(- Dxf - Kpp)
Ham(j+6 ,j+6 ) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec*(- Dxf + Kpp)
Ham(j+7 ,j+7 ) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec* Kpp
Ham(j+8 ,j+8 ) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec* Kpp
Ham(j+9 ,j+9 ) = Eg(k) + Eg(l) + 4.d0 * ( - Ham_dir(k,l) + Ham_ex(k,l) ) - elec*( - Dxf - Kpp)
Ham(j+10,j+10) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec*( - Dxf + Kpp)
Ham(j+11,j+11) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec*( - Dxf - Dso1/3._dp + Kpp)
Ham(j+12,j+12) = Eg(k) + Eg(l) + 4.d0 * ( - Ham_dir(k,l) + Ham_ex(k,l) ) - elec*( - 2._dp*Dxf)
Ham(j+13,j+13) = Eg(k) + Eg(l) + 4.d0 * ( - Ham_dir(k,l) + Ham_ex(k,l) ) - elec*( - Kpp)
Ham(j+14,j+14) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec*(   Kpp)
Ham(j+15,j+15) = Eg(k) + Eg(l) + 4.d0 * (   Ham_ex(k,l)                ) - elec*( - Dxf - Dso1/3._dp + Kpp)

Ham(j+1 ,j+4 ) = elec*Kpp
Ham(j+4 ,j+1 ) = Ham(j+1,j+4)

Ham(j+1 ,j+12 ) = elec*Kpp
Ham(j+12 ,j+1 ) = Ham(j+1,j+12)

Ham(j+4 ,j+12 ) = elec*Kpp
Ham(j+12 ,j+4 ) = Ham(j+4,j+12)

Ham(j+1 ,j+2 ) = elec*sqrt(2._dp)*Dso1/3._dp
Ham(j+2 ,j+1 ) = Ham(j+1,j+2)
Ham(j+3 ,j+4 ) = Ham(j+1,j+2)
Ham(j+4 ,j+3 ) = Ham(j+1,j+2)
Ham(j+11,j+12) = Ham(j+1,j+2)
Ham(j+12,j+11) = Ham(j+1,j+2)
Ham(j+12,j+15) = Ham(j+1,j+2)
Ham(j+15,j+12) = Ham(j+1,j+2)

Ham(j+5 ,j+6 ) = elec*Dso1/3._dp
Ham(j+6 ,j+5 ) = Ham(j+5,j+6)
Ham(j+6 ,j+7 ) = Ham(j+5,j+6)
Ham(j+7 ,j+6 ) = Ham(j+5,j+6)
Ham(j+8 ,j+10) = Ham(j+5,j+6)
Ham(j+10,j+8 ) = Ham(j+5,j+6)
Ham(j+9 ,j+10) = Ham(j+5,j+6)
Ham(j+10,j+9 ) = Ham(j+5,j+6)
Ham(j+11,j+13) = Ham(j+5,j+6)
Ham(j+13,j+11) = Ham(j+5,j+6)
Ham(j+11,j+14) = Ham(j+5,j+6)
Ham(j+14,j+11) = Ham(j+5,j+6)
Ham(j+15,j+13) = Ham(j+5,j+6)
Ham(j+13,j+15) = Ham(j+5,j+6)
Ham(j+14,j+15) = Ham(j+5,j+6)
Ham(j+15,j+14) = Ham(j+5,j+6)

Ham(j+5 ,j+7 ) = elec*Dso1/(-3._dp)
Ham(j+7 ,j+5 ) = Ham(j+5 ,j+7 )
Ham(j+8 ,j+9 ) = Ham(j+5 ,j+7 )
Ham(j+9 ,j+8 ) = Ham(j+5 ,j+7 )

Ham(j+13,j+14) = elec*2._dp*Dso1/3._dp
Ham(j+14,j+13) = Ham(j+13,j+14)

!TransHam(0 ,j+1 ) = abs(TDM(k))
!TransHam(0 ,j+4 ) = abs(TDM(k))
!TransHam(0 ,j+5 ) = abs(TDM(k))
!TransHam(0 ,j+9 ) = abs(TDM(k))
!TransHam(0 ,j+12) = abs(TDM(k))
!TransHam(0 ,j+13) = abs(TDM(k))

j = j + 15

   enddo
enddo

endif

!H12A
Ham(1 ,13) = Ham_ex(1,2)
Ham(2 ,14) = Ham_ex(1,2)
Ham(3 ,15) = - Ham_dir(1,2) + Ham_ex(1,2)
Ham(4 ,16) = Ham_ex(1,2)
Ham(5 ,17) = Ham_ex(1,2)
Ham(6 ,18) = Ham_ex(1,2)
Ham(7 ,19) = - Ham_dir(1,2) + Ham_ex(1,2)
Ham(8 ,20) = Ham_ex(1,2)
Ham(9 ,21) = Ham_ex(1,2)
Ham(10,22) = - Ham_dir(1,2) + Ham_ex(1,2)
Ham(11,23) = Ham_ex(1,2)
Ham(12,24) = Ham_ex(1,2)

TransHam(1 ,13) = abs(TransDip_Ana_h1h2A)
TransHam(2 ,14) = abs(TransDip_Ana_h1h2A)
TransHam(3 ,15) = abs(TransDip_Ana_h1h2A)
TransHam(4 ,16) = abs(TransDip_Ana_h1h2A)
TransHam(5 ,17) = abs(TransDip_Ana_h1h2A)
TransHam(6 ,18) = abs(TransDip_Ana_h1h2A)
TransHam(7 ,19) = abs(TransDip_Ana_h1h2A)
TransHam(8 ,20) = abs(TransDip_Ana_h1h2A)
TransHam(9 ,21) = abs(TransDip_Ana_h1h2A)
TransHam(10,22) = abs(TransDip_Ana_h1h2A)
TransHam(11,23) = abs(TransDip_Ana_h1h2A)
TransHam(12,24) = abs(TransDip_Ana_h1h2A)

if ( n .gt. nQDA+nQDB ) then
!H12B
Ham(25,37) = Ham_ex(3,4)
Ham(26,38) = Ham_ex(3,4)
Ham(27,39) = - Ham_dir(3,4) + Ham_ex(3,4)
Ham(28,40) = Ham_ex(3,4)
Ham(29,41) = Ham_ex(3,4)
Ham(30,42) = Ham_ex(3,4)
Ham(31,43) = - Ham_dir(3,4) + Ham_ex(3,4)
Ham(32,44) = Ham_ex(3,4)
Ham(33,45) = Ham_ex(3,4)
Ham(34,46) = - Ham_dir(3,4) + Ham_ex(3,4)
Ham(35,47) = Ham_ex(3,4)
Ham(36,48) = Ham_ex(3,4)

TransHam(25,37) = abs(TransDip_Ana_h1h2B)
TransHam(26,38) = abs(TransDip_Ana_h1h2B)
TransHam(27,39) = abs(TransDip_Ana_h1h2B)
TransHam(28,40) = abs(TransDip_Ana_h1h2B)
TransHam(29,41) = abs(TransDip_Ana_h1h2B)
TransHam(30,42) = abs(TransDip_Ana_h1h2B)
TransHam(31,43) = abs(TransDip_Ana_h1h2B)
TransHam(32,44) = abs(TransDip_Ana_h1h2B)
TransHam(33,45) = abs(TransDip_Ana_h1h2B)
TransHam(34,46) = abs(TransDip_Ana_h1h2B)
TransHam(35,47) = abs(TransDip_Ana_h1h2B)
TransHam(36,48) = abs(TransDip_Ana_h1h2B)
endif

if ( Von .eq. 'y' ) then 

!H1A1B
Ham(1 ,25) = Ham_ex(1,3)
Ham(2 ,26) = Ham_ex(1,3)
Ham(3 ,27) = - Ham_dir(1,3) + Ham_ex(1,3)
Ham(4 ,28) = Ham_ex(1,3)
Ham(5 ,29) = Ham_ex(1,3)
Ham(6 ,30) = Ham_ex(1,3)
Ham(7 ,31) = - Ham_dir(1,3) + Ham_ex(1,3)
Ham(8 ,32) = Ham_ex(1,3)
Ham(9 ,33) = Ham_ex(1,3)
Ham(10,34) = - Ham_dir(1,3) + Ham_ex(1,3)
Ham(11,35) = Ham_ex(1,3)
Ham(12,36) = Ham_ex(1,3)

!H1A2B
Ham(1 ,37) = Ham_ex(1,4)
Ham(2 ,38) = Ham_ex(1,4)
Ham(3 ,39) = - Ham_dir(1,4) + Ham_ex(1,4)
Ham(4 ,40) = Ham_ex(1,4)
Ham(5 ,41) = Ham_ex(1,4)
Ham(6 ,42) = Ham_ex(1,4)
Ham(7 ,43) = - Ham_dir(1,4) + Ham_ex(1,4)
Ham(8 ,44) = Ham_ex(1,4)
Ham(9 ,45) = Ham_ex(1,4)
Ham(10,46) = - Ham_dir(1,4) + Ham_ex(1,4)
Ham(11,47) = Ham_ex(1,4)
Ham(12,48) = Ham_ex(1,4)

!H2A1B
Ham(13,25) = Ham_ex(2,3)
Ham(14,26) = Ham_ex(2,3)
Ham(15,27) = - Ham_dir(2,3) + Ham_ex(2,3)
Ham(16,28) = Ham_ex(2,3)
Ham(17,29) = Ham_ex(2,3)
Ham(18,30) = Ham_ex(2,3)
Ham(19,31) = - Ham_dir(2,3) + Ham_ex(2,3)
Ham(20,32) = Ham_ex(2,3)
Ham(21,33) = Ham_ex(2,3)
Ham(22,34) = - Ham_dir(2,3) + Ham_ex(2,3)
Ham(23,35) = Ham_ex(2,3)
Ham(24,36) = Ham_ex(2,3)

!H2A2B
Ham(13,37) = Ham_ex(2,4)
Ham(14,38) = Ham_ex(2,4)
Ham(15,39) = - Ham_dir(2,4) + Ham_ex(2,4)
Ham(16,40) = Ham_ex(2,4)
Ham(17,41) = Ham_ex(2,4)
Ham(18,42) = Ham_ex(2,4)
Ham(19,43) = - Ham_dir(2,4) + Ham_ex(2,4)
Ham(20,44) = Ham_ex(2,4)
Ham(21,45) = Ham_ex(2,4)
Ham(22,46) = - Ham_dir(2,4) + Ham_ex(2,4)
Ham(23,47) = Ham_ex(2,4)
Ham(24,48) = Ham_ex(2,4)

endif

do i=0,nstates-1
do j=i+1,nstates-1
Ham(j,i) = Ham(i,j)
TransHam(j,i) = TransHam(i,j)
enddo
enddo

Ham = Ham/Energ_au
TransHam = TransHam/D_to_au

end subroutine make_Ham_FS_FO

end module Make_Ham
