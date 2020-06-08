module Variables_au

use omp_lib
use Constants_au
use Normal
use Vectors

implicit none

   character*5  :: vers
   character*6  :: pgeom
   character*64 :: popc, hmti, norm, tdmM, hmt0, outputdir, norm_ei, popc_ei, Re_c_ei, Im_c_ei, integ,model, line, dummy
   character*64 :: Re_c, Im_c, syst_n, Re_c_l, Im_c_L, cov, cov2, Pop_c_L
   character*1  :: o_Norm, o_Over, o_Coul, o_DipS, o_Osci, o_Exti, o_DipD, dyn, hamilt, get_ei, finest, get_sp
   character*1  :: Dyn_0, Dyn_ei, inbox, Dyn_L,doFT,CEP1,CEP2,CEP3,nofiles, doCovar,doAbs
   character*1  :: rdm_ori, noMat, Dec_L, Von
   logical :: isit
   integer,dimension(10) :: matrices
   integer :: Pulse_f,Tmat_0_f,Tmat_ei_f,Tmat_x_f,Tmat_y_f,Tmat_z_f,H_0_f,H_dir_f,H_ex_f,H_JK_f,TransAbs, DipSpec_conv_f
   integer :: popc_0_f,popc_ei_f,norm_0_f,norm_ei_f,Re_c_ei_f,Im_c_ei_f,Re_c_0_f,Im_c_0_f,TDip_ei_f,tmp, nbands,Liou_f, getMat
   integer :: Re_c_L_f,Im_c_L_f,H_ei_f,Etr_0_f,Etr_ei_f,Abs_imp_f, t2, DipSpec, pol, npol, P_Match_f, DipSpec_R_f,DipSpec_NR_f
   character*64 :: form_mat,form_arr,form_abs,form_pop,form_com,form_TDM,form1,form_com_L, form_DipSpec, form_pop_L
   integer :: syst, ndots, n, rmin, rmax, nsys, npulses, nstates, ntime,i,j,k,l,t,lwork, info, idlink, threads, nQDA, nQDB
   integer :: nhomoA,nhomoB,nhetero,totsys,ndim,nQD, EminID,EmaxID,Estep,nmax,io,abso, kc, kl, nstates2, DipSpec_s, pol2
   integer :: TransAbs_NR, TransAbs_R, TransAbs_s, nFT, FTpow, t_ana, f_ana, Pop_c_L_f, TransAbs_P, ierr, liworku, lworku
   integer :: scos_ssin,TransAbs_7_f,TransAbs_17_f,TransAbs_33_f,TransAbs_39_f,TransAbs_41_f,TransAbs_44_f,DipSpec_7_f,DipSpec_17_f
   integer :: DipSpec_33_f,DipSpec_39_f,DipSpec_41_f,DipSpec_44_f, sphere, DipSpec_43_f, TransAbs_43_f
   integer :: Tmat_avg_f, Tmat_avgx_f, Tmat_avgy_f, Tmat_avgz_f, Etr_avg_f, popc_avg_f, Re_c_avg_f, Im_c_avg_f, TransAbs_avg
   integer :: DipSpec_avg, P_Match_avg, TransAbs_avg_7_f, TransAbs_avg_17_f, TransAbs_avg_33_f, TransAbs_avg_39_f, TransAbs_avg_41_f
   integer :: TransAbs_avg_43_f, TransAbs_avg_44_f, DipSpec_avg_7_f, DipSpec_avg_17_f, DipSpec_avg_33_f, DipSpec_avg_39_f
   integer :: DipSpec_avg_41_f, DipSpec_avg_43_f, DipSpec_avg_44_f, error, label_0, sigfit, Dmat, rhofit, nt, nt3
   integer,allocatable :: seed(:),icol(:,:),irow(:,:),iwork2(:),zero(:), values(:)
   real(dp) :: a13_1d_he,a13_2d_he,a13_3d_he,a13_4d_he,a15_1d_he,a15_2d_he,a15_3d_he,a15_4d_he,a17_1d_he,a17_2d_he,a17_3d_he,&
               a17_4d_he,a24_1d_he,a24_2d_he,a24_3d_he,a24_4d_he,a26_1d_he,a26_2d_he,a26_3d_he,a26_4d_he,a28_1d_he,a28_2d_he,&
               a28_3d_he,a28_4d_he,a35_1d_he,a35_2d_he,a35_3d_he,a35_4d_he,a37_1d_he,a37_2d_he,a37_3d_he,a37_4d_he,a46_1d_he,&
               a46_2d_he,a46_3d_he,a46_4d_he,a48_1d_he,a48_2d_he,a48_3d_he,a48_4d_he,a55_1d_he,a55_2d_he,a55_3d_he,a55_4d_he,&
               a57_1d_he,a57_2d_he,a57_3d_he,a57_4d_he,a66_1d_he,a66_2d_he,a66_3d_he,a66_4d_he,a68_1d_he,a68_2d_he,a68_3d_he,&
               a68_4d_he,a77_1d_he,a77_2d_he,a77_3d_he,a77_4d_he,a88_1d_he,a88_2d_he,a88_3d_he,a88_4d_he,a14_1d_he,a14_2d_he,&
               a14_3d_he,a14_4d_he,a16_1d_he,a16_2d_he,a16_3d_he,a16_4d_he,a18_1d_he,a18_2d_he,a18_3d_he,a18_4d_he,a23_1d_he,&
               a23_2d_he,a23_3d_he,a23_4d_he,a25_1d_he,a25_2d_he,a25_3d_he,a25_4d_he,a27_1d_he,a27_2d_he,a27_3d_he,a27_4d_he,&
               a36_1d_he,a36_2d_he,a36_3d_he,a36_4d_he,a38_1d_he,a38_2d_he,a38_3d_he,a38_4d_he,a45_1d_he,a45_2d_he,a45_3d_he,&
               a45_4d_he,a47_1d_he,a47_2d_he,a47_3d_he,a47_4d_he,a56_1d_he,a56_2d_he,a56_3d_he,a56_4d_he,a58_1d_he,a58_2d_he,&
               a58_3d_he,a58_4d_he,a67_1d_he,a67_2d_he,a67_3d_he,a67_4d_he,a78_1d_he,a78_2d_he,a78_3d_he,a78_4d_he
   real(dp) :: a13_1e_he,a13_2e_he,a13_3e_he,a13_4e_he,a15_1e_he,a15_2e_he,a15_3e_he,a15_4e_he,a17_1e_he,a17_2e_he,a17_3e_he,&
               a17_4e_he,a24_1e_he,a24_2e_he,a24_3e_he,a24_4e_he,a26_1e_he,a26_2e_he,a26_3e_he,a26_4e_he,a28_1e_he,a28_2e_he,&
               a28_3e_he,a28_4e_he,a35_1e_he,a35_2e_he,a35_3e_he,a35_4e_he,a37_1e_he,a37_2e_he,a37_3e_he,a37_4e_he,a46_1e_he,&
               a46_2e_he,a46_3e_he,a46_4e_he,a48_1e_he,a48_2e_he,a48_3e_he,a48_4e_he,a55_1e_he,a55_2e_he,a55_3e_he,a55_4e_he,&
               a57_1e_he,a57_2e_he,a57_3e_he,a57_4e_he,a66_1e_he,a66_2e_he,a66_3e_he,a66_4e_he,a68_1e_he,a68_2e_he,a68_3e_he,&
               a68_4e_he,a77_1e_he,a77_2e_he,a77_3e_he,a77_4e_he,a88_1e_he,a88_2e_he,a88_3e_he,a88_4e_he,a14_1e_he,a14_2e_he,&
               a14_3e_he,a14_4e_he,a16_1e_he,a16_2e_he,a16_3e_he,a16_4e_he,a18_1e_he,a18_2e_he,a18_3e_he,a18_4e_he,a23_1e_he,&
               a23_2e_he,a23_3e_he,a23_4e_he,a25_1e_he,a25_2e_he,a25_3e_he,a25_4e_he,a27_1e_he,a27_2e_he,a27_3e_he,a27_4e_he,&
               a36_1e_he,a36_2e_he,a36_3e_he,a36_4e_he,a38_1e_he,a38_2e_he,a38_3e_he,a38_4e_he,a45_1e_he,a45_2e_he,a45_3e_he,&
               a45_4e_he,a47_1e_he,a47_2e_he,a47_3e_he,a47_4e_he,a56_1e_he,a56_2e_he,a56_3e_he,a56_4e_he,a58_1e_he,a58_2e_he,&
               a58_3e_he,a58_4e_he,a67_1e_he,a67_2e_he,a67_3e_he,a67_4e_he,a78_1e_he,a78_2e_he,a78_3e_he,a78_4e_he
   real(dp) :: tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8, time_ana, aR_avgA, aR_avgB, start, finish, debug,slope
   real(dp) :: tpx1,tpx2,tpx3,tpx4,tpx5,tpx6,tpx7,tpx8,sig1,sig2,sig3,sig4,rho12,rho13,rho14,rho23,rho24,rho34
   real(dp) :: tpy1,tpy2,tpy3,tpy4,tpy5,tpy6,tpy7,tpy8
   real(dp) :: tpz1,tpz2,tpz3,tpz4,tpz5,tpz6,tpz7,tpz8, re1, im1, re2, im2, re3, im3, re4, im4
   real(dp) :: aA, aB, me, mh, eps, epsout, V0, omegaLO, V0eV, minr, maxr, rsteps, side, link
   real(dp) :: sigma_conv,Emin,Emax,x,w1,w2,w,wstep,time_FFT, FTscale, psirot, phirot, thetarot
   real(dp) :: vertex, zbase, alphae, alphah1, alphah2, betae, betah1, betah2, rdm1, rdm2
   real(dp) :: dispQD, displink, rdmlinker, rdmQDA, rdmQDB, t01, t02, t03, timestep, totaltime, distQD, Kpp, Dsop, Kp
   real(dp) :: omega01, omega02, omega03, phase01, phase02, phase03, width01, width02, width03, Ed01, Ed02, Ed03
   real(dp) :: pulse1, pulse2, pulse3, test, time, cnorm, cnormabs, cnormconj, cnorm2, cnorm2_ei, Kas, Kbs, Kcs, Dso1, Dso2, Dxf
   real(dp) :: Eeh1A,Eeh2A,AeA,Ah1A,Ah2A,BeA,Bh1A,Bh2A,kineA,kinh1A,kinh2A,kouteA,kouth1A ,kouth2A,TransDip_Ana_h1eA
   real(dp) :: TransDip_Ana_h2eA,TransDip_Ana_h1h2A,aRA,aRB,V0eA,V0eB,V0hA,V0hB
   real(dp) :: Eeh1B,Eeh2B,AeB,Ah1B,Ah2B,BeB,Bh1B,Bh2B,kineB,kinh1B,kinh2B,kouteB,kouth1B,kouth2B,TransDip_Ana_h1eB
   real(dp) :: TransDip_Ana_h2eB,TransDip_Ana_h1h2B
   real(dp),allocatable :: aR(:), epsin(:), epsR(:), V0e(:), V0h(:), linker(:), Eg(:), TDM(:), eTDM(:,:,:), sigma(:)
   real(dp),allocatable :: epsinA(:), epsinB(:), epsRA(:), epsRB(:), l1(:), l2(:), l3(:)
   real(dp),allocatable :: Cb_eh1(:), Cb_eh2(:), Norm_Ana_e(:), Norm_Ana_h1(:), Norm_Ana_h2(:), Eeh1(:), Eeh2(:)
   real(dp),allocatable :: OverlapAna_h1e(:), OverlapAna_h2e(:), Cb_Num_eh1(:), Cb_Num_eh1_eh2(:), Cb_Num_eh2(:)
   real(dp),allocatable :: TransDip_Num_h1e(:), TransDip_Num_h2e(:), work(:), lambda(:)
   real(dp),allocatable :: minEe(:,:),minEh(:,:), minEeA(:),minEeB(:)
   real(dp),allocatable :: minEhA(:),minEhB(:), pown(:), pulsesn(:)
   real(dp),allocatable :: TransDip_Ana_h1e(:), TransDip_Ana_h2e(:), Oscillator_Ana_h1e(:), Oscillator_Ana_h2e(:), Transvec(:)
   real(dp),allocatable :: ExctCoef_h1e(:), ExctCoef_h2e(:), Ham(:,:), E0(:), c0(:), TransHam(:,:), Hamt(:,:,:), c(:,:)
   real(dp),allocatable :: TransMat_ei(:,:), TransHam0(:,:), Ham_0(:), Ham_dir(:,:), Ham_ex(:,:), Ham_ei(:,:), Ham_l(:,:), haml(:,:)
   real(dp),allocatable :: TransDip_Ana_h1h2(:), TransHam_ei(:,:), Mat(:,:), QDcoor(:), Dcenter(:), Pe1(:), Pe2(:), Pe3(:)
   real(dp),allocatable :: TransHam_d(:,:,:), TransHam_l(:,:,:), TransHam_ei_l(:,:,:), k_1(:), k_2(:), k_3(:), work1(:), work2(:)
   real(dp),allocatable :: Matx(:,:), Maty(:,:), Matz(:,:),spec(:),dipole(:,:), Ham0_avg(:,:), sigD(:,:), sigDiag(:), powtemp(:)
   real(dp),allocatable :: pow(:),pow_gaus(:),pulses(:), pow_s(:,:), pow_gaus_s(:,:), pulses_FFT(:)
   real(dp),allocatable :: Scov(:,:), pop(:,:),merge_diag(:,:),merge_odiag(:,:), scos(:), ssin(:), lfield(:,:), rho(:,:)
   real(dp),allocatable :: TransHam_avg(:,:),TransHam_avg_l(:,:,:),lambda_avg(:), bloblo(:), blabla(:),maxid(:), rotmat(:,:)
 complex(8) :: ct1, ct2, ct3, ct4, xt01, xt02, xt03, xhbar, im, xwidth, xomega , xEd, xh, xphase, xtime, xhbar_au
 complex(8) :: integPol, integPol_diff, integPolconv, pow_41, pow_41_conv
 complex(8),allocatable :: xHam(:,:) , xHamt(:,:,:), xTransHam(:,:), xE0(:), xHamtk2(:,:,:), xHamtk3(:,:,:), xHamtk4(:,:,:)
 complex(8),allocatable :: xc0(:), xc(:,:), xc_ei(:,:), xcnew(:,:), k1(:), k2(:), k3(:) , k4(:), xHam_ei(:,:), pow_pol_conv(:)
 complex(8),allocatable :: k1_L(:), k2_L(:), k3_L(:) , k4_L(:),  k5_L(:), k6_L(:), k7_L(:) , k8_L(:), xpow_pol(:,:),xliou(:,:,:,:)
 complex(8),allocatable :: dk1(:), dk2(:), dk3(:) , dk4(:), k5(:), k6(:), k7(:) , k8(:), pow_pol(:,:), pow_pol_gaus(:,:)
 complex(8),allocatable :: xc_ei_av(:,:), xctemp(:),xlfield(:,:),xc_L(:,:), xpow_gaus(:),xpulse(:),wft(:),wftp(:),wftf(:),xhugo(:)
 complex(8),allocatable :: wft_pol(:,:),wftf_pol(:,:), pow_pol_diff(:), wft_s(:,:), wftf_s(:,:), xpow_gaus_s(:), xpulse2(:)
 complex(8),allocatable :: k1_rho(:,:), k2_rho(:,:), k3_rho(:,:), k4_rho(:,:), k5_rho(:,:), k6_rho(:,:), k7_rho(:,:), k8_rho(:,:)
 complex(8),allocatable :: xc_rho(:,:,:), wftf_t1(:), pow_poln(:,:)

contains 

subroutine getVariables

NAMELIST /outputs/   inbox,rdm_ori,get_sp,get_ei,Dyn_0,Dyn_ei,Dyn_L,Dec_L,&
                     doAbs,doFT,nofiles,noMat,doCovar
NAMELIST /elecSt/    model,Von,me,mh,eps,epsout,V0eV,omegaLO,slope,side
NAMELIST /fineStruc/ Kas,Kbs,Kcs,Kpp,Dso1,Dso2,Dxf
NAMELIST /pulses/    integ,npulses,t01,t02,t03,timestep,totaltime,omega01,omega02,omega03,phase01,phase02,phase03,&
                     width01,width02,width03,Ed01,Ed02,Ed03,CEP1,CEP2,CEP3,pgeom,vertex
NAMELIST /syst/      nQDA,nQDB,nhomoA,nhomoB,nhetero,dispQD,idlink,aA,aB     
NAMELIST /FT/        FTpow

open(150,file='QD_quest.def',form='formatted')
read(150,NML=outputs)
read(150,NML=elecSt)
read(150,NML=fineStruc)
read(150,NML=pulses)
read(150,NML=syst)
read(150,NML=FT)

im         = dcmplx(0._dp,1._dp)
me         = me*m0
mh         = mh*m0
V0         = V0eV*elec
npol       = 146
call init_random_seed()

timestep   =  timestep  * 1.e-15_dp/t_au  
!compute total length of simulation
selectcase(npulses)
case(0) ; totaltime = (totaltime)     * 1.e-15_dp/t_au
case(1) ; totaltime = (totaltime+t01) * 1.e-15_dp/t_au
case(2) ; totaltime = (totaltime+t02) * 1.e-15_dp/t_au
case(3) ; totaltime = (totaltime+t03) * 1.e-15_dp/t_au
endselect
t01        =  t01       * 1.e-15_dp/t_au       
t02        =  t02       * 1.e-15_dp/t_au       
t03        =  t03       * 1.e-15_dp/t_au       
width01    =  width01   * 1.e-15_dp/t_au   
width02    =  width02   * 1.e-15_dp/t_au   
width03    =  width03   * 1.e-15_dp/t_au   
omega01    =  omega01   * t_au*(2._dp*pi)  
omega02    =  omega02   * t_au*(2._dp*pi)  
omega03    =  omega03   * t_au*(2._dp*pi)  
Ed01       =  Ed01/E_au        
Ed02       =  Ed02/E_au        
Ed03       =  Ed03/E_au        
xh         =  dcmplx(timestep,0._dp)
ntime      =  nint(totaltime/timestep)
nt         =  nint((totaltime-t03)/timestep)
nt3        =  nint(t03/timestep)

selectcase(CEP1)
case('r') ; call random_number(phase01) ; phase01 = phase01 * 2._dp * pi
case('p') ; phase01 = pi
case('0') ; phase01 = 0._dp
endselect
selectcase(CEP2)
case('r') ; call random_number(phase02) ; phase02 = phase02 * 2._dp * pi
case('p') ; phase02 = pi
case('0') ; phase02 = 0._dp
endselect
selectcase(CEP3)
case('r') ; call random_number(phase03) ; phase03 = phase03 * 2._dp * pi
case('p') ; phase03 = pi
case('0') ; phase03 = 0._dp
endselect

selectcase(npulses)
case(3) ; pulse1 = 1._dp ; pulse2 = 1._dp ; pulse3 = 1._dp
case(2) ; pulse1 = 1._dp ; pulse2 = 1._dp ; pulse3 = 0._dp
case(1) ; pulse1 = 1._dp ; pulse2 = 0._dp ; pulse3 = 0._dp
endselect

if ( inbox .eq. 'y' ) then

allocate(k_1(3),k_2(3),k_3(3),Pe1(3),Pe2(3),Pe3(3),source=0._dp)

selectcase(pgeom)

case('boxcar')

zbase = 1000._dp*sqrt(2._dp)*sqrt(1-cos(pi*vertex/180._dp))/sqrt(cos(pi*vertex/180._dp))
k_1(1) =(-2._dp * pi / (cl/(omega01/t_au))) * ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))
k_1(2) = (2._dp * pi / (cl/(omega01/t_au))) * ( 1000._dp      ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))
k_1(3) =(-2._dp * pi / (cl/(omega01/t_au))) * ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))
k_2(1) =(-2._dp * pi / (cl/(omega02/t_au))) * ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))  
k_2(2) = (2._dp * pi / (cl/(omega02/t_au))) * ( 1000._dp      ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))
k_2(3) = (2._dp * pi / (cl/(omega02/t_au))) * ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))
k_3(1) = (2._dp * pi / (cl/(omega03/t_au))) * ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))
k_3(2) = (2._dp * pi / (cl/(omega03/t_au))) * ( 1000._dp      ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))
k_3(3) = (2._dp * pi / (cl/(omega03/t_au))) * ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2 + (zbase / 2._dp)**2))

Pe1(1) =  (-zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))
Pe1(2) =  ( 1000._dp      ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))
Pe1(3) =  (-zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))
Pe2(1) =  (-zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))  
Pe2(2) =  ( 1000._dp      ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))
Pe2(3) =  ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))
Pe3(1) =  ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))
Pe3(2) =  ( 1000._dp      ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))
Pe3(3) =  ( zbase / 2._dp ) / (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase / 2._dp)**2))

case('triang')

zbase= 1000._dp*sqrt(2._dp*(1._dp-cos(pi*vertex/180._dp))/(1._dp-(1-cos(pi*vertex/180._dp))/(1-cos(pi*120._dp/180._dp)))) 
k_1(1) = (2._dp*pi/(cl/(omega01/t_au)))*(zbase / 2._dp ) /        (sqrt((zbase / 2._dp)**2+(1000._dp)**2+&
(zbase/(2._dp*sqrt(3._dp)))**2))
k_1(2) = (2._dp*pi/(cl/(omega01/t_au)))*(1000._dp      ) /        (sqrt((zbase / 2._dp)**2+(1000._dp)**2+&
(zbase/(2._dp*sqrt(3._dp)))**2))
k_1(3) = (2._dp*pi/(cl/(omega01/t_au)))*(zbase/(2._dp*sqrt(3._dp)))/(sqrt((zbase / 2._dp)**2+(1000._dp)**2+&
(zbase/(2._dp*sqrt(3._dp)))**2))
k_2(1) =(-2._dp*pi/(cl/(omega02/t_au)))*(zbase / 2._dp ) /        (sqrt((zbase / 2._dp)**2+(1000._dp)**2+&
(zbase/(2._dp*sqrt(3._dp)))**2))
k_2(2) = (2._dp*pi/(cl/(omega02/t_au)))*(1000._dp      ) /        (sqrt((zbase / 2._dp)**2+(1000._dp)**2+&
(zbase/(2._dp*sqrt(3._dp)))**2))
k_2(3) = (2._dp*pi/(cl/(omega02/t_au)))*(zbase/(2._dp*sqrt(3._dp)))/(sqrt((zbase / 2._dp)**2+(1000._dp)**2+&
(zbase/(2._dp*sqrt(3._dp)))**2))
k_3(1) = ( 2._dp*pi/(cl/(omega03/t_au)))*0._dp 
k_3(2) = ( 2._dp*pi/(cl/(omega03/t_au)))*(1000._dp      ) /        (sqrt((1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))
k_3(3) = (-2._dp*pi/(cl/(omega03/t_au)))*(zbase/(2._dp*sqrt(3._dp)))/(sqrt((1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))

Pe1(1) =  ( zbase / 2._dp ) /        (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))
Pe1(2) =  ( 1000._dp      ) /        (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))
Pe1(3) =  ( zbase/(2._dp*sqrt(3._dp)))/(sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))
Pe2(1) =  ( zbase / 2._dp ) /        (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))
Pe2(2) =  ( 1000._dp      ) /        (sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))
Pe2(3) =  ( zbase/(2._dp*sqrt(3._dp)))/(sqrt((zbase / 2._dp)**2+(1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))
Pe3(1) =  0._dp 
Pe3(2) =  ( 1000._dp      ) /        (sqrt((1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))
Pe3(3) =  ( zbase/(2._dp*sqrt(3._dp)))/(sqrt((1000._dp)**2+(zbase/(2._dp*sqrt(3._dp)))**2))

case('straig')

k_1(1) = 0._dp 
k_1(2) = 0._dp
k_1(3) = 0._dp
k_2(1) = 0._dp
k_2(2) = 0._dp
k_2(3) = 0._dp
k_3(1) = 0._dp
k_3(2) = 0._dp
k_3(3) = 0._dp

Pe1(1) = 1._dp/sqrt(3._dp) 
Pe1(2) = 1._dp/sqrt(3._dp)
Pe1(3) = 1._dp/sqrt(3._dp)
Pe2(1) = 1._dp/sqrt(3._dp)
Pe2(2) = 1._dp/sqrt(3._dp)
Pe2(3) = 1._dp/sqrt(3._dp)
Pe3(1) = 1._dp/sqrt(3._dp)
Pe3(2) = 1._dp/sqrt(3._dp)
Pe3(3) = 1._dp/sqrt(3._dp)

endselect

endif

nQD    = nQDA+nQDB
ndim   = nhomoA+nhomoB+nhetero
nsys   = nQDA+nQDB+nhomoA+nhomoB+nhetero     
totsys = nsys+nhomoA+nhomoB+nhetero

!if ( inbox .eq. 'y' ) then
!open(newunit=sphere,file='sphere.xyz')
!write(sphere,*) nint(totsys/2._dp) 
!write(sphere,*) 
!Dcenter(:) = vectorin(1.e16_dp) 
!Dcenter(:) = (Dcenter(:) + 1e16_dp**(1._dp/3._dp)) * 1.e-9_dp
!write(sphere,*) 'H', Dcenter(1), Dcenter(2), Dcenter(3) 
!endif

end subroutine getVariables
   
end module Variables_au
