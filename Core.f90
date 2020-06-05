Ham_ei = Ham
allocate(lambda(0:nstates-1),source = 0._dp)
allocate(iwork2(3+5*nstates),source=0)
allocate(work1(6*nstates),work2(1+6*nstates+2*nstates*nstates),source=0._dp)
lworku=1+6*nstates+2*nstates*nstates
liworku=3+5*nstates
ierr=0
call dsyevd('v','u',nstates, Ham_ei(0:nstates-1,0:nstates-1),nstates,lambda,work2,lworku,iwork2,liworku,ierr)
deallocate(work1)
deallocate(work2)
deallocate(iwork2)

do i=0,nstates-1
Ham_l(i,i) = lambda(i)
enddo

!!!Make eigenstate TDM
if ( rdm_ori .eq. "n" ) then
Mat(:,:) = matmul(TransHam(:,:),Ham_ei(:,:))
TransHam_ei(:,:) = matmul(transpose(Ham_ei(:,:)),Mat(:,:))
elseif ( rdm_ori .eq. "y" ) then
Matx(:,:) = matmul(TransHam_l(:,:,1),Ham_ei(:,:))
Maty(:,:) = matmul(TransHam_l(:,:,2),Ham_ei(:,:))
Matz(:,:) = matmul(TransHam_l(:,:,3),Ham_ei(:,:))
TransHam_ei_l(:,:,1) = matmul(transpose(Ham_ei(:,:)),Matx(:,:))
TransHam_ei_l(:,:,2) = matmul(transpose(Ham_ei(:,:)),Maty(:,:))
TransHam_ei_l(:,:,3) = matmul(transpose(Ham_ei(:,:)),Matz(:,:))
TransHam_ei = sqrt(TransHam_ei_l(:,:,1)**2 + TransHam_ei_l(:,:,2)**2 + TransHam_ei_l(:,:,1)**2)
endif

!!!!!labelling of eigenstates
do i=0,nstates-1
maxid(i) = maxval(abs(Ham_ei(:,i)))
enddo
do i = 0,nstates-1
do j = 0,nstates-1
if (abs(Ham_ei(i,j)) .eq. maxid(j)) then
        zero(j) = i
        exit
endif
enddo
enddo

if ( Dyn_L .eq. 'y' ) then

write(form1,'("(i4,i4,i4,i4,1x,100(e24.16,1x))")')

sigD = 0._dp

if ( Dec_L .eq. "y" ) then

sig1=0.00533526_dp
sig2=0.00576011_dp
sig3=0.0131895_dp
sig4=0.0137016_dp
rho12= 0.563941_dp
rho13=0.946482_dp
rho14=0.272191_dp
rho23=0.304623_dp
rho24=0.925269_dp
rho34=0._dp

sigDiag(1 ) = sig1
sigDiag(2 ) = sig2
sigDiag(3 ) = sig3
sigDiag(4 ) = sig4
sigDiag(5 ) = sig1
sigDiag(7 ) = sqrt(sig1**2 + sig2**2 - 2._dp*rho12*sig1*sig2)
sigDiag(8 ) = sqrt(sig1**2 + sig3**2 - 2._dp*rho13*sig1*sig3)
sigDiag(9 ) = sqrt(sig1**2 + sig4**2 - 2._dp*rho14*sig1*sig4)
sigDiag(10) = sig2
sigDiag(11) = sqrt(sig1**2 + sig2**2 - 2._dp*rho12*sig1*sig2)
sigDiag(13) = sqrt(sig2**2 + sig3**2 - 2._dp*rho23*sig2*sig3)
sigDiag(14) = sqrt(sig2**2 + sig4**2 - 2._dp*rho24*sig2*sig4)
sigDiag(15) = sig3
sigDiag(16) = sqrt(sig1**2 + sig3**2 - 2._dp*rho13*sig1*sig3)
sigDiag(17) = sqrt(sig2**2 + sig3**2 - 2._dp*rho23*sig2*sig3)
sigDiag(19) = sqrt(sig3**2 + sig4**2 - 2._dp*rho34*sig3*sig4)
sigDiag(20) = sig4
sigDiag(21) = sqrt(sig1**2 + sig4**2 - 2._dp*rho14*sig1*sig4)
sigDiag(22) = sqrt(sig2**2 + sig4**2 - 2._dp*rho24*sig2*sig4)
sigDiag(23) = sqrt(sig3**2 + sig4**2 - 2._dp*rho34*sig3*sig4)

do i=1,nstates2-1
do j=i,nstates2-1
sigD(i ,j ) = 1._dp/((2.d0*pi*6.582119570e-16_dp/sqrt(sigDiag(i)*sigDiag(j))/t_au))
sigD(i ,j ) = sigD(j ,i )
enddo
enddo

sigD = (2._dp*pi*sigD)**2

endif

do i=0,nstates2-1
do j=0,nstates2-1
merge_diag(i,j)  = real(merge(1,0,i.eq.j),kind=dp)
merge_odiag(i,j) = real(merge(0,1,i.eq.j),kind=dp)
enddo
enddo

do i=0,nstates-1
do j=i+1,nstates-1
haml(i,j)=TransHam_ei(i,j)
haml(j,i)=haml(i,j)
enddo
enddo

do i=0,nstates-1
haml(i,i)=lambda(i)
enddo

!!!!!Liouvillian
do k=0,nstates-1
do l=0,nstates-1
do i=0,nstates-1
!!!!!Commutator [H,Eij]
xliou(k,l,i,l) = xliou(k,l,i,l) + dcmplx(haml(i,k),0._dp)
enddo
do j=0,nstates-1
xliou(k,l,k,j) = xliou(k,l,k,j) - dcmplx(haml(l,j),0._dp)
enddo
enddo
enddo

kl = 0 ; kc = 0

!!!!Renumber xLiou
kl = -1
do i=0,nstates-1
do j=0,nstates-1
kl = kl + 1
irow(kl,1) = i
irow(kl,2) = j
kc = -1
do k=0,nstates-1
do l=0,nstates-1
kc = kc + 1
icol(kc,1) = k
icol(kc,2) = l
lfield(kl,kc) = xliou(i,j,k,l)
enddo
enddo
enddo
enddo

!!!!!Print xLiouvillian
do i=0,nstates-1
do j=0,nstates-1
do k=0,nstates-1
do l=0,nstates-1
write(Liou_f,form1) i,j,k,l,xliou(i,j,k,l)
enddo
enddo
enddo
enddo

endif


!! file writing
write(form_mat,'("(",i0,"ES16.5E2)")') nstates
write(form_TDM,'("(",i0,"f16.8)")') nstates
if (n .le. nQDA+nQDB) then
write(form_arr,'("(f16.8,16x,",i0,"f14.8)")') nstates+1
elseif (n .gt. nQDA+nQDB) then
write(form_arr,'("(",i0,"f16.12)")') nstates+3
endif
write(form_abs,'("(2f16.8,2x,i0)")')
write(form_pop,'("(",i0,"ES17.8E3,ES25.16E3)")') nstates+1
write(form_com,'("(ES12.5E3,",i0,"ES18.8E3)")') nstates+1
write(form_com_L,'("(ES12.5E3,",i0,"ES18.8E3,ES28.15)")') nstates2+1
write(form_pop_L,'("(ES12.5E3,",i0,"ES18.8E3,f28.15)")') nstates
write(form_DipSpec,'("(ES12.5E3,",i0,"ES18.8E3)")') nstates+2

if ( noMat .eq. "n" ) then
do i=1,size(matrices)
if (n .le. nQDA+nQDB) then
write(matrices(i),'(i0,f16.8)') n , aRA*1.d9
elseif (n .gt. nQDA+nQDB) then
write(matrices(i),'(i0,2f16.8)') n , aRA*1.d9, aRB*1.d9
endif
enddo
do i=0,nstates-1
write(H_0_f    ,form_mat) (Ham(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_dir_f  ,form_mat) (Ham_dir(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_ex_f   ,form_mat) (Ham_ex(i,j)*Energ_au/elec, j=0,nstates-1)
write(H_JK_f   ,form_mat) ((-1.d0*Ham_dir(i,j) + Ham_ex(i,j))*Energ_au/elec, j=0,nstates-1)
write(Tmat_0_f ,form_TDM) (TransHam(i,j), j=0,nstates-1)
write(H_ei_f   ,form_mat) (Ham_ei(i,j), j=0,nstates-1)
write(Tmat_ei_f,form_TDM) (TransHam_ei(i,j), j=0,nstates-1)
if ( inbox .eq. "y" ) then
write(Tmat_x_f,form_TDM) (TransHam_ei_l(i,j,1), j=0,nstates-1)
write(Tmat_y_f,form_TDM) (TransHam_ei_l(i,j,2), j=0,nstates-1)
write(Tmat_z_f,form_TDM) (TransHam_ei_l(i,j,3), j=0,nstates-1)
endif
enddo
do i=1,size(matrices)
write(matrices(i),*)
enddo
endif

do i=0,nstates-1
write(Abs_imp_f,form_abs) lambda(i)*Energ_au/elec, (TransHam_ei(0,i))**2 ,i
enddo

write(label_0,*) n, (zero(j), j=0,nstates-1)

if (n .le. nQDA+nQDB) then
write(Etr_0_f,form_arr) aRA*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)
write(Etr_ei_f,form_arr) aRA*1.d9, (lambda(i)*Energ_au/elec, i=0,nstates-1)
elseif (n .gt. nQDA+nQDB) then
write(Etr_0_f,form_arr) aRA*1.d9, aRB*1.d9, (Ham(i,i)*Energ_au/elec, i=0,nstates-1)
write(Etr_ei_f,form_arr) aRA*1.d9, aRB*1.d9, (lambda(i)*Energ_au/elec, i=0,nstates-1)
endif

!Proceed with dynamics
if ( ( Dyn_0 .eq. 'y' ) .or. ( Dyn_ei .eq. 'y' ) .or. ( Dyn_L .eq. 'y' ) ) then
call RK_0_ei
endif

deallocate(TransHam,TransHam_ei_l,TransHam_l,TransHam_d,TransHam_ei,Mat,Matx,Maty,Matz,Ham,Ham_l,Ham_0,Ham_dir,Ham_ex,Ham_ei,haml)
deallocate(Transvec,TransMat_ei,lambda,xc,k1,k2,k3,k4,k5,k6,k7,k8,c0,xc_ei,xc_L,xc0,pop,sigDiag)
deallocate(k1_L,k2_L,k3_L,k4_L,k5_L,k6_L,k7_L,k8_L,xliou,lfield)
deallocate(merge_diag,merge_odiag,icol,irow,maxid,zero,sigD)
