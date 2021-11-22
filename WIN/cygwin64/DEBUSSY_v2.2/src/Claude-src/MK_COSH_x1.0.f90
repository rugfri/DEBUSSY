program Paint_It_Black
use nano_deftyp
use ATOMIX
use paper_blood
use GEODIS
use fast_PLACE
use GURKEN_2012
use Sultans_Of_Swing
use INPUT_CLUMKFILE

implicit real(DP)(a-h,o-z),integer(I4B)(i-n)
character(11),parameter :: file_inp0 = 'coshmkQ.ini'
real(DP) :: DVds(3,3),DVrs(3,3),LVds(3),Diam_max_nm,MTds(3,3),cellCM(3),eigMT(3)
integer(I4B) :: NCL2MK(2),allocOK, lobo(3),upbo(3),usedlin
integer(I4B) :: Ncube1(3),Ncube2(3),mmax2(3),mmax(3),mmax_K(3),mmax_S(3)
real(DP)     :: DCL2MK(2)
real(DP) :: parapar(9),wl_tt(2)
real(DP) :: aux_latvec(3),aux_Dvec(3),xyzP(3)
integer(I4B) :: dr_step=1
character(3) :: clushape
character(512) :: pwd,rlx,inpath, file_inp,namestr, rl
character(8) :: pearsy
Logical      :: DB_exist, isordnot, orthoXY, eqatpair, zerolay, &
                isorigin, isoriginP, isoriginL, output_xyz, &
                nonzerowid, allsizes, nzw0
character(8)  :: para_model = 'WelbIsot'
logical :: occs_are_one = .true., done_once=.false.

real(DP),allocatable :: Cube1(:,:,:,:),Cube1r(:,:,:),Cube2(:,:,:,:)
integer(I4B),allocatable  :: Cube1i(:,:,:)
real(DP),allocatable      :: radat(:), mulat(:),mulats(:)
real(DP),allocatable      :: patte_len(:,:,:),patte_invmul(:,:,:)
integer(I4B) :: aux1(3),aux2(3), nlatpoi_clu(2),nacu(2)
real(DP) :: Raux1(3),Raux2(3), ooo1(2),ooo2(2),soccp(3), xlatpoi_clu(2), ooo_latV(3)
character(len=14),dimension(:,:),allocatable :: finu
real(DP),allocatable,save     :: xnat0(:,:,:),summul0(:,:,:),termcon0(:,:,:)
integer(I4B),allocatable,save :: nat0(:,:,:),ndi0(:,:,:)
real(DP),allocatable,save     :: xncelKS(:,:,:),olastV(:)
integer(I4B),allocatable,save :: ncelKS(:,:,:),ndidi0V(:),popd0V(:)
character(len=777) :: fxyz,namef_multi,path_cels
character(len=7):: wnu
character(len=14):: ccccc

call def_eps


num_grow_dir = 2

DB_exist = .true.

call CPU_TIME(ttot0)
call GET_PWD(pwd=pwd,lpwd=lpwd)

call me_Al(fn_clumkMULT=file_inp0,usedlin=usedlin,name_files=namef_multi,path_cel_files=path_cels)

call openread_clumkX(clumkfn=file_inp0, &
                     rlx=rlx,lrlx=lrlx, &
                     itra=itra,nclu=NCL2MK,dclu=DCL2MK, &
                     clusha=clushape, &
                     parapar=parapar,smoothing_width=smowid, wlenttmax=wl_tt, &
                     do_allsizes=allsizes, sample_allsteps=multiple_steps, &
                     step_base = dr_step, para_model=para_model,occ1=occs_are_one, do_the_xyz=output_xyz, &
                     nskip=usedlin)
if (clushape /= 'CSH') then
  print*,'Shape given is '//clushape//' but I deal only with shape SPH. Bye bye...'
  stop 'I deal only with shape SPH. Bye bye...'
endif

print*,'MK_COSH: after openread_clumkX: ncelltypes is ',ncelltypes
allocate(numcells_byphase(ncelltypes))

working_wavelength=wl_tt(1)
working_ttmax=wl_tt(2)
if (working_wavelength == zero) working_wavelength=wl_ourqmax

!___________________________________________ Paths and filename of .cel/.cely
file_inp=''
namef_multi=trim(adjustl(namef_multi))
lnamef_multi=len_trim(namef_multi)
print*,' &&&&&&    Built output filename root as  ',trim(namef_multi)

path_cels=trim(adjustl(path_cels))
lpath_cels=len_trim(path_cels)
if (path_cels(lpath_cels:lpath_cels)/=separator) then
  lpath_cels=lpath_cels+1
  path_cels(lpath_cels:lpath_cels) = separator
endif
file_inp(1:lnamef_multi+len_trim(path_cels))=trim(path_cels)//namef_multi(1:lnamef_multi)
namestr=trim(namef_multi)//'_'
Lnamestr=len_trim(namef_multi)+1
linpath=lpath_cels
inpath=trim(adjustl(path_cels))
print'(a)',' INPATH:   '//inpath(1:linpath)
print'(a)',' NAMESTR:  '//namestr(1:Lnamestr)
print'(a)',' FILE_INP: '//trim(file_inp)
!___________________________________________ Paths and filename of .cel/.cely DONE


!_________________________ PARACRYSTALLINE
IF (para_model == 'WelbAnys') then
  rcorra=parapar(1)
  scorra=parapar(2)
  rcorrb=parapar(3)
  scorrb=parapar(4)
  cocoxy=parapar(5)
  sigmaxa=parapar(6)
  sigmaxb=parapar(7)
  B000_c=parapar(8)
  sigma_min_sph=parapar(9)
!!!!!!!!!!!
  var_targ = (sigmaxa**2+sigmaxb**2)
  r_targ=20.d0
  slope_var = var_targ/r_targ
!!!!!!!!!!!
  variance_min = sigma_min_sph**2
  corr_coe = max(-one+s4eps_DP,min(one-s4eps_DP,cocoxy))
  sigma_c_sq = B000_c/(pi2*pi2)
  max_intrinsic_sigma = max(zero, sigmaxa,sigmaxb)
  nzw0 = (max_intrinsic_sigma > scepsDP)
else IF (para_model == 'WelbIsot') then
  max_intrinsic_sigma=parapar(1) ! limit variance is  2 * max_intrinsic_sigma**2
  decorr_len = parapar(2)

!!!!!!!!!!!
  var_targ = two*(max_intrinsic_sigma**2)
  r_targ=decorr_len
  slope_var = var_targ/r_targ
!!!!!!!!!!!
  variance_min = zero
  corr_coe = zero
  nzw0 = (max_intrinsic_sigma > scepsDP)
endif

print*, 'Working Wavelength = ',working_wavelength
print*, '2*Theta Max = ',working_ttmax
print*,'MK_COSH: before READALLO_CELY: ncelltypes is ',ncelltypes

!!!!!!!!!!!!!!!!!!!! DONE reading name root of input file and its path
call READALLO_CELY(cely_fn_V=celfinames,lcely_fn_V=len_trim(celfinames), &
                   ncl2mk=ncl2mk,dcl2mk=dcl2mk,itra=itra,geom=clushape, &
                   diace=Delta_R,cdx=cdx,c_lattice=c, occ1=occs_are_one, ncelltypes1=ncelltypes )
print*,'MK_COSH: after READALLO_CELY: ncelltypes is ',ncelltypes
print*, 'NCL2MK =', NCL2MK
print'(a,i3,a,i3,a)', 'Building clusters with 0 to ',NCL2MK(1)-1,' KERN/CORE layers PLUS 0 to ',NCL2MK(2)-1,' SCHALE/SHELL layers'
DCL2MK = (NCL2MK-1)*Delta_R*0.2d0
print'(2(a,1x,g12.4),a)', 'Building clusters with 0 to ',DCL2MK(1),' nm KERN/CORE PLUS 0 to ',DCL2MK(2),' nm SCHALE/SHELL'
print*, 'clushape = ',clushape
kcelp=0
do icel1=1,ncelltypes
  do icel2=icel1,ncelltypes
    kcelp=kcelp+1
    write(wnu,'(i3.3,"-",i3.3)') icel1,icel2
    iup=find_unit()
    print'(a)','WRITING FILE: '//trim(file_inp)//'_'//wnu//'.pat'
    open(iup,status='replace',file=trim(file_inp)//'_'//wnu//'.pat')
    write(iup,'(a,i4)')'# N. atom pairs (global) : ',napair_all
    write(iup,'(a,i4," [",2i4,"]")')'# This cell pair : ',kcelp,icel1,icel2
    write(iup,'(a,i4,a,i6," / ",11i6)')"# Cell ",icel1," : atoms: (sum / by species) ",sum(Celty(icel1)%nat),Celty(icel1)%nat
    write(iup,'(a,i4,a,i6," / ",11i6)')"# Cell ",icel2," : atoms: (sum / by species) ",sum(Celty(icel2)%nat),Celty(icel2)%nat
    do ipa=1,napair_all
      is1=glopair(1,ipa)
      is2=glopair(2,ipa)
      is1l=return_address_atsp(is1,icel1)
      is2l=return_address_atsp(is2,icel2)
      if (min(is1l,is2l)==0) cycle
      izz1=Celty(icel1)%Z_at(is1l)
      izz2=Celty(icel2)%Z_at(is2l)
      sopair = sum(patte(4,1:numpatv(ipa,kcelp),ipa,kcelp))
    !  if (is1==is2) sopair=sopair-Celty(1)%nat(is1)
      write(iup,*)is1,is2,Celty(icel1)%nat(is1l),Celty(icel2)%nat(is2l)
      write(iup,*)izz1,izz2,numpatv(ipa,kcelp),sopair
      ss1=zero; ss2=zero
      do jx=1,numpatv(ipa,kcelp)
        zzz=maxval(abs(patte(1:3,jx,ipa,kcelp)))
        fk=one
        if (icel1==icel2.and.zzz>sceps_DP) fk=two
        ss1=ss1+patte(4,jx,ipa,kcelp); ss2=ss2+patte(4,jx,ipa,kcelp)*fk
        write(iup,*) patte(1:4,jx,ipa,kcelp),patte(4,jx,ipa,kcelp)*fk
      enddo
      write(iup,*)'                SUM ',ss1,ss2
    enddo
    close(iup)
  enddo
enddo

abcabg_loc=abcabgy
sampled_folder='DISTANCES'; lsampled_folder=9
basepath_sampled=inpath(1:linpath); lbasepath_sampled=linpath
if (Ndelta1==Ndelta2) then
  print'(a,/,a,i3.3,a)','Sampled distances .smp files will be written in ',&
  inpath(1:linpath)//separator//trim(sampled_folder)//separator//'SAMPTO',nint(1000*Deltas(Ndelta1)),separator
else
  print'(a,/,a,/,a,i3.3," ... ",i3.3)','Sampled distances .smp files will be written in ',&
  inpath(1:linpath)//separator//trim(sampled_folder)//separator//'SAMPTO***'//separator, &
  ' where *** = ',nint(1000*Deltas(Ndelta1)),nint(1000*Deltas(Ndelta2))
endif

DVds=zero
ttrig=( cos(degrees_to_radians*abcabgy(4))-cos(degrees_to_radians*abcabgy(5))*cos(degrees_to_radians*abcabgy(6)) ) / &
           sin(degrees_to_radians*abcabgy(6))
DVds(1,1:3)=[abcabgy(1), abcabgy(2)*cos(degrees_to_radians*abcabgy(6)), abcabgy(3)*cos(degrees_to_radians*abcabgy(5))]
DVds(2,2:3)=[abcabgy(2)*sin(degrees_to_radians*abcabgy(6)), abcabgy(3) * ttrig ]
DVds(3,3) = abcabgy(3) * sqrt( -(ttrig**2) + (sin(degrees_to_radians*abcabgy(5)))**2 )
where(abs(DVds)<sceps_DP) DVds=zero
print*,'DVds : '
print*,DVds(1,:)
print*,DVds(2,:)
print*,DVds(3,:)
DVrs=zero
DVrs(1,1)=one/DVds(1,1)
DVrs(2,2)=one/DVds(2,2)
DVrs(3,3)=one/DVds(3,3)
DVrs(2,1)=-DVrs(1,1)*DVrs(2,2)*DVds(1,2)
DVrs(3,1)= DVrs(1,1)*DVrs(2,2)*DVrs(3,3)*(DVds(1,2)*DVds(2,3)-DVds(1,3)*DVds(2,2))
DVrs(3,2)=-DVrs(3,3)*DVrs(2,2)*DVds(2,3)
where(abs(DVrs)<sceps_DP) DVrs=zero
do i=1,3
  LVds(i)=sqrt(sum(DVrs(i:3,i)**2))
enddo
where(abs(LVds)<sceps_DP) LVds=zero


MTds(1,:) = [abcabgy(1)**2, abcabgy(1)*abcabgy(2)*cos(degrees_to_radians*abcabgy(6)), &
                           abcabgy(1)*abcabgy(3)*cos(degrees_to_radians*abcabgy(5))]
MTds(2,:) = [abcabgy(2)*abcabgy(1)*cos(degrees_to_radians*abcabgy(6)), abcabgy(2)**2, &
                           abcabgy(2)*abcabgy(3)*cos(degrees_to_radians*abcabgy(4))]
MTds(3,:) = [abcabgy(3)*abcabgy(1)*cos(degrees_to_radians*abcabgy(5)),  &
                           abcabgy(3)*abcabgy(2)*cos(degrees_to_radians*abcabgy(4)), &
                           abcabgy(3)**2]
where(abs(MTds)<sceps_DP) MTds=zero
do i=1,3
  print'("MTds [",i2,",:] = ",3(1x,g24.16))',i,MTds(i,:)
enddo

eigMT = EigMTDS(al_C=abcabgy(4),be_C=abcabgy(5),ga_C=abcabgy(6), a_C=abcabgy(1),b_C=abcabgy(2),c_C=abcabgy(3))
eigMT = sqrt(eigMT)
print*,'Metric tensor eigenvalues square roots: ',eigMT
print*,'NCL2MK: core/shell = ',NCL2MK
print*,'DCL2MK: core/shell = ',DCL2MK
NL_Kern=NCL2MK(1)-1
NL_Schale=NCL2MK(2)-1
allocate(xncelKS(2,0:NL_Kern,0:NL_Schale),ncelKS(2,0:NL_Kern,0:NL_Schale))
do ik=0,NL_Kern
  ncelKS(1,ik,:) = ik**3
  do is=0,NL_Schale
    ncelKS(2,ik,is) = (is+ik)**3 - ncelKS(1,ik,is)
  enddo
enddo
xncelKS=real(ncelKS,DP)

NCL2MK_MAX = sum(NCL2MK)-2
print*,'*********** NCL2MK_MAX   ',NCL2MK_MAX
print*,'*********** Delta_R,LVds ',Delta_R,LVds
Ncube1 = 1+ceiling((half+real(NCL2MK_MAX+2,DP))*Delta_R*LVds)
allocate(popd0V(0:NCL2MK_MAX),ndidi0V(0:NCL2MK_MAX),olastV(0:NCL2MK_MAX))

print'(a,2i4,3x,2(1x,"0 :",i4),i6)','NCL2MK(1:2), #_Kern, #_Schale, Layers_Max: ',NCL2MK,NCL2MK-1,NCL2MK_MAX

NCube2 = NCube1 * 2
print *,'Max cube dimensions (cells) : ',Ncube1,Ncube2

allocate(Cube1(ncelltypes,-Ncube1(1):Ncube1(1),-Ncube1(2):Ncube1(2),-Ncube1(3):Ncube1(3)),stat = allocOK)
if (allocOK/=0) then
  print*,'Cannot allocate Cube1, too many nodes! ',2*Ncube1+1,product(2*Ncube1+1)
  stop 'Cannot allocate Cube1, too many nodes! '
endif
Cube1=zero

allocate(Cube1r(-Ncube1(1):Ncube1(1),-Ncube1(2):Ncube1(2),-Ncube1(3):Ncube1(3)),stat = allocOK)
if (allocOK/=0) then
  print*,'Cannot allocate Cube1r, too many nodes! ',2*Ncube1+1,product(2*Ncube1+1)
  stop 'Cannot allocate Cube1r, too many nodes! '
endif
Cube1r=zero

allocate(Cube1i(-Ncube1(1):Ncube1(1),-Ncube1(2):Ncube1(2),-Ncube1(3):Ncube1(3)),stat = allocOK)
if (allocOK/=0) then
  print*,'Cannot allocate Cube1i, too many nodes! ',2*Ncube1+1,product(2*Ncube1+1)
  stop 'Cannot allocate Cube1i, too many nodes! '
endif
Cube1i=0

allocate(Cube2(ncellpairs,-Ncube2(1):Ncube2(1),-Ncube2(2):Ncube2(2),0:Ncube2(3)),stat = allocOK)
if (allocOK/=0) then
  print*,'Cannot allocate Cube2, too many nodes! ',[2,2,1]*Ncube2+1,(Ncube2(3)+1)*(2*Ncube2(1)+1)*(2*Ncube2(2)+1)
  stop 'Cannot allocate Cube2, too many nodes! '
endif
Cube2=zero
RRmax_a = Delta_R*real(NCL2MK_MAX,DP)
print*,'Max.radius [Angstroem] : ',RRmax_a
volu_clus_max = quter*Pi*RRmax_a*RRmax_a*RRmax_a
nu_dis_max = size(Cube2) / ncellpairs
print'(a,/,20x,3(1x,g18.10),i4)','Max. sphere volume, cell volume, WS-cell volume, # centers ',&
        volu_clus_max,cell_volume,cell_volume/real(ncenters_cell,DP),ncenters_cell
print'(a,/,20x,2(1x,g18.10))','N. cells in largest sphere = Max. sphere volume/cell volume, '//&
                              'Max. sphere volume/WS-cell volume', &
                               volu_clus_max/cell_volume, ncenters_cell*volu_clus_max/cell_volume
print*,'nu_dis_max = ',nu_dis_max,volu_clus_max/cell_volume,nu_dis_max/(volu_clus_max/cell_volume)




allocate(radat(nu_dis_max),mulat(nu_dis_max),mulats(0:nu_dis_max))

allocate( finu(0:NCL2MK(1)-1,0:NCL2MK(2)-1) )
do j_k=0,NCL2MK(1)-1; do j_s=0,NCL2MK(2)-1
  finu(j_k,j_s)(:)=''
  write(finu(j_k,j_s)(1:14),'("_k",i3.3,"_s",i3.3,"_",a3)')j_k,j_s,clushape
 !! j_k = # effective kern layers; j_s = # effective schale layers
enddo;enddo

! fill boxes
call CPU_TIME(ts0)
do i1=-NCube1(1),NCube1(1)
  aux1(1)=i1
  do i2=-NCube1(2),NCube1(2)
    aux1(2)=i2
    do i3=-NCube1(3),NCube1(3)
      aux1(3)=i3
      ddd2 = sum(aux1*matmul(MTds,aux1))
      ddd=sqrt(ddd2)
      Cube1r(i1,i2,i3)=ddd
    enddo
  enddo
enddo
call CPU_TIME(ts1)
print'(a,1x,f20.6," sec.")','Time for Cube1r: ',ts1-ts0
print*,'Number of atom species/cell, all species, all pairs: ',n_sp_atom_V,naspe_all,napair_all
do ice=1,ncelltypes
  print'(a,i3,5x,15i4)','Species address - cell # ',ice,address_atsp(:,ice)
enddo

diamax = two*(half+real(sum(NCL2MK)+1,DP))*Delta_R  !ten*DCL2MK(1)+two*maxval(eigMT)
dimens = max(nbeta,N1cont)+CEILING(diamax/Deltas)
dimemax=maxval(dimens(Ndelta1:Ndelta2))
print*,'samv : ',Ndelta1,Ndelta2,Deltas(Ndelta1:Ndelta2),dimens(Ndelta1:Ndelta2),dimemax,Delta_R
if (allocated(samv)) deallocate(samv)
allocate(samv(dimemax,napair_all,Ndelta1:Ndelta2,0:NCL2MK(2)-1))
samv=zero


!_____________ ALLOCATE LOCAL VECTORS
 
 if (allocated(nat0)) deallocate(nat0)
 if (allocated(xnat0)) deallocate(xnat0)
 if (allocated(ndi0)) deallocate(ndi0)
 if (allocated(summul0)) deallocate(summul0)
 if (allocated(termcon0)) deallocate(termcon0)
 if (allocated(patte_len)) deallocate(patte_len)
 if (allocated(patte_invmul)) deallocate(patte_invmul)
 mmmp=maxval(numpatv(1:napair_all,1:ncellpairs))
 allocate(nat0(1:naspe_all,0:NCL2MK(1)-1,0:NCL2MK(2)-1),xnat0(1:naspe_all,0:NCL2MK(1)-1,0:NCL2MK(2)-1),& 
          ndi0(1:napair_all,0:NCL2MK(1)-1,0:NCL2MK(2)-1),&
          summul0(1:napair_all,0:NCL2MK(1)-1,0:NCL2MK(2)-1), termcon0(napair_all,0:NCL2MK(1)-1,0:NCL2MK(2)-1),&
          patte_len(mmmp,napair_all,ncellpairs),patte_invmul(mmmp,napair_all,ncellpairs))
patte_len=zero
patte_invmul=zero
!__ Fill up patte_* local vectors
kcelp=0
do icel1=1,ncelltypes; do icel2=icel1,ncelltypes
  kcelp=kcelp+1
  do ipa=1,napair_all
    do ix=1,numpatv(ipa,kcelp)
      xlens = sum(patte(1:3,ix,ipa,kcelp)*matmul(MTds,patte(1:3,ix,ipa,kcelp)))
      patte_len(ix,ipa,kcelp) = xlens
      fk=two
      if ((abs(xlens)<sceps_DP .and. icel1==icel2).and.is1==is2) fk=one
      patte_invmul(ix,ipa,kcelp) = fk
    enddo
  enddo
enddo; enddo

do ikern=0,NCL2MK(1)-1; do ischale=0,NCL2MK(2)-1
  summul0(:,ikern,ischale)=zero
  termcon0(:,ikern,ischale)=zero
  do k=1,napair_all
    nc_kern   = (ikern)**3
    nc_schale = (ikern+ischale)**3 - nc_kern
    ndi0(k,ikern,ischale)=8*(nc_kern*numpatv(k,1)+nc_schale*numpatv(k,3)) ! very approximately, anyway useless
  enddo
enddo; enddo


  print*,'Global pairs: ',napair_all
  print*,glopair(1,:)
  print*,glopair(2,:)
  print*,Z_at_glo(glopair(1,:))
  print*,Z_at_glo(glopair(2,:))
do k=1,napair_all
  print*,'******* PAIR ',k
  ipa=k
  is1=glopair(1,k)
  is2=glopair(2,k)
  izz1=Z_at_glo(is1)
  izz2=Z_at_glo(is2)
  print*,'Pair ',k,' : Species index    ',is1,is2
  print*,'Pair ',k,' : Species Zs       ',izz1,izz2
  print*,'Pair ',k,' : Species # at.    ',nat_glo(is1),nat_glo(is2)
  print*,'Pair ',k,' : Species sum occ. ',xnat_glo(is1),xnat_glo(is2)
enddo
!stop

!finqi

yoll=sceps_DP
ndidi=0
radat=zero
mulat=zero
mulats=zero
Cube1i = 0
call CPU_TIME(ts0)
do i1=-Ncube1(1),Ncube1(1)
  aux1(1)=i1
  do i2=-Ncube1(2),Ncube1(2)
    aux1(2)=i2
    do i3=-Ncube1(3),Ncube1(3)
      aux1(3)=i3
      ddd = Cube1r(i1,i2,i3)
      IF (ndidi<1) THEN
        ndidi = 1
        radat(ndidi) = ddd
        mulat(ndidi) = one
        print*,'First',1,ddd,Ncube1
        cycle
      ENDIF
      
      call PLACER(ddd,radat(1:ndidi),n2place,isiteq,yoll)
      
      IF (isiteq==1) THEN
        mulat(n2place) = mulat(n2place)+one
      ELSE
        IF (ndidi==nu_dis_max) THEN
          print*,'nu_dis_max not enough! ',nu_dis_max
          stop 'nu_dis_max not enough! '
        ENDIF
        radat(n2place+2:ndidi+1) = radat(n2place+1:ndidi)
        mulat(n2place+2:ndidi+1) = mulat(n2place+1:ndidi)
        radat(n2place+1) = ddd
        mulat(n2place+1) = one
        ndidi=ndidi+1
        print*,'New',ndidi,ddd
      ENDIF
    enddo
  enddo
enddo
call CPU_TIME(ts1)
print'(a,1x,f20.6," sec.")','Time for radat/mulat: ',ts1-ts0
print*,'nu_dis_max = ',nu_dis_max,ndidi,volu_clus_max/cell_volume,nu_dis_max/(volu_clus_max/cell_volume)
if (verbose) print*,'Allocating check ',ndidi,nu_dis_max,size(Cube1)
where(abs(mulat)<=yoll) mulat = zero
Cube1i=0
do i=1,ndidi
  ddd=radat(i)
  where(abs(Cube1r-ddd)<yoll) Cube1i=i
enddo
n00=count(Cube1i==0)
if (n00>0) stop 'Unfilled index box'
mulats=zero
do ja=1,ndidi
  mulats(ja)=mulats(ja-1)+mulat(ja)
enddo
ndidi0V=0
popd0V=0
olastV=zero
print*,NCL2MK_MAX,size(olastV),size(ndidi0V)
do i=1,NCL2MK_MAX
  xwish = real(i**3,DP) / real(ncenters_cell,DP)
  do ja=1,ndidi
    if (abs(mulats(ja)-xwish) <= sceps_DP) then
      ndidi0V(i)=ja
      olastV(i)=one
      popd0V(i)=mulat(ja)
      print*,'Total R ',i,i*Delta_R,ndidi0V(i),olastV(i)
      do jx=1,ja
        print*,jx,mulats(jx),xwish
      enddo
      exit
    else if (mulats(ja)-xwish > sceps_DP) then
      ndidi0V(i)=ja
      excess=mulats(ja)-xwish
      olastV(i) = (xwish-mulats(ja-1))/mulat(ja)
      popd0V(i)=mulat(ja)
      print*,'Total R ',i,i*Delta_R,ndidi0V(i),olastV(i)
      do jx=1,ja
        print*,jx,mulats(jx),xwish
      enddo
      exit
    endif
  enddo
enddo
call CPU_TIME(ts2)
print'(a,1x,f20.6," sec.")','Time for ndidi0V/olastV: ',ts2-ts1



Cube1=zero
jlower_K=NCL2MK(1)-1
jupper_K=NCL2MK(1)-1

if (allsizes) then
  jlower_K=0
  if (modulo(NCL2MK(1)-1,dr_step)>0) then
    jupper_K=(1+((NCL2MK(1)-1)/dr_step))*dr_step
  endif
endif
print*,'Doing j_KERN from ',jlower_K,' to ',jupper_K,' by ',dr_step
jlower_S=NCL2MK(2)-1
jupper_S=NCL2MK(2)-1

if (allsizes) then
  jlower_S=0
  if (modulo(NCL2MK(2)-1,dr_step)>0) then
    jupper_S=(1+((NCL2MK(2)-1)/dr_step))*dr_step
  endif
endif
print*,'Doing j_SCHALE from ',jlower_S,' to ',jupper_S,' by ',dr_step




!istop=1
!if (istop==1) stop

call CPU_TIME(tstart)
t0=tstart
DO k_kern=jlower_K,jupper_K,dr_step
  samv=zero
  
  radius_K = real(k_kern,DP)*Delta_R
  radius2_K = radius_K**2+sceps_DP
  diameter_K = two*radius_K
  nlap_K = (k_kern)**3
  xlap_K = real(nlap_K,DP)
  mmax_K = 1+ceiling((half+real(k_kern,DP))*Delta_R*LVds)
    !!min(NCube1,1+ceiling( (radius + sceps_DP) / minval(eigMT) ))
  xwish_K = real((k_kern)**3,DP) / real(ncenters_cell,DP) ! # of whole unit cells
  
  Rmax_K=two*radius_K
  Rmax2_K=Rmax_K**2
  
  
  Cube1=zero
  Cube2=zero
  
  ndidi0_k = ndidi0V(k_kern)
  xlast_K  = olastV(k_kern)
  nlast_K  = popd0V(k_kern)
  where (Cube1i==ndidi0_k) 
    Cube1(1,:,:,:) = xlast_K
  elsewhere (Cube1i>ndidi0_K) 
    Cube1(1,:,:,:) = zero
  elsewhere (Cube1i<ndidi0_K) 
    Cube1(1,:,:,:) = one
  end where
  if (verbose) then
    xxx=sum(Cube1(1,:,:,:))
    print*,'Cutting Kern : ',k_kern,ndidi0_k, xxx,nint(xxx*ncenters_cell)
  endif
!_____________ Make dist  
  ! Shell occulpancy
  xlast2_K = xlast_k**2
  C_xlast2_K = (one-xlast_k)**2
  P_xlast2_K = (one-xlast_k)*xlast_k*two
  xcuq = xlast2_K * nlast_K + real(count(Cube1i<ndidi0_K),DP)
  nlatpoi_clu(1) = count(Cube1(1,:,:,:)>yoll)
  xlatpoi_clu(1) = xlap_K
  
  DO k_schale=jlower_S,jupper_S,dr_step
    k_all=k_schale+k_kern
    radius_S = real(k_all,DP)*Delta_R
    radius2_S = radius_S**2+sceps_DP
    diameter_S = two*radius_S
    nlap_S = (k_all**3)-(k_kern**3)
    xlap_S = real(nlap_S,DP)
    mmax_S = 1+ceiling((half+real(k_all+2,DP))*Delta_R*LVds)
      !!min(NCube1,1+ceiling( (radius + sceps_DP) / minval(eigMT) ))
    xwish_S = real((k_all)**3-k_kern**3,DP) / real(ncenters_cell,DP) ! # of whole unit cells
    print'("k_kern,k_schale,k_all,mmax_S : ",3i6,2x,3i5, 4(1x,g22.16))',k_kern,k_schale,k_all,mmax_S, &
       real(k_all,DP)*Delta_R, real(mmax_S,DP)*abcabg_loc(1)
    
    Rmax_S=two*radius_S
    Rmax2_S=Rmax_S**2
    
    ndidi0_S = ndidi0V(k_all)
    xlast_S  = olastV(k_all)
    nlast_S  = popd0V(k_all)
    print'("k_kern,k_schale,k_all,xlast_S,xlast_K,nlast_S,nlast_K : ",3i6,2x,2(1x,g22.16),4i5 )', &
           k_kern,k_schale,k_all,xlast_S,xlast_K,nlast_S,nlast_K,ndidi0_K,ndidi0_S
    Cube2=zero
    Cube2(1,0,0,0) = xcuq
    if (k_schale==0) then
      Cube1(2,:,:,:) = zero
      Cube2(2:3,0,0,0)=zero
    else
      if (ndidi0_S>ndidi0_K) then
        where (Cube1i==ndidi0_S)
          Cube1(2,:,:,:) = xlast_S
        elsewhere (Cube1i>ndidi0_S) 
          Cube1(2,:,:,:) = zero
        elsewhere (Cube1i<ndidi0_S)
          Cube1(2,:,:,:) = one-Cube1(1,:,:,:)
        end where
        xlast2_S = xlast_S**2
        Cube2(3,0,0,0) = xlast2_S * nlast_S + C_xlast2_K * nlast_K + real(count(Cube1i<ndidi0_S.and.Cube1i>ndidi0_K),DP)
        Cube2(2,0,0,0) = P_xlast2_K * nlast_K
      else if (ndidi0_S==ndidi0_K) then
        where (Cube1i==ndidi0_S)
          Cube1(2,:,:,:) = xlast_S-xlast_K
        elsewhere (Cube1i/=ndidi0_S) 
          Cube1(2,:,:,:) = zero
        end where
        xlast2_S = (xlast_S-xlast_K)**2
        Cube2(3,0,0,0) = xlast2_S * nlast_S
        Cube2(2,0,0,0) = (xlast_S-xlast_K) * xlast_K * nlast_S * two
      endif
    endif
    if (verbose) then
      xxx=sum(Cube1(1,:,:,:))
      xxx2=sum(Cube1(2,:,:,:))
      print*,'Cutting Schale : ',k_kern,k_schale,k_all,ndidi0_S, xxx,xxx2,nint([xxx,xxx2,(xxx+xxx2)]*ncenters_cell)
    endif
    
    nlatpoi_clu(2) = count(Cube1(2,:,:,:)>yoll)
    xlatpoi_clu(2) = xlap_S
    print'(a,2i4,1x,2(1x,f16.8))','#Core-#Shell : Sum lattice point occ.s = ',k_kern,k_schale,xlatpoi_clu(:)
    xnat0(:,k_kern,k_schale)=zero
    nat0(:,k_kern,k_schale)=0
    do k=1,naspe_all
      do icel=1,ncelltypes
        kloc = return_address_atsp(k,icel)
        if (kloc==0) cycle
        xnat0(k,k_kern,k_schale) = xnat0(k,k_kern,k_schale) + &
                                   xlatpoi_clu(icel) * sum(Celty(icel)%ato(1:Celty(icel)%nat(kloc),kloc)) &
                                   / ncenters_cell
        nat0(k,k_kern,k_schale)  = nat0(k,k_kern,k_schale) + nlatpoi_clu(icel) * Celty(icel)%nat(kloc)
      enddo
    enddo
    if (k_kern==0.and.k_schale==0) then
      do k=1,naspe_all
        print*,'Species / all species, #Core, #Shell',k,naspe_all,k_kern,k_schale
        do icel=1,ncelltypes
          kloc = return_address_atsp(k,icel)
          if (kloc==0) cycle
          print'(a,i2,a,i3,a,i2,1x,f16.8,i6)','Cell type ',icel,' : Global species #',k,' is mapped onto ',kloc,&
                                               xlatpoi_clu(icel), Celty(icel)%nat(kloc)
          xnat0(k,k_kern,k_schale) = xnat0(k,k_kern,k_schale) + &
                                     xlatpoi_clu(icel) * sum(Celty(icel)%ato(1:Celty(icel)%nat(kloc),kloc)) &
                                     / ncenters_cell
          nat0(k,k_kern,k_schale)  = nat0(k,k_kern,k_schale) + nlatpoi_clu(icel) * Celty(icel)%nat(kloc)
          ssssx=0.d0
          do iiii=1,Celty(icel)%nat(kloc)
            ssssx=ssssx+Celty(icel)%ato(iiii,kloc)
            print'(i4,2i2,2(1x,f16.8))',iiii,icel,kloc,Celty(icel)%ato(iiii,kloc),ssssx
          enddo
        enddo
      enddo
    endif
  
    call CPU_TIME(t1)
    !print*,'*** ',two*radius/ten,k_kern,NCL2MK(1),c000a1,c000b1,xlast,Cube2(0,0,0),t1-t0,t1-tstart
    t0=t1
    mmax2=2*mmax_S
    mmax=mmax_S
    do i3=0,mmax2(3)
      aux1(3)=i3
      ii2=-mmax2(2)
      if (i3==0) ii2=0
      lobo(3)=max(-mmax(3),-mmax(3)-i3)
      upbo(3)=min( mmax(3), mmax(3)-i3)
      do i2=ii2,mmax2(2)
        aux1(2)=i2
        ii1=-mmax2(1)
        if (i3==0.and.i2==0) ii1=1
        lobo(2)=max(-mmax(2),-mmax(2)-i2)
        upbo(2)=min( mmax(2), mmax(2)-i2)
        do i1=ii1,mmax2(1)
          aux1(1)=i1
          iimax=maxval(abs(aux1))
          lobo(1)=max(-mmax(1),-mmax(1)-i1)
          upbo(1)=min( mmax(1), mmax(1)-i1)
          soccp=zero
          do j1=lobo(1),upbo(1)
            do j2=lobo(2),upbo(2)
              do j3=lobo(3),upbo(3)
                ooo1=Cube1(:,j1,j2,j3)
                if (maxval(abs(ooo1))<=yoll) cycle
                ooo2=Cube1(:,j1+i1,j2+i2,j3+i3)
                if (maxval(abs(ooo2))<=yoll) cycle
                soccp(1)=soccp(1)+ooo2(1)*ooo1(1)
                soccp(2)=soccp(2)+half*(ooo2(1)*ooo1(2) + ooo2(2)*ooo1(1))
                soccp(3)=soccp(3)+ooo2(2)*ooo1(2)
              enddo
            enddo
          enddo
          Cube2(:,i1,i2,i3)=Cube2(:,i1,i2,i3)+soccp
        enddo
      enddo
    enddo
    call CPU_TIME(t222)
    do iqb=1,3
      ssssx=sum(Cube2(iqb,:,:,:))
      ssss0=Cube2(iqb,0,0,0)
      if (modulo(iqb,2)==1) then
        soccp(iqb) = ssssx-half*ssss0
      else
        soccp(iqb) = half*(ssssx+ssss0)
      endif
      print'(a,3i3,3(1x,f16.8))','Summing: soccp',k_kern,k_schale,iqb,32*[ssssx,ssss0,soccp(iqb)]
    enddo
    soccp=soccp * two
    print'(a,2i3,7(1x,f17.6))','Sums qbe',k_kern,k_schale,four*[sum(Cube1(1,:,:,:)),sum(Cube1(2,:,:,:))],&
       16*soccp,16*sum(soccp),4*sqrt(sum(soccp))
    if (output_xyz) then
      lfile_inp=len_trim(file_inp)
      fxyz=''
      fxyz(1:lfile_inp)=file_inp(1:lfile_inp)
      print'(a,2x,2i6)',fxyz(1:lfile_inp),lfile_inp,-(lfile_inp-3)+lfile_inp+9+1
      write(ccccc,'("_k",i3.3,"_s",i3.3,".xyz")')k_kern+1,k_schale+1
      !write(fxyz(lfile_inp-3:lfile_inp+9),'("_k",i3.3,"_s",i3.3,".xyz")')k_kern+1,k_schale+1
      iux=find_unit()
      open(iux,status='replace',file=fxyz(1:lfile_inp)//ccccc)
      write(iux,*) sum(nat0(:,k_kern,k_schale))
      write(iux,*) fxyz(1:lfile_inp),k_kern,k_schale
      do icel=1,ncelltypes
        nacu(icel) = count(Cube1(icel,:,:,:)>yoll)
      enddo
      print'(a,2i4,6x,2i6)','Inventory: ',k_kern,k_schale,nacu
      do k=1,naspe_all
        print*,'Species list # k ',k
        do icel=1,ncelltypes
          kloc = return_address_atsp(k,icel)
          if (kloc==0) then
            print*,'k,icel,kloc, Z, N_Z, NLat',k,icel,kloc,0,0,nacu(icel)
            cycle
          else
            izz = Celty(icel)%Z_at(kloc)
            print*,'k,icel,kloc, Z, N_Z, NLat',k,icel,kloc,izz,Celty(icel)%nat(kloc),nacu(icel)
          endif
          do i1=-Ncube1(1),Ncube1(1)
            aux1(1)=i1
            do i2=-Ncube1(2),Ncube1(2)
              aux1(2)=i2
              do i3=-Ncube1(3),Ncube1(3)
                aux1(3)=i3
                ooo=Cube1(icel,i1,i2,i3)
                if (abs(ooo)<=yoll) cycle
                do ia=1,Celty(icel)%nat(kloc)
                  Raux1=aux1+Celty(icel)%atx(:,ia,kloc)
                  oooa=ooo*Celty(icel)%ato(ia,kloc)
                  xyzP=matmul(DVds,Raux1)
                  write(iux,'(a2,1x,3(1x,g14.8),1x,g12.6)')symb_of_Z(izz),xyzP,oooa
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      close(iux)
    endif
    
    rlatmax=zero
    do i3=0,mmax2(3)
      aux1(3)=i3
      Raux1(3)=real(i3,DP)
      ii2=-mmax2(2)
      if (i3==0) ii2=0
      do i2=ii2,mmax2(2)
        Raux1(2)=real(i2,DP)
        aux1(2)=i2
        ii1=-mmax2(1)
        if (i3==0.and.i2==0) ii1=0
        do i1=ii1,mmax2(1)
          aux1(1)=i1
          Raux1(1)=real(i1,DP)
          ooo_latV=Cube2(:,i1,i2,i3)
    !      if (maxval(abs(ooo_latV))<yoll) cycle
          aux_latvec=matmul(MTds,Raux1)
          d0lat=sum(Raux1*aux_latvec) ! RL**2
          rlatmax=max(d0lat,rlatmax)
          aux_latvec=aux_latvec*two
          isoriginL = ALL(aux1==0)
          fffmuL=two
          if (isoriginL) fffmuL=one
          iiimuL=nint(fffmuL)
  !       
          PAICE:do ipce=1,ncellpairs
            if (abs(ooo_latV(ipce)) < yoll) cycle PAICE
            icel1=paolopa(1,ipce)
            icel2=paolopa(2,ipce)
            BLACKMORE:DO ipa=1,napair_all
              is1=glopair(1,ipa)
              is2=glopair(2,ipa)
              nlord=1
              if ((is1/=is2).and.(icel1/=icel2)) nlord=2
              LORD:do klord=1,nlord
                if (klord > 1) then
                  ixc=icel1
                  icel1=icel2
                  icel2=ixc
                endif
                isl1=return_address_atsp(is1,icel1)
                isl2=return_address_atsp(is2,icel2)
                if (min(isl1,isl2)==0) then
                  cycle LORD
                  done_once = .false.
                else
                  done_once = .true.
                endif
                izz1=Celty(icel1)%Z_at(isl1)
                izz2=Celty(icel2)%Z_at(isl2)
                xminallow = at_radii(izz1)+at_radii(izz2)
                GILLAN:do jx=1,numpatv(ipa,ipce)
                  ooo_this = patte(4,jx,ipa,ipce) * ooo_latV(ipce)
                  if (abs(ooo_this)<sceps_DP) cycle GILLAN
                  fffmuP = patte_invmul(jx,ipa,ipce)
                  isoriginP = (abs(patte_len(jx,ipa,ipce))<sceps_DP).and.(icel1==icel2) !(nint(fffmuP)==1)
                  ooo_this = ooo_this *fffmuL  !* fffmuP!max(fffmuP,fffmuL)
                  d000a=d0lat + patte_len(jx,ipa,ipce)
                  d000b=sum(aux_latvec*patte(1:3,jx,ipa,ipce))
   !               if (icel1/=icel2.and.is1/=is2) ooo_this=ooo_this*two
                  full_dis_sq1 = d000a + d000b
                  
                  d0s=sqrt(max(zero,full_dis_sq1))
                  if (d0s < xminallow ) then            ! too short distance, set it to 0
                     termcon0(ipa,k_kern,k_schale) = termcon0(ipa,k_kern,k_schale) + ooo_this
                  else
                    call SAM_ONE_TOALL(d0=d0s)
                    summul0(ipa,k_kern,k_schale)=summul0(ipa,k_kern,k_schale)+ooo_this
                    do ksam=Ndelta1,Ndelta2
                      jsum1=max(1,centersum(ksam)-N1cont)
                      jrr1=jsum1-centersum(ksam)
                      jsum2=min(dimemax,centersum(ksam)+N1cont)
                      jrr2=jsum2-centersum(ksam)
                      samv(jsum1:jsum2,ipa,ksam,k_schale) = samv(jsum1:jsum2,ipa,ksam,k_schale) + &
                                                       for_samv_P(jrr1:jrr2,ksam) * ooo_this
                      if (bnd_N(2,ksam) <= 0) cycle
                      samv(bnd_N(1,ksam):bnd_N(2,ksam),ipa,ksam,k_schale) = &
                        samv(bnd_N(1,ksam):bnd_N(2,ksam),ipa,ksam,k_schale) + &
                           for_samv_N(bnd_N(1,ksam):bnd_N(2,ksam),ksam) * ooo_this
                    enddo
                  endif
          !______________________ difference equals sum if one is zero
          
  !                if ((isoriginL.or.isoriginP).and.(icel1==icel2)) cycle
   !               if (((isoriginL.and.isoriginP).and.(is1==is2)).and.(icel1==icel2)) cycle
                  if (((isoriginP).and.(is1==is2)).and.(icel1==icel2)) cycle
   !!               if (((isoriginP).and.(is1==is2))) cycle
                  
                  full_dis_sq1 = d000a - d000b
                  
                  d0s=sqrt(max(zero,full_dis_sq1))
      !            if (k_kern==0) print*,'in',2,ipa,d0s,ooo_this,ooo_lat,patte(4,jx,ipa,1)
                  if (d0s < xminallow) then            ! too short distance, set it to 0
                     termcon0(ipa,k_kern,k_schale) = termcon0(ipa,k_kern,k_schale)+ooo_this
                  else
                    call SAM_ONE_TOALL(d0=d0s)
                    summul0(ipa,k_kern,k_schale)=summul0(ipa,k_kern,k_schale)+ooo_this
                    do ksam=Ndelta1,Ndelta2
                      jsum1=max(1,centersum(ksam)-N1cont)
                      jrr1=jsum1-centersum(ksam)
                      jsum2=min(dimemax,centersum(ksam)+N1cont)
                      jrr2=jsum2-centersum(ksam)
                      samv(jsum1:jsum2,ipa,ksam,k_schale) = samv(jsum1:jsum2,ipa,ksam,k_schale) + &
                                                       for_samv_P(jrr1:jrr2,ksam) * ooo_this
                      if (bnd_N(2,ksam) <= 0) cycle
                      samv(bnd_N(1,ksam):bnd_N(2,ksam),ipa,ksam,k_schale) = &
                        samv(bnd_N(1,ksam):bnd_N(2,ksam),ipa,ksam,k_schale) + &
                          for_samv_N(bnd_N(1,ksam):bnd_N(2,ksam),ksam) * ooo_this
                    enddo
                  endif
                enddo GILLAN
                if (done_once) exit LORD
              enddo LORD
            enddo BLACKMORE
          enddo PAICE          
        enddo
      enddo
    enddo
    print*,'sumqbe2 :',k_kern,k_schale,sum(Cube2(1,:,:,:)),sum(Cube2(2,:,:,:)),sum(Cube2(3,:,:,:))
    print*,'000qbe2 :',k_kern,k_schale,Cube2(:,0,0,0)
    print*,'sumqbe1 :',k_kern,k_schale,sum(Cube1(1,:,:,:)),sum(Cube1(2,:,:,:))
    print*,'sumqbe1s:',k_kern,k_schale,sum(Cube1(1,:,:,:)**2),sum(Cube1(2,:,:,:)**2)
!  Celty(1)%nat=nat0(:,k_kern)
!  Celty(1)%xnat=xnat0(:,k_kern)
!  Celty(1)%ndi=ndi0(:,k_kern,k_schale)
!  Celty(1)%summul=summul0(:,k_kern,k_schale)
!  Celty(1)%termcon_all  = termcon0(:,k_kern,k_schale)
!  Celty(1)%termcon_ineq = termcon0(Celty(1)%point_neqpair(:),k_kern,k_schale)
!  Celty(1)%termcon      = termcon0(Celty(1)%point_eqpair(:),k_kern,k_schale)
!  
!  Z_at_glo = Celty(1)%Z_at

    nat_glo = nat0(:,k_kern,k_schale)
    xnat_glo = xnat0(:,k_kern,k_schale)
    zappa_glo = glopair
    ndi_glo = ndi0(:,k_kern,k_schale)
    termcon_all_glo = termcon0(:,k_kern,k_schale)
    ksp=0
    kin=0
    do iss1=1,naspe_all
      MORSE:do iss2=iss1,naspe_all
        ksp=ksp+1
        if (iss1/=iss2) then
          kin=kin+1
          termcon_ineq_glo(kin) = termcon0(ksp,k_kern,k_schale)
        else
          termcon_glo(iss1) = termcon0(ksp,k_kern,k_schale)
        endif
      enddo MORSE
    enddo
    summul_glo = summul0(:,k_kern,k_schale)
    numcells_byphase(1:2) = [nlap_K,nlap_S]
    
    Size_Abscissa=real([k_kern,k_schale],DP)
  
    Actual_Diam = [growing_diam(k_kern+1,1),growing_diam(k_schale+1,2)]
    
    call OUTSAMP(isam=k_schale,finumb=namestr(1:lnamestr)//finu(k_kern,k_schale)(1:14)//'.smp', Bpre=zero)
    print*,'Rlat max ',rlatmax
  !if (k_kern==0) stop
  
  enddo
enddo
iucon=find_unit()
open(iucon,status='replace',file='check_termcon.out')
write(iucon,*)'#',NCL2MK(1),naspe_all,napair_all
DO k_kern=jlower_K,jupper_K,dr_step
  DO k_schale=jlower_S,jupper_S,dr_step
    write(iucon,*)k_kern,k_schale,nat0(:,k_kern,k_schale),xnat0(:,k_kern,k_schale)
    write(iucon,*)' A   ',summul0(:,k_kern,k_schale)
    write(iucon,*)' B  ',termcon0(:,k_kern,k_schale)
    write(iucon,*)' C  ',termcon0(:,k_kern,k_schale)+summul0(:,k_kern,k_schale),xnat0(:,k_kern,k_schale)**2
    write(iucon,*)' D  ',Celty(1)%point_neqpair, Celty(1)%point_eqpair
    ipa=0
    do iss1=1,naspe_all
      do iss2=iss1,naspe_all
        ipa=ipa+1
        if (iss1==iss2) then
          sadd=xnat0(iss1,k_kern,k_schale)**2
          ckval=termcon0(ipa,k_kern,k_schale)-sadd
        else
          sadd=xnat0(iss1,k_kern,k_schale)*xnat0(iss2,k_kern,k_schale)*two
          ckval=termcon0(ipa,k_kern,k_schale)-sadd
        endif
        write(iucon,'(a,2i5,3i3,1x,g18.8,5x,4(1x,f14.4))')'soll0: ',k_kern,k_schale,ipa,iss1,iss2,ckval, &
              sadd,termcon0(ipa,k_kern,k_schale),xnat0(iss1,k_kern,k_schale),xnat0(iss2,k_kern,k_schale)
      enddo
    enddo
  enddo
enddo
close(iucon)
end program Paint_It_Black
