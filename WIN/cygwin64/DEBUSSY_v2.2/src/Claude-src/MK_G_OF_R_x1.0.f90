program dopdf1
use nano_deftyp
use rhapsody_in_blue
implicit real(DP)(a-h,o-z),integer(I4B)(i-n)
integer(I4B),parameter :: n2r = 3
real(DP),allocatable :: i_expt(:),tt(:),sig(:),ibroad(:),scale_I_vec(:)
real(DP),allocatable :: r(:),sigGofR(:,:),GofR(:,:)
real(DP),allocatable :: fractional_comp(:),qcuts(:)
integer(I4B),allocatable :: z_comp(:),lfn(:)
real(DP)     :: CoeCol(2,n2r),pola(2)
integer(I4B) :: col_index(2,n2r)

integer(I4B)             :: fi,nfil,fk,xray1_neutron2_el3,rescaling_mode,iotf,degbk,type_bkg
real(DP) :: r_minmaxstep(3),tt_minmax(2),ubroad,x_y_e(n2r), tail_fraction=0.05d0,r_degbk, bkg_lowr
character(999) :: rl,smp,tai,gug
character(len=1) :: flag_scal,flag_rad
character(999),allocatable :: fn(:)
logical                    :: sigI, do_normalize, scale_abs, ext_rminstepmax, ext_ttminmax, &
                              is_difference
logical :: readcolfile=.false., given_mass_dens=.false.
character(len=1) :: letter_bkg,y_soper_tophat


real(DP),allocatable :: Iofq(:,:), aux(:)
logical,allocatable :: essex(:)
!!__________________________


print*, '                                         ------------------------------------'
print*,'                                                 DebUsSy Suite v2.X   '
print*, '                                         ------------------------------------'
print*, ' '
print*,'     Running MK_G_OF_R Program    '
print*, ' '



sigI = .false.
do_normalize = .false.
scale_abs = .false.
ubroad=zero
pola(1)=zero
pola(2)=zero

fi= 0
iu11=find_unit()
open(iu11,status='old',action='read',file='dopdf.inp',iostat=ierr)
    IF (ierr /=0) THEN
        print*, ' Error opening input file: dopdf.inp '
       STOP
     ENDIF
binwidth=zero
r_minmaxstep=-one
extscale = one
ext_rminstepmax=.false.
xray1_neutron2_el3=1
tt_minmax=[zero, 180.d0]
ncut_q = 1
is_cut_rel=1
nfil=0
letter_bkg='c'
degbk=0
type_bkg=1
bkg_lowr=0.d0

!___ default : read only 2 cols. x and y
ncolfile=2
col_index=0; CoeCol=zero

col_index(1,:)=[1,2,0]
col_index(2,:)=[0,0,0]
CoeCol(1,:) =[ one, one,zero]
CoeCol(2,:) =[zero,zero,zero]


is_difference=.false.
do
  read(iu11,'(a)',iostat=io)rl
  if (io/=0) exit
  call CLEAN_line(rl)
  rl=trim(adjustl(rl))
  ll=len_trim(rl)
  if (ll==0.or.rl(1:1)=='!'.or.rl(1:1)=='>'.or.rl(1:1)=='%') cycle
  if (rl(1:4)=='NFIL') then
    read(rl(5:ll),*)nfil
    allocate(fn(nfil),lfn(nfil))
  else if (rl(1:4)=='SIGI') then  !__________ SIGI is now redundant
    sigI = .true.
    !___ default for this case
    ncolfile=3
    col_index(1,:)=[1,2,3]
    col_index(2,:)=[0,0,0]
    CoeCol(1,:) =[ one, one, one]
    CoeCol(2,:) =[zero,zero,zero]
  else if (rl(1:5)=='NCOLS') then
    read(rl(6:ll),*)ncolfile
    ! ncolfile is the number of real column in input file
    readcolfile=.true.
  else if (rl(1:5)=='TTCOL') then
    call RCXXX(a=rl(6:ll),numc=col_index(:,1),coec=CoeCol(:,1))
    !read(rl(6:ll),*) colfile(1)
  else if (rl(1:5)=='I_COL') then
    call RCXXX(a=rl(6:ll),numc=col_index(:,2),coec=CoeCol(:,2))
    !read(rl(6:ll),*) colfile(2)
  else if (rl(1:5)=='E_COL') then
    call RCXXX(a=rl(6:ll),numc=col_index(:,3),coec=CoeCol(:,3))
  else if (rl(1:4)=='FILE') then
    fi= fi+1
    fn(fi)=trim(adjustl(rl(5:ll)))
    lfn(fi)=len_trim(fn(fi))
    print'("*",a,"*")',fn(fi)(1:lfn(fi))
  else if (rl(1:4)=='XRNT') then
    tai=trim(adjustl(rl(5:ll)))
    ltai=len_trim(tai)
    read(tai(1:ltai),*) xray1_neutron2_el3
    if (xray1_neutron2_el3==1) then
      flag_rad='x'
    else if (xray1_neutron2_el3==2) then
      flag_rad='n'
    else if (xray1_neutron2_el3==3) then
      flag_rad='e'
    endif
  else if (rl(1:4)=='RAYS') then
    tai=trim(adjustl(rl(5:ll)))
    ltai=len_trim(tai)
    flag_rad=tai(1:1)
  else if (rl(1:4)=='POLA') then
    read(rl(5:ll),*)pola
  else if (rl(1:5)=='VALEN') then
    gug=trim(adjustl(rl(6:ll)))
    mod_valence = (gug(1:1)=='y'.or.(gug(1:1)=='Y'.or.(gug(1:1)=='1')))
    lgug=len_trim(gug)
    if (lgug>1) then
      !!___radius_fraction is used to make a Gaussian valence broadening with HWHM = radius_fraction * at_radii
      read(gug(2:lgug),*,iostat=iotf) radius_fraction
      if (iotf/=0) then
        radius_fraction=2.0d0 !!__default value
      else
        radius_fraction=max(0.001d0,radius_fraction)
      endif
    endif
    if (mod_valence) then
      print*,'Form factors will be modified for valence electrons of atoms with Z = 5 to Z = 9!'
      print*, '   with valence broadening (at. radii units) ',radius_fraction
    endif
  else if (rl(1:5)=='ARANG') then
    tai=trim(adjustl(rl(6:ll)))
    ltai=len_trim(tai)
    read(tai(1:ltai),*,iostat=iorr) tt_minmax,binwidth
    if (iorr/=0) then
      print*,"error reading 2theta min / max / step [degrees] from ARANG line - see below: "
      print'(a)',rl(1:ll)
      print*,'STOP'
      stop
    endif
    ext_ttminmax=.true.
  else if (rl(1:5)=='BROAD') then
    read(rl(6:ll),*)ubroad
  else if (rl(1:5)=='RRANG') then
    tai=trim(adjustl(rl(6:ll)))
    ltai=len_trim(tai)
    read(tai(1:ltai),*,iostat=iorr) r_minmaxstep
    if (iorr/=0) then
      print*,"error reading R min / max / step [Angstroem] from RRANG line - see below: "
      print'(a)',rl(1:ll)
      print*,'STOP'
      stop
    endif
    ext_rminstepmax=.true.
  else if (rl(1:5)=='DONOR') then
    do_normalize = .true.
    norm_mode=1
    tai=trim(adjustl(rl(6:ll)))
    ltai=len_trim(tai)
    if (ltai == 3) then
      if (tai(1:3) == 'fa2') then
        norm_mode=1
      else if (tai(1:3) == 'f2a') then
        norm_mode=2
      else if (tai(1:3) == 'Za2') then
        norm_mode=3
      else if (tai(1:3) == 'Z2a') then
        norm_mode=4
      endif
    endif
  else if (rl(1:5)=='SCALE') then
    tai=trim(adjustl(rl(6:ll)))
    ltai=len_trim(tai)
    flag_scal=tai(1:1)
    if (flag_scal == 'm') then
      read(tai(2:ltai),*) extscale
    else if (flag_scal == 't' .or. flag_scal == 'z') then
      read(tai(2:ltai),*,iostat=iotf) tail_fraction
      if (iotf/=0) then
        tail_fraction=0.05d0
      else
        tail_fraction=max(0.001d0,tail_fraction)
      endif
    endif
    !_______ Input new, AC15.08.2018
  else if (rl(1:5)=='PACKF') then
    read(rl(6:ll),*) Packing_Frac
  else if (rl(1:5)=='DOHAT'.or.rl(1:5)=='SOPER') then
    tai=trim(adjustl(rl(6:ll)))
    y_soper_tophat = tai(1:1)
    if (y_soper_tophat=='y'.or.(y_soper_tophat=='Y'.or.(y_soper_tophat=='1'))) then
      Soper_action=.true.
    else
      Soper_action=.false.
    endif
  else if (rl(1:5)=='QTHAT'.or.rl(1:5)=='QTSOP') then
    read(rl(6:ll),*) Soper_QT
  else if (rl(1:5)=='RFLAT') then
    read(rl(6:ll),*)Rmin_line
  else if (rl(1:5)=='NAVOL') then
    read(rl(6:ll),*) num_at_vol_den
  else if (rl(1:5)=='MDENS') then
    read(rl(6:ll),*) Mass_Density_gcm3
    given_mass_dens=.true.
!!!!!!!!!!!!!! kkzz?
  else if (rl(1:5)=='DEGBK') then
    read(rl(6:ll),*,iostat=iox)letter_bkg,r_degbk
    if (iox/=0) then
      letter_bkg='c'
      read(rl(6:ll),*)r_degbk
    endif
    if (letter_bkg=='c' .or. letter_bkg=='C') then
      type_bkg=1
      degbk=nint(r_degbk)
      bkg_lowr=0.d0
    else
      type_bkg=2
      bkg_lowr=r_degbk
    endif
  else if (rl(1:4)=='WLEN') then
    read(rl(5:ll),*)wlX
  else if (rl(1:5)=='NQCUT') then
    read(rl(6:ll),*)ncut_q, is_cut_rel
    allocate(qcuts(ncut_q))
  else if (rl(1:5)=='QCUTV') then
    if (.not.ALLOCATED(qcuts)) then
      print*,' ERROR: line QCUTV must be before line CUTV!'
      stop 'line QCUT must be before line NQCUT!'
    endif
    read(rl(6:ll),*,iostat=ioq)qcuts
    if (ioq/=0) then
      print'(a,i3,a)','Unable to read ',ncut_q,' values of Q cutoffs from line CUTV; stop'
      stop 'Unable to read values of Q cutoffs from line CUTV; stop'
    endif
  else if (rl(1:4)=='QMIN') then
    read(rl(5:ll),*)qcuts_min
  else if (rl(1:4)=='BINW') then
    read(rl(5:ll),*)binwidth
  else if (rl(1:5)=='INCOH') then
    read(rl(6:ll),*)Iuse_incoh
  else if (rl(1:4)=='NATO') then
    read(rl(5:ll),*)nato
    allocate(z_comp(nato),fractional_comp(nato))
    z_comp=1
  else if (rl(1:4)=='ZELE') then
    if (nato==0) then
      print*,' ERROR: line NATO must be before line ZELE!'
      stop 'line NATO must be before line ZELE!'
    endif
    read(rl(5:ll),*)z_comp
  else if (rl(1:4)=='CHEM') then
    if (nato==0) then
      print*,' ERROR: line ZELE must be before line CHEM!'
      stop 'line ZELE must be before line CHEM!'
    endif
    read(rl(5:ll),*)fractional_comp
    fractional_comp = fractional_comp/sum(fractional_comp)
    av_atwei=sum(fractional_comp*[(atwei(z_comp(iz)),iz=1,nato)])
    print'(a,1x,f16.6)','Computed average atomic weight                = ',av_atwei
    if (given_mass_dens) then
      num_at_vol_den = (Mass_Density_gcm3/av_atwei)*N_A*1.0d-24
      print'(a,1x,f16.6)','GIVEN mass density [g cm^-3]                  = ',Mass_Density_gcm3
      print'(a,1x,f16.6)','Computed atomic number density [Angstroem^-3] = ',num_at_vol_den
    else
      Mass_Density_gcm3 = num_at_vol_den*av_atwei*(1.0d24/N_A)
      print'(a,1x,f16.6)','Computed mass density [g cm^-3]               = ',Mass_Density_gcm3
      print'(a,1x,f16.6)','GIVEN atomic number density [Angstroem^-3]    = ',num_at_vol_den
    endif
  endif
enddo
1 close(iu11)

if (nfil==0) then 
  print*, 'ERROR: No PD pattern in the input file. The program stops!'
  STOP
endif


if (binwidth<sceps_DP) binwidth=-one

rmax=r_minmaxstep(2)
dr=r_minmaxstep(3)
rmin=r_minmaxstep(1)
npr=1+ceiling((rmax-rmin)/dr)
allocate(r(npr),sigGofR(npr,ncut_q),GofR(npr,ncut_q))
r = rmin+[(i*dr,i=0,npr-1)]

is_difference = ( ANY( CoeCol(1:2,2) < -sceps_DP ) .or. flag_scal=='n' )

ttmn= 99.99d99
ttmx=-99.99d99

do fk=1,nfil
  iu1=find_unit()
  smp=trim(adjustl((fn(fk))))
  open(iu1,status='old',action='read',file=trim(smp))
!finqi

  np = 0
  !!__using Qmin to cut, most restrictive bound
  ttcuts_min = asin(qcuts_min*wlX/(four*pi))/duet2r
  tt_minmax(1)=max(tt_minmax(1),ttcuts_min)
  do
    x_y_e = COLSREAD(nc=ncolfile,cols=col_index,colco=CoeCol,iunit=iu1,nr=n2r,ios=io2,cyc=icyc)
    if (io2 /= 0) exit
    if (icyc==1) cycle
    if (x_y_e(1)<tt_minmax(1).or.x_y_e(1)>tt_minmax(2)) cycle
    np = np + 1
    ttmx=max(ttmx,x_y_e(1))
    ttmn=min(ttmn,x_y_e(1))
  enddo
  if (is_cut_rel==1) then 
    Q_max_data=pi2*two*sin(ttmx*duet2r)/wlX
    qcuts=qcuts*Q_max_data
  endif
  is_cut_rel=0
  !!__using Qcuts to cut, most restrictive bound
  tt_minmax(2)=min( tt_minmax(2), asin(maxval(qcuts)*wlX/(four*pi))/duet2r )
  
  rewind(iu1)
  np = 0
  do
    x_y_e = COLSREAD(nc=ncolfile,cols=col_index,colco=CoeCol,iunit=iu1,nr=n2r,ios=io2,cyc=icyc)
    if (io2 /= 0) exit
    if (icyc==1) cycle
    if (x_y_e(1)<tt_minmax(1).or.x_y_e(1)>tt_minmax(2)) cycle
    np = np + 1
  enddo
  
  rewind(iu1)
  np1=np
  print*,'# of 2theta points:',np1
  if (allocated(i_expt)) deallocate(i_expt)
  if (allocated(ibroad)) deallocate(ibroad)
  if (allocated(sig)) deallocate(sig)
  if (allocated(tt)) deallocate(tt)
  allocate(i_expt(np1),tt(np1),sig(np1),ibroad(np1),scale_I_vec(np1))
  ibroad=one
  scale_I_vec = one
  np = 0
  
  do
    x_y_e = COLSREAD(nc=ncolfile,cols=col_index,colco=CoeCol,iunit=iu1,nr=n2r,ios=io2,cyc=icyc)
    if (io2 /= 0) exit
    if (icyc==1) cycle
    if (x_y_e(1)<tt_minmax(1).or.x_y_e(1)>tt_minmax(2)) cycle
    np = np + 1
    tt(np)=x_y_e(1); i_expt(np)=x_y_e(2); sig(np)=x_y_e(3)
  enddo
  close(iu1)
  print*,'Call:',z_comp,size(z_comp)
  print*,'# of 2theta points: ',np
  if (ubroad > sceps_DP) then
    xxx=pi2*ubroad*two/wlX
    ibroad = exp(-half * ( (xxx*sin(tt*duet2r))**2 ) )
  else 
    ibroad=one
  endif
  
!   if (iotf/=0 .and. size(qcuts)==1 .and. is_cut_rel==0) then
!     tail_fraction=((asin(qcuts*wlX/(four*pi))/duet2r - tt_minmax(1)) / binwidth) / np
!   endif

  if (binwidth>sceps_DP) then
  !______ UNIFORM 2THETA STEP, the step dx is a scalar
    if (maxval(abs(ibroad-one))>sceps_DP) then
      call GofR_from_start(wlen=wlX, rvals=r,Qmax_cut=qcuts, is_Qmax_rel=(is_cut_rel==1),Grvals=GofR,EGrvals=sigGofR, &
                         yarr=i_expt,e_yarr=sig,sc_mode=flag_scal, &
                         comp_X=fractional_comp,comp_Z=z_comp,rad_type=flag_rad,use_incoh=(Iuse_incoh==1), &
                         use_Z_to_scale_Xray=(norm_mode>2), &
                         xarr=tt,dx_scal=binwidth, twotheta1_q2_bigq3=1, scale_val=extscale, &
                         force_f2a=(modulo(norm_mode,2)==0),dont_subtract=is_difference,brobb=ibroad,&
                         tail_fr=tail_fraction,bkg_deg=degbk,bkg_lowr=bkg_lowr,beam_pol=pola)
    else
      call GofR_from_start(wlen=wlX, rvals=r,Qmax_cut=qcuts, is_Qmax_rel=(is_cut_rel==1),Grvals=GofR,EGrvals=sigGofR, &
                         yarr=i_expt,e_yarr=sig,sc_mode=flag_scal, &
                         comp_X=fractional_comp,comp_Z=z_comp,rad_type=flag_rad,use_incoh=(Iuse_incoh==1), &
                         use_Z_to_scale_Xray=(norm_mode>2), &
                         xarr=tt,dx_scal=binwidth, twotheta1_q2_bigq3=1, scale_val=extscale, &
                         force_f2a=(modulo(norm_mode,2)==0),dont_subtract=is_difference,&
                         tail_fr=tail_fraction,bkg_deg=degbk,bkg_lowr=bkg_lowr,beam_pol=pola)
    endif
  else
  !______ NON-UNIFORM 2THETA STEP, the step dx is a vector!! MISSING dx_vec
    if (maxval(abs(ibroad-one))>sceps_DP) then
      call GofR_from_start(wlen=wlX, rvals=r,Qmax_cut=qcuts, is_Qmax_rel=(is_cut_rel==1),Grvals=GofR,EGrvals=sigGofR, &
                         yarr=i_expt,e_yarr=sig,sc_mode=flag_scal, &
                         comp_X=fractional_comp,comp_Z=z_comp,rad_type=flag_rad,use_incoh=(Iuse_incoh==1), &
                         use_Z_to_scale_Xray=(norm_mode>2), &
                         xarr=tt, twotheta1_q2_bigq3=1, scale_val=extscale, &
                         force_f2a=(modulo(norm_mode,2)==0),dont_subtract=is_difference,brobb=ibroad,&
                         tail_fr=tail_fraction,bkg_deg=degbk,bkg_lowr=bkg_lowr,beam_pol=pola)
    else
      call GofR_from_start(wlen=wlX, rvals=r,Qmax_cut=qcuts, is_Qmax_rel=(is_cut_rel==1),Grvals=GofR,EGrvals=sigGofR, &
                         yarr=i_expt,e_yarr=sig,sc_mode=flag_scal, &
                         comp_X=fractional_comp,comp_Z=z_comp,rad_type=flag_rad,use_incoh=(Iuse_incoh==1), &
                         use_Z_to_scale_Xray=(norm_mode>2), &
                         xarr=tt, twotheta1_q2_bigq3=1, scale_val=extscale, &
                         force_f2a=(modulo(norm_mode,2)==0),dont_subtract=is_difference,&
                         tail_fr=tail_fraction,bkg_deg=degbk,bkg_lowr=bkg_lowr,beam_pol=pola)
    endif
  endif

  iu2=find_unit()
  idot=index(smp(1:lfn(fk)),'.',.true.)
  if (idot==0) then
    idot=lfn(fk)+1
  endif
  open(iu2,status='replace',file=smp(1:idot-1)//'.rpdfn')
  print*,'npr = ',npr
  do jr=1,npr
    write(iu2,*)r(jr),GofR(jr,:),sigGofR(jr,:)
  enddo
  close(iu2)
  
  scale_I_vec(:) = one * scale_factor_I
  print*, 'MK_G_OF_R, scale_factor_I ', scale_factor_I, scale_I_vec(1), scale_I_vec(2)
  
  iu4=find_unit()
  open(iu4,status='replace',file=smp(1:idot-1)//'.sqobs')
  do iq=1,np1
!     write(iu4,*)Qvec(iq),dQvec(iq)*SofQ_arr(iq),Qvec(iq)*dQvec(iq)*(SofQ_arr(iq)-1),i_expt(iq),&
!           dQvec(iq)*Bsub_bkg(iq),dQvec(iq),fofa2_inc(iq),ff_AverageSquared(iq),ff_SquaredAverage(iq)
    write(iu4,*)Qvec(iq),SofQ_arr(iq),Qvec(iq)*(SofQ_arr(iq)-1),i_expt(iq),scale_I_vec(iq),&
          Bsub_bkg(iq),dQvec(iq),fofa2_inc(iq),ff_AverageSquared(iq),ff_SquaredAverage(iq)
  enddo
  close(iu4)

enddo
! dq = 0.00052536373507057542d0
! 
! npmq = 2+ceiling((Qvec(np1) )/dq)
! allocate(Iofq(npmq,2),essex(npmq),aux(npmq))
! essex = .false.
! 
! IofQ=zero
! print*,dq,sum(dQvec(1:np1))/np1,minval(dQvec(1:np1)),maxval(dQvec(1:np1))
! do iq=1,np1
!   x = Qvec(iq)
!   y = dQvec(iq)*SofQ_arr(iq)
!   z = one
!   q1 = x - dQvec(iq)/2.0d0
!   q2 = x + dQvec(iq)/2.0d0
!   yy = y*dq/(q2-q1)
!   zz = z*dq/(q2-q1)
!   ozz=one/(zz**2)
!   if (ozz<eps_DP.or.abs(ozz)>one/eps_DP.or.isnan(ozz)) &
!     print*,'error ', z,zz,ozz,dq,q2-q1
!   ii1=floor(q1/dq)-1
!   ii2=ceiling(q2/dq)+1
!   do ii=ii1,ii2
!     beta = max(0.d0, min(q2,dq*ii)-max(q1,dq*(ii-1)))/dq
!     if (beta<eps_dp) cycle
!     IofQ(ii,1) = IofQ(ii,1)+beta*ozz
!     IofQ(ii,2) = IofQ(ii,2)+beta*yy*ozz
!     essex(ii) = .true.
!   enddo
! enddo
! 
! aux = 0.d0
! where(essex(:))
!   aux(:)=IofQ(:,2)/IofQ(:,1)
!   IofQ(:,2) = 1.d0/sqrt(IofQ(:,1))
!   IofQ(:,1) = aux(:)
! end where
! open(iu4,status='replace',file=smp(1:idot-1)//'.sqobs_rebin')
! do iw=1,npmq
!   if (essex(iw)) write(iu4,'(1x,f15.6,2(1x,f15.8))') iw*dq, IofQ(iw,:)
! enddo
! close(iu4)



print*, '******* JOB PAIR CORRELATION DISTRIBUTION DONE! *******'

end program dopdf1
