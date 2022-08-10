program SCH_3D_ArbitraryQW_2
    implicit none
!*****************************************************************************************
!***************************************** Variables *************************************
!   ---------------Constant---------------------------------------------------------------
    double precision :: muo,pi,epso,co,qe,h,h_rad,m0
!   ---------------Domain and Iteration---------------------------------------------------
    integer          :: i,m,n,o,nmax,EN,S_max_itr,NFFT,NF,iMov
    integer          :: ie,ib,je,jb,ke,kb
    integer          :: iAir1_a,iAir1_o,iAir2_a,iAir2_o
    integer          :: jAir1_a,jAir1_o,jAir2_a,jAir2_o
    integer          :: kAir1_a,kAir1_o,kAir2_a,kAir2_o
    integer          :: jCirc
    integer          :: slice_x,slice_y,slice_z

    double precision :: dt,ALAMBDA,x0,xp,y0,yp,z0,zp
    double precision :: xx,ddx,yy,ddy,zz,ddz
    double precision :: x,y,z

    double precision :: tx_GaAs,ty_GaAs,tz_GaAs
    double precision :: t_Air_top,t_Air_bot,t_Air_lef,t_Air_rig,t_Air_fro,t_Air_bac
!   ---------------PML--------------------------------------------------------------------
    integer          :: pml_grid,minv
    double precision :: sigma0
    complex          :: Rpml
!   ---------------Material Properties----------------------------------------------------
    double precision :: me_GaAs,eps_GaAs,Eg_GaAs,X_GaAs,Ec_GaAs
    real             :: mat_PML,mat_Air,mat_GaAs
    real             :: rad_Circ
    double precision :: dist2p
!   ---------------Schrodinger Equ.-------------------------------------------------------
    double precision :: ra,psi_norm,phi_norm,Fs,arg
!   ---------------Test Function----------------------------------------------------------
    double precision :: Lx2,Ly2,Lz2

!*****************************************************************************************
!***************************************** Universal Constants ***************************
    parameter (pi=4.0*atan(1.0),        muo=4.0*pi*1.0d-7,       epso=8.854187817d-12)
    parameter (co=1.0/dsqrt(muo*epso),  qe=1.602176487d-19,      m0=9.10938215d-31)
    parameter (h_rad=1.054571628d-34,   h=6.62606896d-34)

!******************************************** Material Parameters ************************
!   ---------------GaAs-------------------------------------------------------------------
    parameter (me_GaAs = 0.067*m0, eps_GaAs = 15.0)
    parameter (Eg_GaAs = 1.90*qe, X_GaAs = 11.1*qe, Ec_GaAs = -X_GaAs)

!******************************************** Computational Domain ***********************
!   ---------------Iteration--------------------------------------------------------------
    parameter (nmax=1, iMov=8)
    parameter (NF=2**6)
    parameter (NFFT=18, S_max_itr=2**NFFT, EN=5)
    parameter (pml_grid = 0)
!   ---------------Domain-----------------------------------------------------------------
    parameter (tx_GaAs = 8.0d-9,   ty_GaAs = 13.0d-9,     tz_GaAs = 3.5d-9)
!    parameter (tx_GaAs = 25.0d-9,   ty_GaAs = 35.0d-9,     tz_GaAs = 20d-9)

    parameter (t_Air_top = 20.0d-9, t_Air_bot = 20.0d-9)
    parameter (t_Air_lef = 20.0d-9, t_Air_rig = 20.0d-9)
    parameter (t_Air_fro = 20.0d-9, t_Air_bac = 20.0d-9)

    parameter (xx = tx_GaAs, ddx = 0.25d-9)
    parameter (yy = ty_GaAs, ddy = 0.25d-9)
    parameter (zz = tz_GaAs, ddz = 0.25d-9)

    parameter (ie = nint((xx + t_Air_fro + t_Air_bac)/ddx) + 2*pml_grid, ib = ie + 1)
    parameter (je = nint((yy + t_Air_lef + t_Air_rig)/ddy) + 2*pml_grid, jb = je + 1)
    parameter (ke = nint((zz + t_Air_top + t_Air_bot)/ddz) + 2*pml_grid, kb = ke + 1)

    ! X AXIS
    parameter (iAir1_a = pml_grid + 1,              iAir1_o = iAir1_a + nint(t_Air_bac/ddx))
    parameter (iAir2_a = iAir1_o + nint(xx/ddx),    iAir2_o = iAir2_a + nint(t_Air_fro/ddx))
    ! Y AXIS
    parameter (jAir1_a = pml_grid + 1,              jAir1_o = jAir1_a + nint(t_Air_lef/ddy))
    parameter (jAir2_a = jAir1_o + nint(yy/ddy),    jAir2_o = jAir2_a + nint(t_Air_rig/ddy))
    parameter (jCirc = jAir1_o + nint(yy/(2*ddy)))
    ! Z AXIS
    parameter (kAir1_a = pml_grid + 1,              kAir1_o = kAir1_a + nint(t_Air_bot/ddz))
    parameter (kAir2_a = kAir1_o + nint(zz/ddz),    kAir2_o = kAir2_a + nint(t_Air_top/ddz))

!   ---------------Initial Wave-----------------------------------------------------------
    parameter (xp = 0.06*8.0d-9, x0=5.0d-9)
    parameter (yp = 0.06*13.0d-9, y0=5.0d-9)
    parameter (zp = 0.06*3.0d-9, z0=5.0d-9)
    parameter (ALAMBDA=20d-9)

!   ---------------PML--------------------------------------------------------------------
    parameter (sigma0=0.005, Rpml=cmplx(cos(pi/4),sin(pi/4)))

!*****************************************************************************************
!*************************************************** Vectors *****************************
    double precision, dimension(S_max_itr)  :: Han_Win,psirt,psimt,f
    double precision, dimension(EN)         :: eig_energy,eig_freq
    integer, dimension(3)                   :: Probe

    double precision, allocatable :: Ec(:,:,:),U(:,:,:),me(:,:,:),eps(:,:,:),mat(:,:,:)
    double precision, allocatable :: psir(:,:,:),psim(:,:,:),phir(:,:,:),phim(:,:,:),phirf(:,:,:),phimf(:,:,:)
    double precision, allocatable :: ca_psim(:,:,:),cb_psim(:,:,:),cc_psim(:,:,:),ca_psir(:,:,:),cb_psir(:,:,:),cc_psir(:,:,:)
    double precision, allocatable :: sigmax(:,:,:),gamma_re(:,:,:),gamma_im(:,:,:),delta_psir(:,:,:),delta_psim(:,:,:)
    double precision, allocatable :: psir_init(:,:,:)

    allocate(Ec(ib,jb,kb),U(ib,jb,kb),me(ib,jb,kb),eps(ib,jb,kb),mat(ib,jb,kb), &
        psir(ib,jb,kb),psim(ib,jb,kb),phir(ib,jb,kb),phim(ib,jb,kb),phirf(ib,jb,kb),phimf(ib,jb,kb), &
        ca_psim(ib,jb,kb),cb_psim(ib,jb,kb),cc_psim(ib,jb,kb),ca_psir(ib,jb,kb),cb_psir(ib,jb,kb),cc_psir(ib,jb,kb), &
        sigmax(ib,jb,kb),gamma_re(ib,jb,kb),gamma_im(ib,jb,kb),delta_psir(ib,jb,kb),delta_psim(ib,jb,kb), &
        psir_init(ib,jb,kb))

!**************************************************** Open Files *************************
    open(21,file = 'SCH_3D_ArbitQW3_param',status='replace')
    open(22,file = 'SCH_3D_ArbitQW3_U',status='replace')
    open(23,file = 'SCH_3D_ArbitQW3_phirf',status = 'replace')
    open(24,file = 'SCH_3D_ArbitQW3_phimf',status = 'replace')
    open(25,file = 'SCH_3D_ArbitQW3_phir',status = 'replace')

    open(26,file = 'SCH_3D_ArbitQW3_psirt',status = 'replace')
    open(27,file = 'SCH_3D_ArbitQW3_psimt',status = 'replace')
    open(28,file = 'SCH_3D_ArbitQW3_psirtf',status = 'replace')
    open(29,file = 'SCH_3D_ArbitQW3_psir_init',status = 'replace')

    open(91,file = 'SCH_3D_ArbitQW3_psir_movX',access = 'stream', form = 'unformatted',status = 'replace')
    open(911,file = 'SCH_3D_ArbitQW3_psir_movY',access = 'stream', form = 'unformatted',status = 'replace')
    open(912,file = 'SCH_3D_ArbitQW3_psir_movZ',access = 'stream', form = 'unformatted',status = 'replace')

    open(92,file = 'SCH_3D_ArbitQW3_psim_mov',access = 'stream', form = 'unformatted',status = 'replace')

!*****************************************************************************************
!**************************************************** Variables **************************
    me = m0;    Ec = 0;    eps = epso;
    U = 0;
    rad_circ = 2.5e-9

    mat_PML = 1.0; mat_Air = 0.50; mat_GaAs = 0.0;
    mat = mat_PML;
    slice_X = int(ie/2); slice_Y = int(je/2); slice_Z = int(ke/2);

    dt = 0.2d-17
    ra = (0.5*h_rad/me_GaAs)*(dt/ddx**2)
    Fs = 1/dt

    write(21,*)ib,jb,kb,ddx,ddy,ddz,nmax,dt,slice_X,slice_Y,slice_Z,iMov,S_max_itr

!*********************************************** PML Initialization  **********************************************
    sigmax = 0; gamma_re = 0; gamma_im = 0;
!    do m = 1,pml_grid
!        minv = ib - m
!        sigmax(m,:,:) = sigma0*(float(m - pml_grid))**2
!        sigmax(minv,:,:) = sigma0*(float(m - pml_grid))**2
!
!        minv = jb - m
!        sigmax(:,m,:) = sigma0*(float(m - pml_grid))**2
!        sigmax(:,minv,:) = sigma0*(float(m - pml_grid))**2
!
!        minv = kb - m
!        sigmax(:,:,m) = sigma0*(float(m - pml_grid))**2
!        sigmax(:,:,minv) = sigma0*(float(m - pml_grid))**2
!    enddo
!
!    gamma_re = real((1/(1 + Rpml*sigmax))**2)
!    gamma_im = aimag((1/(1 + Rpml*sigmax))**2)
!
!    open(51,file = 'SCH_3D_ArbitQW3_sigmax',status = 'replace')
!    write(51,*)sigmax(:,:,ke/2);     close(51);
!    open(52,file = 'SCH_3D_ArbitQW3_sigmax1',status = 'replace')
!    write(52,*)sigmax(:,:,ke/2-10);     close(52);
!    open(53,file = 'SCH_3D_ArbitQW3_sigmax2',status = 'replace')
!    write(53,*)sigmax(:,:,ke/2+10);     close(53);
!
!    open(61,file = 'SCH_3D_ArbitQW3_gamma_re',status = 'replace')
!    write(61,*)gamma_re(:,:,ke/2);     close(61);
!    open(62,file = 'SCH_3D_ArbitQW3_gamma_im',status = 'replace')
!    write(62,*)gamma_im(:,:,ke/2);     close(62);
!
!    open(63,file = 'SCH_3D_ArbitQW3_gamma_re1',status = 'replace')
!    write(63,*)gamma_re(:,:,ke/2-10);   close(63);
!    open(64,file = 'SCH_3D_ArbitQW3_gamma_im1',status = 'replace')
!    write(64,*)gamma_im(:,:,ke/2-10);   close(64);

!*********************************************** Material Initialization  ****************
     do m = 1,jb
        if(m.ge.jAir1_a.and.m.le.jAir1_o)then
            me (iAir1_a:iAir2_o,m,kAir1_a:kAir2_o) = m0
            mat(iAir1_a:iAir2_o,m,kAir1_a:kAir2_o) = mat_Air
            U  (iAir1_a:iAir2_o,m,kAir1_a:kAir2_o) = 0

        elseif(m.gt.jAir1_o.and.m.lt.jAir2_a)then
            me (iAir1_a:iAir1_o,m,kAir1_a:kAir2_o) = m0
            mat(iAir1_a:iAir1_o,m,kAir1_a:kAir2_o) = mat_Air
            U  (iAir1_a:iAir1_o,m,kAir1_a:kAir2_o) = 0

            mat(iAir2_a:iAir2_o,m,kAir1_a:kAir2_o) = m0
            mat(iAir2_a:iAir2_o,m,kAir1_a:kAir2_o) = mat_Air
            U  (iAir2_a:iAir2_o,m,kAir1_a:kAir2_o) = 0

            me (iAir1_a:iAir2_o,m,kAir1_a:kAir1_o) = m0
            mat(iAir1_a:iAir2_o,m,kAir1_a:kAir1_o) = mat_Air
            U  (iAir1_a:iAir2_o,m,kAir1_a:kAir1_o) = 0

            me (iAir1_a:iAir2_o,m,kAir2_a:kAir2_o) = m0
            mat(iAir1_a:iAir2_o,m,kAir2_a:kAir2_o) = mat_Air
            U  (iAir1_a:iAir2_o,m,kAir2_a:kAir2_o) = 0

            me (iAir1_o+1:iAir2_a-1,m,kAir1_o+1:kAir2_a-1) = me_GaAs
            Ec (iAir1_o+1:iAir2_a-1,m,kAir1_o+1:kAir2_a-1) = -X_GaAs
            mat(iAir1_o+1:iAir2_a-1,m,kAir1_o+1:kAir2_a-1) = mat_GaAs

            U  (iAir1_o+1:iAir2_a-1,m,kAir1_o+1:kAir2_a-1) = 0
            U  (iAir1_o+1:iAir2_a-1,m,kAir1_o+1) = 5*qe
            U  (iAir1_o+1:iAir2_a-1,m,kAir1_o+1) = 5*qe
            U  (iAir1_o+1,m,kAir1_o+1:kAir2_a-1) = 5*qe
            U  (iAir2_a-1,m,kAir1_o+1:kAir2_a-1) = 5*qe
            U  (iAir1_o+1:iAir2_a-1,jAir1_o+1,kAir1_o+1:kAir2_a-1) = 5*qe
            U  (iAir1_o+1:iAir2_a-1,jAir2_a-1,kAir1_o+1:kAir2_a-1) = 5*qe


            if (m.gt.(jCirc - int(rad_Circ/ddx)).and.m.lt.(jCirc + int(rad_Circ/ddx))) then
                do n = iAir1_o,iAir1_o + nint(rad_Circ/ddx)
                    dist2p = sqrt((float(n) - float(iAir1_o))**2 + (float(m) - float(jCirc))**2)
                    if (dist2p.le.(rad_Circ/ddx))then
                        i = i+1
                        me (n,m,kAir1_a:kAir2_o) = m0
                        mat(n,m,kAir1_a:kAir2_o) = mat_Air
                        U  (n,m,kAir1_a:kAir2_o) = 0
                    endif
                enddo
                do n = iAir2_a - nint(rad_Circ/ddx), iAir2_a
                    dist2p = sqrt((float(n) - float(iAir2_a))**2 + (float(m) - float(jCirc))**2)
                    if (dist2p.le.(rad_Circ/ddx))then
                        me (n,m,kAir1_a:kAir2_o) = m0
                        mat(n,m,kAir1_a:kAir2_o) = mat_Air
                        U  (n,m,kAir1_a:kAir2_o) = 0
                    endif
                enddo
            endif
        elseif(m.ge.jAir2_a.and.m.le.jAir2_o)then
            me (iAir1_a:iAir2_o,m,kAir1_a:kAir2_o) = m0
            mat(iAir1_a:iAir2_o,m,kAir1_a:kAir2_o) = mat_Air
            U  (iAir1_a:iAir2_o,m,kAir1_a:kAir2_o) = 0
        endif
    enddo

    open(41,file = 'SCH_3D_ArbitQW3_Mat_YZ', access = 'stream',FORM = 'unformatted')
    write(41)mat(slice_X,1:jb,1:kb);     close(41);

    open(42,file = 'SCH_3D_ArbitQW3_Mat_XZ', access = 'stream',FORM = 'unformatted')
    write(42)mat(1:ib,slice_Y,1:kb);     close(42);
    open(420,file = 'SCH_3D_ArbitQW3_Mat_XZ2', access = 'stream',FORM = 'unformatted')
    write(420)mat(1:ib,slice_Y+40,1:kb);     close(420);
    open(421,file = 'SCH_3D_ArbitQW3_Mat_XZ3', access = 'stream',FORM = 'unformatted')
    write(421)mat(1:ib,slice_Y-40,1:kb);     close(421);

    open(43,file = 'SCH_3D_ArbitQW3_Mat_XY', access = 'stream',FORM = 'unformatted')
    write(43)mat(1:ib,1:jb,slice_Z);     close(43);
    open(430,file = 'SCH_3D_ArbitQW3_Mat_XY2', access = 'stream',FORM = 'unformatted')
    write(430)mat(1:ib,1:jb,slice_Z-20);     close(430);

    open(44,file = 'SCH_3D_ArbitQW3_Mat', access = 'stream',FORM = 'unformatted')
    write(44)mat(1:ib,1:jb,1:kb);     close(44);

!    open(45,file = 'SCH_3D_ArbitQW_Mat2',status = 'replace')
!    write(45,*)mat;     close(45);

!*********************************************** Schrodinger Initialization  *************
    ca_psir = 1.0; cb_psir = -h_rad*dt/(2*me*ddx**2); cc_psir =  dt/h_rad;
    ca_psim = 1.0; cb_psim =  h_rad*dt/(2*me*ddx**2); cc_psim = -dt/h_rad;
    psir_init = 0.0; phirf = 0.0 ; phimf = 0.0

!*********************************************** Source **********************************
    Lx2 = ie/2*ddx; Ly2 = je/2*ddy; Lz2 = ke/2*ddz;
    Probe = (/int(Lx2/ddx),int(Ly2/ddy),int(Lz2/ddz)/)
    print*,'Probe x',Probe(1),Probe(1)*ddx/1.0d-9,Lx2
    print*,'Probe y',Probe(2),Probe(2)*ddy/1.0d-9,Ly2
    print*,'Probe z',Probe(3),Probe(3)*ddz/1.0d-9,Lz2

    do m = 1,ib
        x = float(m)*ddx
        do n = 1,jb
            y = float(n)*ddy
            do o = 1,kb
                z = float(o)*ddz
                psir_init(m,n,o) = exp(-((x-Lx2)/xp)**2 - ((y-Ly2)/yp)**2 -((z-Lz2)/zp)**2)
            enddo
        enddo
    enddo

!    do o = 1,kb
!        z = float(o)*ddz
!        psir_init(int(ie/2),int(je/2),o) = exp(-((z-Lz2)/zp)**2)
!    enddo

    psir = 0.0 ; psim = 0.0;
    psir = psir_init
    psi_norm = sqrt(sum(psir**2 + psim**2))
    psir = psir/psi_norm;
    psim = psim/psi_norm;
    psi_norm = sqrt(sum(psir**2 + psim**2))
    print*,'Psi Norm: ',psi_norm
    write(29,*)psir

!*********************************************** Hanning Window **************************
    do m = 1,S_max_itr
        Han_Win(m) = .5*(1-cos(2*pi*m/S_max_itr))
    enddo
    open(31,file = 'SCH_3D_ArbitQW3_HanWin',status = 'replace')
    write(31,*)Han_Win;     close(31);

!*********************************************** Frequency Axis ************************
    do m = 1,S_max_itr
        f(m) = Fs*float(m-1)/S_max_itr - Fs/2.0
    enddo
    write(28,*)f*h/qe

!*********************************************** Potential Energy ************************
!    has been defined in the material initialization
!    U = 0;

!!*****************************************************************************************
!!*********************************************** Printing ******************************
    print*, '------------------------------------------------------'
    print*, 'dt ',dt
    print*, 'ra ',ra
    print*, '------------------------------------------------------'
    print*, 'Ie:Je:Ke ',ie,je,ke, ' points'
    print*, 'Ib:Jb:Kb ',ib,jb,kb, ' points'
    print*, '------------------------------------------------------'
    print*, 'x : ',xx,  ' m'
    print*, 'y : ',yy,  ' m'
    print*, 'z : ',zz,  ' m'
    print*, '------------------------------------------------------'
    print*, 'dx: ',ddx, ' m'
    print*, 'dy: ',ddy, ' m'
    print*, 'dz: ',ddz, ' m'
    print*, '------------------------------------------------------'
    print*, 'iAir1_a:o ',iAir1_a,iAir1_o
    print*, 'iAir2_a:o ',iAir2_a,iAir2_o
    print*, '------------------------------------------------------'
    print*, 'jAir1_a:o ',jAir1_a,jAir1_o
    print*, 'jAir2_a:o ',jAir2_a,jAir2_o
    print*, '------------------------------------------------------'
    print*, 'kAir1_a:o ',kAir1_a,kAir1_o
    print*, 'kAir2_a:o ',kAir2_a,kAir2_o
    print*, '------------------------------------------------------'
    print*, 'ca:cb:cc',ca_psir(ie/2,je/2,ke/2),cb_psir(ie/2,je/2,ke/2),cc_psir(ie/2,je/2,ke/2)
    print*, 'ca:cb:cc',ca_psim(ie/2,je/2,ke/2),cb_psim(ie/2,je/2,ke/2),cc_psim(ie/2,je/2,ke/2)
    print*, '------------------------------------------------------'

!*****************************************************************************************
!*********************************************** Main Loop *******************************
    do m=1,S_max_itr
        if(mod(m,32).eq.0) write(6,*)m

!*****************************************************************************************

!        delta_psim(2:ie,2:je,2:ke) = -6*psim(2:ie,2:je,2:ke)   + &
!                                        psim(1:ie-1,2:je,2:ke) + psim(3:ib,2:je,2:ke) + &
!                                        psim(2:ie,1:je-1,2:ke) + psim(2:ie,3:jb,2:ke) + &
!                                        psim(2:ie,2:je,1:ke-1) + psim(2:ie,2:je,3:kb)
!        delta_psir(2:ie,2:je,2:ke) = -6*psir(2:ie,2:je,2:ke)   + &
!                                        psir(1:ie-1,2:je,2:ke) + psir(3:ib,2:je,2:ke) + &
!                                        psir(2:ie,1:je-1,2:ke) + psir(2:ie,3:jb,2:ke) + &
!                                        psir(2:ie,2:je,1:ke-1) + psir(2:ie,2:je,3:kb)
!
!        psir(2:ie,2:je,2:ke) = psir(2:ie,2:je,2:ke) + &
!                               cb_psir(2:ie,2:je,2:ke)*(delta_psim(2:ie,2:je,2:ke)*gamma_re(2:ie,2:je,2:ke) + &
!                                                        delta_psir(2:ie,2:je,2:ke)*gamma_im(2:ie,2:je,2:ke))+ &
!                               cc_psir(2:ie,2:je,2:ke)*U(2:ie,2:je,2:ke)*psim(2:ie,2:je,2:ke)
!
!        psim(2:ie,2:je,2:ke) = psim(2:ie,2:je,2:ke) + &
!                               cb_psim(2:ie,2:je,2:ke)*(delta_psir(2:ie,2:je,2:ke)*gamma_re(2:ie,2:je,2:ke) - &
!                                                        delta_psim(2:ie,2:je,2:ke)*gamma_im(2:ie,2:je,2:ke))+ &
!                               cc_psim(2:ie,2:je,2:ke)*U(2:ie,2:je,2:ke)*psir(2:ie,2:je,2:ke)

!*****************************************************************************************

        psir(2:ie,2:je,2:ke) = psir(2:ie,2:je,2:ke) + &
                               cb_psir(2:ie,2:je,2:ke)*(-6 * psim(2:ie,2:je,2:ke) + &
                                   psim(1:ie-1,2:je,2:ke) + psim(3:ib,2:je,2:ke) + &
                                   psim(2:ie,1:je-1,2:ke) + psim(2:ie,3:jb,2:ke) + &
                                   psim(2:ie,2:je,1:ke-1) + psim(2:ie,2:je,3:kb))+ &
                               cc_psir(2:ie,2:je,2:ke)*U(2:ie,2:je,2:ke)*psim(2:ie,2:je,2:ke)

        psim(2:ie,2:je,2:ke) = psim(2:ie,2:je,2:ke) + &
                               cb_psim(2:ie,2:je,2:ke)*(-6 * psir(2:ie,2:je,2:ke) + &
                                   psir(1:ie-1,2:je,2:ke) + psir(3:ib,2:je,2:ke) + &
                                   psir(2:ie,1:je-1,2:ke) + psir(2:ie,3:jb,2:ke) + &
                                   psir(2:ie,2:je,1:ke-1) + psir(2:ie,2:je,3:kb))+ &
                               cc_psim(2:ie,2:je,2:ke)*U(2:ie,2:je,2:ke)*psir(2:ie,2:je,2:ke)

!*****************************************************************************************

        psirt(m) =  psir(Probe(1),Probe(2),Probe(3))*Han_Win(m)
        psimt(m) =  psim(Probe(1),Probe(2),Probe(3))*Han_Win(m)

        write(26,*)psirt(m)
        write(27,*)psimt(m)

        if(mod(m,iMov).eq.0) write(91 )psir(Probe(1),:,:);
        if(mod(m,iMov).eq.0) write(911)psir(:,Probe(2),:);
        if(mod(m,iMov).eq.0) write(912)psir(:,:,Probe(3));
    enddo

!**************************************************** Eigen Energies *********************
    call FFTC(psirt,psimt,S_max_itr,NFFT)
    psirt = dsqrt(psirt**2 + psimt**2)
    write(28,*)psirt;           close(28);
    eig_freq(1) = 2*pi*f(maxloc(psirt,1))
    eig_energy(1) = abs(h_rad * eig_freq(1))
    arg = eig_freq(1) * dt

    print*,'eigen energy', 1 ,'=' ,eig_energy(1)/qe * 1000.0 , 'meV'
    open(101,file = 'SCH_3D_ArbitQW3_EigEnergy',status = 'replace')
    write(101,*)eig_energy;     close(101);

!**************************************************** Eigen Function *********************
    psir = 0.0 ; psim = 0.0; phir = 0.0 ; phim = 0.0
    psir = psir_init

    psi_norm = sqrt(sum(psir**2 + psim**2))
    psir = psir/psi_norm;
    psim = psim/psi_norm;
    psi_norm = sqrt(sum(psir**2 + psim**2))

    print*,'Psi Norm: ',psi_norm

    do m=1,S_max_itr
        psir(2:ie,2:je,2:ke) = psir(2:ie,2:je,2:ke) + &
                               cb_psir(2:ie,2:je,2:ke)*(-6 * psim(2:ie,2:je,2:ke) + &
                                   psim(1:ie-1,2:je,2:ke) + psim(3:ib,2:je,2:ke) + &
                                   psim(2:ie,1:je-1,2:ke) + psim(2:ie,3:jb,2:ke) + &
                                   psim(2:ie,2:je,1:ke-1) + psim(2:ie,2:je,3:kb))+ &
                               cc_psir(2:ie,2:je,2:ke)*U(2:ie,2:je,2:ke)*psim(2:ie,2:je,2:ke)

        psim(2:ie,2:je,2:ke) = psim(2:ie,2:je,2:ke) + &
                               cb_psim(2:ie,2:je,2:ke)*(-6 * psir(2:ie,2:je,2:ke) + &
                                   psir(1:ie-1,2:je,2:ke) + psir(3:ib,2:je,2:ke) + &
                                   psir(2:ie,1:je-1,2:ke) + psir(2:ie,3:jb,2:ke) + &
                                   psir(2:ie,2:je,1:ke-1) + psir(2:ie,2:je,3:kb))+ &
                               cc_psim(2:ie,2:je,2:ke)*U(2:ie,2:je,2:ke)*psir(2:ie,2:je,2:ke)

        phir(2:ie,2:je,2:ke) = phir(2:ie,2:je,2:ke) + &
                               Han_Win(m)*(psir(2:ie,2:je,2:ke)*cos(arg*m) + psim(2:ie,2:je,2:ke)*sin(arg*m));
        phim(2:ie,2:je,2:ke) = phim(2:ie,2:je,2:ke) + &
                               Han_Win(m)*(psim(2:ie,2:je,2:ke)*cos(arg*m) - psir(2:ie,2:je,2:ke)*sin(arg*m));

        if(mod(m,NF).eq.0) write(25,*)phir
    enddo
    phi_norm = sqrt(sum(phir**2 + phim**2))
    phir = phir/phi_norm;
    phim = phim/phi_norm;
    phi_norm = sqrt(sum(phir**2 + phim**2))
    print*, 'PHI = ', phi_norm

    phirf = phir
    phimf = phim

    write(23,*)phirf
    write(24,*)phimf

!*****************************************************************************************
!*********************************************** The End *********************************
    print*, '-------------Alhamdulillaah-------------'

end program SCH_3D_ArbitraryQW_2

!*****************************************************************************************
!*********************************************** FFT *************************************
SUBROUTINE FFTC (AR,AI,N,M)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: N,M
    INTEGER :: N1,N2,I,J,K,L,L1,L2
    double precision :: PI,A1,A2,Q,U,V
    double precision, INTENT (INOUT), DIMENSION (N) :: AR,AI

    PI = 4.0*ATAN(1.0)
    N2 = N/2
    N1 = 2**M
    IF(N1.NE.N) STOP 'Indices do not match'
!   --------------------------------------------------------------------------------------
!   Rearrange the data to the bit reversed order
    L = 1
    DO K = 1, N-1
        IF (K.LT.L) THEN
            A1    = AR(L)
            A2    = AI(L)
            AR(L) = AR(K)
            AR(K) = A1
            AI(L) = AI(K)
            AI(K) = A2
        END IF
        J   = N2
        DO WHILE (J.LT.L)
            L = L-J
            J = J/2
        END DO
        L = L+J
    END DO
!   --------------------------------------------------------------------------------------
!   Perform additions at all levels with reordered data
    L2 = 1
    DO L = 1, M
        Q  =  0.0
        L1 =  L2
        L2 =  2*L1
        DO K = 1, L1
            U   =  COS(Q)
            V   = -SIN(Q)
            Q   =  Q + PI/L1
            DO J = K, N, L2
                I     =  J + L1
                A1    =  AR(I)*U-AI(I)*V
                A2    =  AR(I)*V+AI(I)*U
                AR(I) =  AR(J)-A1
                AR(J) =  AR(J)+A1
                AI(I) =  AI(J)-A2
                AI(J) =  AI(J)+A2
            END DO
        END DO
    END DO
!   --------------------------------------------------------------------------------------
    AR = cshift(AR,-N2)
    AI = cshift(AI,-N2)
END SUBROUTINE FFTC
