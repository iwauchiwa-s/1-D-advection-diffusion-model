      program moldif
!============================================================================
!-----------------------------------------------------------------------
!  	 MOLDIF : 1-D molecular and eddy diffusion model
!           by S. S.
!-----------------------------------------------------------------------
!
!     Revision history : 
!
!     Feb. 11, 2016 : start developing
!     Feb. 15, 2016 : basic model developed
!     Feb. 22, 2016 : modified for sim. type change
!
!-----------------------------------------------------------------------
!     constants, parameters and variables used in the program
!     ----------------------------------
!  vertical structure
!     zdt          [m]   height of model top
!     dzd1         [m]   default thickness of layer
!     nzd          []    number of vertical layer
!     zdb(nzd+1)   [m]   height at layer boundary
!     zdm(nzd)     [m]   height at mid-point of layer
!     dzd(nzd)     [m]   thickness of layer
!     dzm(nzd-1)   [m]   thickness between mid-points of layer
!  constants
!     gg         [m/s2]     gravitational acceleration
!     rg         [J/mol/K]  gas constant
!     amsa       [kg/mol]   molar mass of air
!     amsg       [kg/mol]   molar mass of target gas
!     dco2       [m2/s]     molecular diffusion coefficient of co2
!     sfd        []         scale factor of diffusion coefficient for target gas
!     tdf        []         thermal diffusion factor of target gas
!     cnp        []         weighting parameter of Crank-Nicolson Method (=0.5)
!  given variables
!     temp(nzd)    [K]       air temperature
!     tmpd(nzd)    [K/m]     vertical temperature gradient, dT/dz
!     na(nzd)      [mol/m3]  total air quantity
!     dmol(nzd)    [m2/s]    molecular diffusion coefficient of target gas
!     dedy(nzd)    [m2/s]    eddy diffusion coefficient
!     w(nzd)       [m/s]     vertical wind 
!
!  integration variable
!     dt           [s]  integration time step
!     c(nzd)       []   relative concentration of target gas
!    initial condition
!     c0(nzd)      []   initial concentration of target gas
!    upper or lower boundary condition
!     cb=f(t)      []   concentration of target gas at top or bottom of model layer
!  
!-----------------------------------------------------------------------
!
!  Expression
!    One-dimensional advection- (molecular and eddy) diffusion equation
!       for relative concentration (c=n/na)
!     dc/dt+(1/na)d(na*w*c)/dz-(1/na)d/dz[na*Dmol*{dc/dz+lambda*c}+na*Dedy*dc/dz]
!      here, 
!       lambda=(amsg-amsa)*gg/rg/temp+tdf*tmpd/temp
!
!     Equation to be solved,
!      dc/dt+delta*d2c/dz2+alpha*dc/dz+beta*c=0
!       here,
!         delta=-(Dmol+Dedy)
!         gamma=w-lambda*Dmol
!         alpha=gamma+(1/na)*d(na*delta)dz
!         beta=(1/na)*d(na*gamma)dz
!
!-----------------------------------------------------------------------
!  Boundary condition
!  Simulation type can be changed by setting isys to be 1 or 2
!
!    isys = 1 : L.B.=zero flux, U.B.=given conc. (for boundary atmosphere)
!    isys = 2 : L.B.=given conc. U.B.=zero flux  (for stratosphere)
!
!   Zero-flux: flux=0 at z=0 or top
!          w*c-{Dmol*(dc/dz+lambda*c)+Dedy*dc/dz}=0
!            i.e. delta*dc/dz+gamma*c=0
!
!   Given conc.: c=f(t) at z=top or 0
!           f(t) is given
!
!-----------------------------------------------------------------------
!============================================================================
!
      implicit none
!=====================================================================
!  	... Declaration
!======================================================================
!-----------------------------------------------------------------------
!  	... Parameters and constants
!-----------------------------------------------------------------------
      double precision, parameter :: cnp = 0.5d0      ! weighting Crank-Nicolson Method []
!      double precision, parameter :: cnp = 1.0d0      ! weighting Crank-Nicolson Method []
      double precision, parameter :: cnp1 = 1.d0-cnp  ! weighting Crank-Nicolson Method []
      double precision, parameter :: rg = 8.314510d0  ! gas constant [J/mol/K]
      double precision, parameter :: gg = 9.80665d0   ! gravitational acceleration [m/s2]
      double precision, parameter :: amsa = 28.966d-3 ! molar mass of air [kg/mol]
      double precision, parameter :: dco2 = 0.1247d-4 ! CO2 mol. diff. coef. [m2/s]
      double precision, parameter :: pres0 = 1.013d5  ! pressure at z=0 [Pa]
!-----------------------------------------------------------------------
!  	...  Main Variables
!-----------------------------------------------------------------------
      integer  ::  i, j, jj, nmxt, nmxd, nmxw
!
      double precision, allocatable, dimension(:) :: zdm, dzd ! (nzd)
      double precision, allocatable, dimension(:) :: zdb  ! (nzd+1)
      double precision, allocatable, dimension(:) :: dzm  ! (nzd-1)
!
      double precision, allocatable, dimension(:) :: temp, tmpd, na  ! (nzd)
      double precision, allocatable, dimension(:) :: dmol, dedy, w  ! (nzd)
!
      double precision, allocatable, dimension(:) :: c  ! (nzd)
!
!-----------------------------------------------------------------------
!  	...  Variables for reporting results
!-----------------------------------------------------------------------
      double precision, allocatable, dimension(:, :) :: cmx  ! (nzd,366) daily results stored in last year
      double precision, allocatable, dimension(:, :) :: csx  ! (nzd,366) daily results stored in first year
      double precision :: cb, cb1   ! upper boundary cb(j) & cb(j+1)
      double precision, allocatable, dimension(:, :) :: cyx  ! (nzd,iyp-iys+1) once per year (Jan. 1st in each year)
!
!     time-dependent atmospheric conditions (daily data)
      double precision, allocatable, dimension(:, :) :: tmx  ! (nzd,366) daily temp profiles
      double precision, allocatable, dimension(:, :) :: dmx  ! (nzd,366) daily dedy profiles
      double precision, allocatable, dimension(:, :) :: wmx  ! (nzd,366) daily w profiles
!
!-----------------------------------------------------------------------
!  	...  Working Variables
!-----------------------------------------------------------------------
      double precision, allocatable, dimension(:) :: alpha, beta, gamma, delta, lambda !(nzd)
      double precision, allocatable, dimension(:) :: b1, b2, b3, r, t ! (nzd) for tridiagonal system
      double precision, allocatable, dimension(:) :: q, s ! (nzd-1) for tridiagonal system
      double precision, allocatable, dimension(:) :: press  !(nzd)  pressure
      double precision, allocatable, dimension(:) :: xx, yy  !(nzd) working
      integer :: info  ! for SUBROUTINE dgtsv return
      character(20) a, ay, asy
      integer :: imp, kmp        ! monthly output
      double precision :: cmix(10)  ! perturbation
!
!-----------------------------------------------------------------------
!  	...  Initialization Variables
!-----------------------------------------------------------------------
      integer :: isys        ! sim. type (boundary conditions)
!
      integer :: nzd         ! max num of layers
      double precision :: zdt   ! height of top [m]
      double precision :: dzd1  ! thickness of layer [m]
      double precision :: dt    ! time step [s]
      double precision :: schk    ! stability check
!
      integer :: irinit ! flag of initial read from file (0=const. 1=read from file)
!
      integer :: irdedy ! flag of Dedy read from file (0=const. 1=read from file)
      double precision :: dedy0 ! if constant
!
      integer :: irtemp ! flag of temp read from file (0=const. 1=read from file)
      double precision :: temp0 ! if constant
!
      integer :: irwind ! flag of w read from file (0=const. 1=read from file)
      double precision :: w0 ! if constant
!
      double precision            :: amsg, sfd, tdf  
              ! molar mass, diff. scale factor, thermal diff. factor of gas [kg/mol],[],[]
!-----------------------------------------------------------------------
!  	...  Variables for boundary (time control)
!-----------------------------------------------------------------------
!  time format
!  idy : number of years during covered period
!  myy(n) : n-th year (ex. myy(1)=2001 A.D.)
!  mnm(n) : number of data in n-th year (ex. mnm(1)=250, 250 data in 2001 A.D.)
!  mdd(n,m) : number of days at m-th data in n-th year (ex. mdd(1,55)=121, 121 days at 55th data in 2001 A.D.)
!  atc(n,m) : conc. at m-th data in n-th year
!
      integer :: idy, iys, iyp, imn, idays ! total years, start year, stop year/month/day
      integer :: iatc ! flag for upper boundary scenario (0:constant, 1:given)
      double precision :: atcc  ! if constant scenario, given concentration
!
      integer :: myy(300), mnm(0:300), mdd(300, 366)
      double precision :: xac(366), xac1(366)   ! daily conc. for 1-year
      double precision :: atc(300, 366)   !
!
      integer :: idm, idm1, idd, id1, id2 ! time control working number
      integer :: iyy, ncc, ncc1, ndy, nds  ! time control working number
      double precision :: dtr, ytime, dyint   ! time control working number
!
!=====================================================================
!  	... Initialization
!======================================================================
!
!
!-----------------------------------------------------------------------
!	... Read parameter file
!-----------------------------------------------------------------------
      open (10, file='para.d')
      read (10, '()')       ! skip header
      read (10, *) isys  ! sym. type (1=boundary atm. 2=stratosphere)
      read (10, *) nzd   ! number of layer
      read (10, *) dzd1  ! layer thickness [m]
      read (10, *) dt    ! time step [s]
      read (10, *) iys, iyp, imn, idays  ! start year, stop year/month/day
      read (10, *) amsg, sfd, tdf  ! molar mass, diff. scale factor, thermal diff. factor of gas
      read (10, *) iatc, atcc  ! flag for upper boundary scenario & constant value
      read (10, *) irinit  ! flag of initial read from file (1=read from file "cinit.d")
      read (10, *) irdedy, dedy0  ! flag of Dedy read from file (1=read from file "dedy.d")
      read (10, *) irtemp, temp0  ! flag of temp read from file (1=read from file "temp.d")
      read (10, *) irwind, w0  ! flag of w read from file (1=read from file "wind.d")
      close (10)
!
      zdt = dble(nzd) * dzd1   ! height of top [m]
!   stability check, dz^2/(2*dt*Keddy) should be larger than 1, Keddy ~ 1.0 [m2/s]
      schk = dzd1 * dzd1 / ( 2.d0 * dt * 1.d0 )
!   print
      print *, 'sim. type', isys
      print *, 'number of layer', nzd
      print *, 'vertical resolution [m]', dzd1
      print *, 'top of atmosphere [m]', zdt
      print *, 'time step [s]', dt
      print *, 'stability check...should be larger than 1'
      print *, schk
!-----------------------------------------------------------------------
!	... allocate variables
!-----------------------------------------------------------------------
!
         allocate( zdm(nzd), dzd(nzd), zdb(nzd+1), dzm(nzd-1) )
         allocate( temp(nzd), tmpd(nzd), na(nzd), dmol(nzd), dedy(nzd), w(nzd), c(nzd) )
         allocate( csx(nzd,366) )
         allocate( cmx(nzd,366), tmx(nzd,366), dmx(nzd,366), wmx(nzd,366) )
         allocate( alpha(nzd), beta(nzd), gamma(nzd), delta(nzd), lambda(nzd) )
         allocate( b1(nzd), b2(nzd), b3(nzd), q(nzd-1), r(nzd), s(nzd-1), t(nzd) )
         allocate( press(nzd), xx(nzd), yy(nzd) )
!         allocate( cyx(nzd,iyp-iys+1) )
         allocate( cyx(nzd,1201) )
!
!-----------------------------------------------------------------------
!	... setting model structure
!-----------------------------------------------------------------------
!
      dzd(:) = dzd1  !  constant layer thickness
!
!  height of each boundary
      zdb(1) = 0.0d0
      do i = 1, nzd
        zdb(i+1) = zdb(i) + dzd(i)
      end do
!  
      zdm(:) = zdb(1:nzd) + 0.5d0 * dzd(:)  ! height at mid-point of layer
!  
      dzm(:) = zdm(2:nzd)-zdm(1:nzd-1)  ! thickness between mid-points of layer
!
!-----------------------------------------------------------------------
!	... upper boundary condition
!-----------------------------------------------------------------------
!  atc setting
      if (iatc == 1) then
        call setatmc (idy, myy, mdd, atc, mnm)
      else if (iatc == 0) then
        atc(:,:) = atcc
      end if
!
!-----------------------------------------------------------------------
!	... initial condition
!-----------------------------------------------------------------------
!    set xac(1) in start year as initial
      if (iatc == 1) then  ! if atc is time-dependent, 
        do i = 1, idy
          if ( myy(i) == iys) then
            idm = i      ! find simulation start year number
          end if
        end do
!
        call dayatc(idm, xac, myy, mdd, mnm, atc)  ! xac (quasi-daily data) in start year, iys
        ! use only xac(1) (first data in starting year) as initial profile
!
      else if (iatc == 0) then  ! if atc is constant
        xac(1) = atcc
      end if
!
!   set vertical conc. distribution, by xac(1) or by reading init file
      if (irinit == 0) then    
        c(:) = xac(1)     ! constant c as initial
      else
        call read_init (nzd, c)  ! read from file
      end if
!
!   atc output (for debug)
      if (iatc == 1) then
        open (30, file='atc_date_out.d')
        do j = 1, idy
          write (30, *) myy(j), mnm(j), (mdd(j,jj), jj=1, mnm(j))
        end do
        close (30)
      end if
!
!
!-----------------------------------------------------------------------
!	... time dependent atmospheric structure
!-----------------------------------------------------------------------
!
!   temperature, temp
      if (irtemp == 0) then
        tmx(:,:) = temp0  ! constant
      else
        call read_temp (nzd, zdm, tmx, nmxt)  ! read from file
      end if
!
!   eddy diffusion coefficient, dedy
      if (irdedy == 0) then
        dmx(:,:) = dedy0  ! constant
      else
        call read_dedy (nzd, zdm, dmx, nmxd)  ! read from file
      end if
!
!   advection flow, w
      if (irwind == 0) then
        wmx(:,:) = w0  ! constant
      else
        call read_w (nzd, zdm, wmx, nmxw)  ! read from file
      end if
!
!
!-----------------------------------------------------------------------
!	... output initialized variable
!-----------------------------------------------------------------------
!
      temp(:) = tmx(:,1)
      dedy(:) = dmx(:,1)
      w(:) = wmx(:,1)
!
!   total air quantity, na
      xx(:) = 1/temp(:)  ! 1/T
      call z_integ(nzd, dzd, xx, yy)  ! integration 1/T
      press(:) = pres0 * dexp( - amsa * gg /rg * yy(:) )
      na(:) = press(:) / rg / temp(:)
!
!   molecular diffusion coefficient, dmol
      dmol(:) = sfd * dco2 &
            * ( (temp(:) / 253.d0 ) ** 1.85d0 )  & ! temp dependent
            * ( pres0 / press(:) )                  ! press dependent
!
!   init. output
      open (20, file='init_out.d')
      write (20, *) 'i zdm na temp dmol dedy w press cinit'
      do i = 1, nzd
        write (20,'(I4, 7e12.4, e18.10)') i, zdm(i), na(i), temp(i) &
                        , dmol(i), dedy(i), w(i), press(i), c(i)
      end do
      close (20)
!
!   debug
!      print *, 'Initialize OK'
!      pause
!
!=====================================================================
!  	... Time integration
!======================================================================
!
      print *, 'Start Integration', iys
!
      dtr = 0.d0   ! the rest of time, will be added to next year
!
!-----------------------------------------------------------------------
      do iyy = iys, iyp - 1   ! integration loop for years
!-----------------------------------------------------------------------
!      time step control for each year
!-----------------------------------------------------------------------
       call leapyear (iyy, ndy )  ! total days number of this year
!
!     setting daily conc., xac(366) in corresponding year
        if (iatc==1) then
          idm1 = idm + iyy - iys
          call dayatc(idm1, xac, myy, mdd, mnm, atc)  ! xac in this year
          call dayatc(idm1+1, xac1, myy, mdd, mnm, atc)  ! xac in next year
        end if
!
!     time-step control
        ytime = dble(ndy) * 24.d0 * 60.d0 * 60.d0 + dtr - dt  ! total time of this year
        ncc = idint( ytime / dt ) + 1   ! total integration number of this year
        dyint = dt / 24.d0 / 60.d0 / 60.d0  ! time step (in unit of day)
        dtr = dmod( ytime, dt )   ! the rest of time, will be added to next year
!       debug
!        print *,  "ncc,dyint,dtr"
!        print *,  ncc,dyint,dtr
!        PAUSE
!
! ++++++ perturbation section start ++++++
!               if ( iyy == iys ) then
!                temp(20:30) = tmx(20:30,1) - 50.d0
!               end if
!               if ( iyy == iys ) then
!                c(25:30) = atcc
!               c(25:30) = ( 1.d-4 + 1.d0) * 0.00364d0 / 0.99636d0 * 0.781d0 * 2.d0  ! 100 permeg
!                cmix(1:10) = c(21:30)
!                c(21:30) = c(31:40)
!                c(31:40) = cmix(1:10)
!               end if
! ++++++ perturbation section end ++++++
!     annual output on Jan. 1st
               if ( iyy == iys ) then
                kmp=1
                cyx(:,kmp) = c(:)
               end if
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
        do idd = 1, ncc   ! integration loop for each time step in one year
!-----------------------------------------------------------------------
!
!    setting of upper boundary
          id1 = idint(dble(idd-1)*dyint) + 1  ! present days number 1 to 365 or 366
          id2 = idint(dble(idd)*dyint) + 1  ! next step days number 1 to 365 or 366
!
!
          if (iatc == 1) then   ! time-dependent
            cb = xac(id1)     ! present data
            if (idd .lt. ncc) then  ! next step data for boundary condition
              cb1 = xac(id2)  ! next step data
            else
              cb1 = xac1(1)   ! first data of next year
            end if
          else if (iatc == 0) then
            cb = atcc      ! constant
            cb1 = atcc      ! constant
          end if
!
      imp = mod(id1, 30)
      if (imp == 0 .and. id1 /= id2 ) then
       kmp = kmp + 1
       cyx(:,kmp) = c(:)
      end if
!-----------------------------------------------------------------------
!  calculate,
!    lambda=(amsg-amsa)*gg/rg/temp+tdf*tmpd/temp
!    delta=-(Dmol+Dedy)
!    gamma=w-lambda*Dmol
!    alpha=gamma+(1/na)*d(na*delta)dz
!    beta=(1/na)*d(na*gamma)dz
!
!   lambda
      call z_differen (nzd, dzd, temp, tmpd)
      lambda(:) = (amsg - amsa) * gg / rg / temp(:) + tdf * tmpd(:) / temp(:)
!   delta
      delta(:) = - (dmol(:) + dedy(:))
!   gamma
      gamma(:) = w(:) - lambda(:) * dmol(:)
!   alpha
      xx(:) = na(:) * delta(:)
      call z_differen (nzd, dzd, xx, yy)
      alpha(:) = gamma(:) + yy(:) / na(:)
!   beta
!      xx(:) = na(:) * gamma(:)
!  !!! d(na*W)/dz=0 , air mass can not be conserved for 1-D model when w is given
      xx(:) = - na(:) * lambda(:) * dmol(:)
      call z_differen (nzd, dzd, xx, yy)
      beta(:) = yy(:) / na(:)
!
!-----------------------------------------------------------------------
!  Tridiagonal System  
!        ( here dzd = dzm)
!    Super-diagonal (DU) element = s(1:nzd-1)
!          Diagonal (D)  element = r(1:nzd)
!      Sub-diagonal (DL) element = q(1:nzd-1)
!
!   (for i=2, nzd-1)
!     q(i-1)*c(j+1,i-1) + r(i)*c(j+1,i) + s(i)*c(j+1,i+1) = t(i)
!       q(i-1) = cnp*b1(i)
!       r(i) = 1/dt + cnp*b2(i)
!       s(i) = cnp*b3(i)
!       t(i) = -cnp1*b1(i)*c(j,i-1) - {-1/dt+cnp1*b2(i)}*c(j,i) - cnp1*b3(i)*c(j,i+1)
!         b1(i) = delta(i)/dzd(i)/dzd(i) - alpha(i)/dzd(i)/2
!         b2(i) = beta(i) - 2*delta(i)/dzd(i)/dzd(i)
!         b3(i) = delta(i)/dzd(i)/dzd(i) + alpha(i)/dzd(i)/2
!
        b1(2:nzd-1) = delta(2:nzd-1) / dzd(2:nzd-1) / dzd(2:nzd-1) &
                    - alpha(2:nzd-1) / dzd(2:nzd-1) / 2.d0
!
        b2(2:nzd-1) = beta(2:nzd-1)  &
                    - 2.d0 * delta(2:nzd-1) / dzd(2:nzd-1) / dzd(2:nzd-1)
!
        b3(2:nzd-1) = delta(2:nzd-1) / dzd(2:nzd-1) / dzd(2:nzd-1) &
                    + alpha(2:nzd-1) / dzd(2:nzd-1) / 2.d0
!
        q(1:nzd-2) = cnp * b1(2:nzd-1)
        r(2:nzd-1) = 1.d0 / dt + cnp * b2(2:nzd-1)
        s(2:nzd-1) = cnp * b3(2:nzd-1)
        t(2:nzd-1) = - cnp1 * b1(2:nzd-1) * c(1:nzd-2) &
                     - ( -1.d0 / dt + cnp1 * b2(2:nzd-1)) * c(2:nzd-1) &
                     - cnp1 * b3(2:nzd-1) * c(3:nzd)
!
!  Boundary elements
      if (isys == 1) then  !   isys = 1 (L.B.=zero-flux, U.B.=given conc.)
!   (for i=1)
!      r(1)*c(j+1,1) + s(1)*c(j+1,2) = t(1)
!        r(1) = -cnp*b2(1)
!        s(1) = cnp*b3(1)
!        t(1) = cnp1*{b2(1)*c(j,1) - b3(1)*c(j,2)}
!          b1(1) not used
!          b2(1) = delta(1)/dzd(1) - gamma(1)
!          b3(1) = delta(1)/dzd(1)
!
        b1(1) = 0.d0
        b2(1) = delta(1) / dzd(1) - gamma(1)
        b3(1) = delta(1) / dzd(1)
!
        r(1) = -cnp * b2(1)
        s(1) = cnp * b3(1)
        t(1) = cnp1 * ( b2(1) * c(1) - b3(1) * c(2) )
!
!   (for i=nzd)
!      q(nzd-1)*c(j+1,nzd-1) + r(nzd)*c(j+1,nzd) = t(nzd)
!        q(nzd-1) = 0
!        r(nzd) = cnp
!        t(nzd) = cnp*cb(j+1) + cnp1*{cb(j)-c(j,nzd)}
!
        q(nzd-1) = 0.d0
        r(nzd) = cnp
        t(nzd) = cnp * cb1 + cnp1 * ( cb - c(nzd) )
!
      else if (isys == 2) then  !   isys = 2 (U.B.=zero-flux, L.B.=given conc.)
!   (for i=1)
!      r(1)*c(j+1,1) + s(1)*c(j+1,2) = t(1)
!        r(1) = cnp
!        s(1) = 0
!        t(1) = cnp*cb(j+1) + cnp1*{cb(j)-c(j,1)}
!
        r(1) = cnp
        s(1) = 0.d0
        t(1) = cnp * cb1 + cnp1 * ( cb - c(1) )
!
!   (for i=nzd)
!      q(nzd-1)*c(j+1,nzd-1) + r(nzd)*c(j+1,nzd) = t(nzd)
!        q(nzd-1) = -cnp*b1(nzd)
!        r(nzd) = cnp*b2(nzd)
!        t(nzd) = cnp1*{b1(nzd)*c(j,nzd-1) - b2(nzd)*c(j,nzd)}
!          b1(nzd) = delta(nzd)/dzd(nzd)
!          b2(nzd) = delta(nzd)/dzd(nzd) + gamma(nzd)
!          b3(1) not used
!
        b1(nzd) = delta(nzd) / dzd(nzd)
        b2(nzd) = delta(nzd) / dzd(nzd) + gamma(nzd)
        b3(nzd) = 0.d0
!
        q(nzd-1) = -cnp * b1(nzd)
        r(nzd) = cnp * b2(nzd)
        t(nzd) = cnp1 * ( b1(nzd) * c(nzd-1) - b2(nzd) * c(nzd) )
!
      end if
!
!-----------------------------------------------------------------------
!    solve (integrate) 
      call dgtsv(nzd,1,q,r,s,t,nzd,info)
      c(:) = t(:)   ! dgtsv solution t(:)
!
!-----------------------------------------------------------------------
!       process after one-day integration
          if (id1 .ne. id2) then
!         print day number and store data
!            print *, iyy, id1
            cmx( : , id1 ) = c(:)
!
!             only for first year
               if ( iyy == iys ) then
                csx( : , id1 ) = c(:)
                nds = id1  ! total days of first year
               end if
!
!          change atmospheric condition
           if ( (ndy == 365 .and. id2 == 366) .or. &
                (ndy == 366 .and. id2 == 367) ) then
            temp(:) = tmx(:,1)
            dedy(:) = dmx(:,1)
            w(:) = wmx(:,1)
           else
            temp(:) = tmx(:,id2)
            dedy(:) = dmx(:,id2)
            w(:) = wmx(:,id2)
           end if
!
!          total air quantity, na
            xx(:) = 1/temp(:)  ! 1/T
            call z_integ(nzd, dzd, xx, yy)  ! integration 1/T
            press(:) = pres0 * dexp( - amsa * gg /rg * yy(:) )
            na(:) = press(:) / rg / temp(:)
!
!          molecular diffusion coefficient, dmol
            dmol(:) = sfd * dco2 &
               * ( (temp(:) / 253.d0 ) ** 1.85d0 )  & ! temp dependent
               * ( pres0 / press(:) )                  ! press dependent
!
          end if
!
!-----------------------------------------------------------------------
        end do  ! end of integration loop for each time step in one year
!-----------------------------------------------------------------------
!
            print *, iyy
!-----------------------------------------------------------------------
      end do   ! end of integration loop for years
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      iyy = iyp   !  last year integration
!-----------------------------------------------------------------------
!      time step control for each year
!-----------------------------------------------------------------------
!
       call leapyear (iyy, ndy )  ! total days number of this year
!
!     setting daily conc., xac(366) in corresponding year
      if (iatc == 1) then
        idm1 = idm + iyy - iys
        call dayatc(idm1, xac, myy, mdd, mnm, atc)  ! xac in this year
      end if
!
!     time-step control
        ytime = dble(ndy) * 24.d0 * 60.d0 * 60.d0 + dtr - dt  ! total time of this year
        ncc = idint( ytime / dt ) + 1   ! total integration number of this year
        dyint = dt / 24.d0 / 60.d0 / 60.d0  ! time step (in unit of day)
        dtr = dmod( ytime, dt )   ! the rest of time, will be added to next year

        ncc1 = idint( dble(ncc) * (dble(imn-1)  &
              * 30.5d0 + dble(idays)) / 366.d0 )  ! approx. integ. number of final year
!
!     ncc1 is just approximation of day number at last year/month/day
!     if stop time is ex. 2000/12/31, runtime error occurs 
!     because ncc1 become larger than ncc
        if (ncc1 .gt. ncc) then
          ncc1 = ncc
        end if
!
!     annual output on Jan. 1st
!      cyx(:,iyy-iys+1) = c(:)
!-----------------------------------------------------------------------
        do idd = 1, ncc1   ! integration loop for each time step in last year
!-----------------------------------------------------------------------
!
!    setting of upper boundary
          id1 = idint(dble(idd-1)*dyint) + 1  ! present days number 1 to 366
          id2 = idint(dble(idd)*dyint) + 1  ! next step days number 1 to 366
!
!
          if (iatc == 1) then   ! time-dependent
            cb = xac(id1)     ! present data
            if (idd .lt. ncc) then  ! next step data for boundary condition
              cb1 = xac(id2)  ! next step data
            else
              cb1 = xac(id1)   ! same data with present step
            end if
          else if (iatc == 0) then
            cb = atcc      ! constant
            cb1 = atcc      ! constant
          end if
!
      imp = mod(id1, 30)
      if (imp == 0 .and. id1 /= id2 ) then
       kmp = kmp + 1
       cyx(:,kmp) = c(:)
      end if
!-----------------------------------------------------------------------
!  calculate,
!    lambda=(amsg-amsa)*gg/rg/temp+tdf*tmpd/temp
!    delta=-(Dmol+Dedy)
!    gamma=w-lambda*Dmol
!    alpha=gamma+(1/na)*d(na*delta)dz
!    beta=(1/na)*d(na*gamma)dz
!
!   lambda
      call z_differen (nzd, dzd, temp, tmpd)
      lambda(:) = (amsg - amsa) * gg / rg / temp(:) + tdf * tmpd(:) / temp(:)
!   delta
      delta(:) = - (dmol(:) + dedy(:))
!   gamma
      gamma(:) = w(:) - lambda(:) * dmol(:)
!   alpha
      xx(:) = na(:) * delta(:)
      call z_differen (nzd, dzd, xx, yy)
      alpha(:) = gamma(:) + yy(:) / na(:)
!   beta
!      xx(:) = na(:) * gamma(:)
!  !!! d(na*W)/dz=0 , air mass can not be conserved for 1-D model when w is given
      xx(:) = - na(:) * lambda(:) * dmol(:)
      call z_differen (nzd, dzd, xx, yy)
      beta(:) = yy(:) / na(:)
!
!-----------------------------------------------------------------------
!  Tridiagonal System  
!        ( here dzd = dzm)
!    Super-diagonal (DU) element = s(1:nzd-1)
!          Diagonal (D)  element = r(1:nzd)
!      Sub-diagonal (DL) element = q(1:nzd-1)
!
!   (for i=2, nzd-1)
!     q(i-1)*c(j+1,i-1) + r(i)*c(j+1,i) + s(i)*c(j+1,i+1) = t(i)
!       q(i-1) = cnp*b1(i)
!       r(i) = 1/dt + cnp*b2(i)
!       s(i) = cnp*b3(i)
!       t(i) = -cnp1*b1(i)*c(j,i-1) - {-1/dt+cnp1*b2(i)}*c(j,i) - cnp1*b3(i)*c(j,i+1)
!         b1(i) = delta(i)/dzd(i)/dzd(i) - alpha(i)/dzd(i)/2
!         b2(i) = beta(i) - 2*delta(i)/dzd(i)/dzd(i)
!         b3(i) = delta(i)/dzd(i)/dzd(i) + alpha(i)/dzd(i)/2
!
        b1(2:nzd-1) = delta(2:nzd-1) / dzd(2:nzd-1) / dzd(2:nzd-1) &
                    - alpha(2:nzd-1) / dzd(2:nzd-1) / 2.d0
!
        b2(2:nzd-1) = beta(2:nzd-1)  &
                    - 2.d0 * delta(2:nzd-1) / dzd(2:nzd-1) / dzd(2:nzd-1)
!
        b3(2:nzd-1) = delta(2:nzd-1) / dzd(2:nzd-1) / dzd(2:nzd-1)  &
                    + alpha(2:nzd-1) / dzd(2:nzd-1) / 2.d0
!
        q(1:nzd-2) = cnp * b1(2:nzd-1)
        r(2:nzd-1) = 1.d0 / dt + cnp * b2(2:nzd-1)
        s(2:nzd-1) = cnp * b3(2:nzd-1)
        t(2:nzd-1) = - cnp1 * b1(2:nzd-1) * c(1:nzd-2)  &
                     - ( -1.d0 / dt + cnp1 * b2(2:nzd-1)) * c(2:nzd-1)  &
                     - cnp1 * b3(2:nzd-1) * c(3:nzd)
!  Boundary elements
      if (isys == 1) then  !   isys = 1 (L.B.=zero-flux, U.B.=given conc.)
!   (for i=1)
!      r(1)*c(j+1,1) + s(1)*c(j+1,2) = t(1)
!        r(1) = -cnp*b2(1)
!        s(1) = cnp*b3(1)
!        t(1) = cnp1*{b2(1)*c(j,1) - b3(1)*c(j,2)}
!          b1(1) not used
!          b2(1) = delta(1)/dzd(1) - gamma(1)
!          b3(1) = delta(1)/dzd(1)
!
        b1(1) = 0.d0
        b2(1) = delta(1) / dzd(1) - gamma(1)
        b3(1) = delta(1) / dzd(1)
!
        r(1) = -cnp * b2(1)
        s(1) = cnp * b3(1)
        t(1) = cnp1 * ( b2(1) * c(1) - b3(1) * c(2) )
!
!   (for i=nzd)
!      q(nzd-1)*c(j+1,nzd-1) + r(nzd)*c(j+1,nzd) = t(nzd)
!        q(nzd-1) = 0
!        r(nzd) = cnp
!        t(nzd) = cnp*cb(j+1) + cnp1*{cb(j)-c(j,nzd)}
!
        q(nzd-1) = 0.d0
        r(nzd) = cnp
        t(nzd) = cnp * cb1 + cnp1 * ( cb - c(nzd) )
!
      else if (isys == 2) then  !   isys = 2 (U.B.=zero-flux, L.B.=given conc.)
!   (for i=1)
!      r(1)*c(j+1,1) + s(1)*c(j+1,2) = t(1)
!        r(1) = cnp
!        s(1) = 0
!        t(1) = cnp*cb(j+1) + cnp1*{cb(j)-c(j,1)}
!
        r(1) = cnp
        s(1) = 0.d0
        t(1) = cnp * cb1 + cnp1 * ( cb - c(1) )
!
!   (for i=nzd)
!      q(nzd-1)*c(j+1,nzd-1) + r(nzd)*c(j+1,nzd) = t(nzd)
!        q(nzd-1) = -cnp*b1(nzd)
!        r(nzd) = cnp*b2(nzd)
!        t(nzd) = cnp1*{b1(nzd)*c(j,nzd-1) - b2(nzd)*c(j,nzd)}
!          b1(nzd) = delta(nzd)/dzd(nzd)
!          b2(nzd) = delta(nzd)/dzd(nzd) + gamma(nzd)
!          b3(1) not used
!
        b1(nzd) = delta(nzd) / dzd(nzd)
        b2(nzd) = delta(nzd) / dzd(nzd) + gamma(nzd)
        b3(nzd) = 0.d0
!
        q(nzd-1) = -cnp * b1(nzd)
        r(nzd) = cnp * b2(nzd)
        t(nzd) = cnp1 * ( b1(nzd) * c(nzd-1) - b2(nzd) * c(nzd) )
!
      end if
!
!-----------------------------------------------------------------------
!    solve (integrate) 
      call dgtsv(nzd,1,q,r,s,t,nzd,info)
      c(:) = t(:)   ! dgtsv solution t(:)
!
!-----------------------------------------------------------------------
!       process after one-day integration
          if (id1 .ne. id2) then
!         print day number and store data
!            print *, iyy, id1
            cmx( : , id1 ) = c(:)
!
!          change atmospheric condition
           if ( (ndy == 365 .and. id2 == 366) .or. &
                (ndy == 366 .and. id2 == 367) ) then
            temp(:) = tmx(:,1)
            dedy(:) = dmx(:,1)
            w(:) = wmx(:,1)
           else
            temp(:) = tmx(:,id2)
            dedy(:) = dmx(:,id2)
            w(:) = wmx(:,id2)
           end if
!
!          total air quantity, na
            xx(:) = 1/temp(:)  ! 1/T
            call z_integ(nzd, dzd, xx, yy)  ! integration 1/T
            press(:) = pres0 * dexp( - amsa * gg /rg * yy(:) )
            na(:) = press(:) / rg / temp(:)
!
!          molecular diffusion coefficient, dmol
            dmol(:) = sfd * dco2 &
               * ( (temp(:) / 253.d0 ) ** 1.85d0 )  & ! temp dependent
               * ( pres0 / press(:) )                  ! press dependent
!
          end if
!
!-----------------------------------------------------------------------
      end do  ! end of integration loop for each time step in last year
!-----------------------------------------------------------------------
!
            print *, iyy
!
!-----------------------------------------------------------------------
!	... output result 
!-----------------------------------------------------------------------
      write (a,'(I3)') id1-1    ! number to text, for format strings
!      write (ay,'(I3)') iyp-iys+1    ! number to text, for format strings
      write (ay,'(I4)') kmp    ! number to text, for format strings
      write (asy,'(I4)') nds    ! number to text, for format strings
!
!   output(final vertical profile)
      open (40, file='result.d')
      write (40, *) nzd, ' zdm c'
      do i = 1, nzd
        ! write (40,'(e14.6, e18.10)') zdm(i), c(i)
        write (40, *) zdm(i), c(i)
      end do
      close (40)
!
!    output (long-term variations of c, once per year, Jan. 1st in each year)
!      open (90, file='cyx.d')
!      write (90, *) nzd, iyp-iys+1, ' zdm cyx'
!      do i = 1, nzd
!        write (90, '('//ay//'E25.16)') cyx(i, 1:iyp-iys+1)
!      end do
!      close (90)
      open (90, file='cyx.d')
      write (90, *) nzd, kmp, ' zdm cyx'
      do i = 1, nzd
        write (90, '('//ay//'E25.16)') cyx(i, 1:kmp)
      end do
      close (90)
!
!    output (daily variations of c, temp, dedy, w in last year)
      open (50, file='cmx.d')
      write (50, *) nzd, id1-1, ' zdm cmx'
      do i = 1, nzd
        write (50, '('//a//'E25.16)') cmx(i, 1:id1-1)
      end do
      close (50)
!
      open (60, file='tmx.d')
      write (60, *) nzd, id1-1, ' zdm tmx'
      do i = 1, nzd
        write (60, '('//a//'E25.16)') tmx(i, 1:id1-1)
      end do
      close (60)
!
      open (70, file='dmx.d')
      write (70, *) nzd, id1-1, ' zdm dmx'
      do i = 1, nzd
        write (70, '('//a//'E25.16)') dmx(i, 1:id1-1)
      end do
      close (70)
!
      open (80, file='wmx.d')
      write (80, *) nzd, id1-1, ' zdm wmx'
      do i = 1, nzd
        write (80, '('//a//'E25.16)') wmx(i, 1:id1-1)
      end do
      close (80)
!
!   first year c data
      open (95, file='csx.d')
      write (95, *) nzd, nds, ' zdm csx'
      do i = 1, nzd
        write (95, '('//asy//'E25.16)') csx(i, 1:nds)
      end do
      close (95)
!
      print *, 'Finish Integration', iyp, imn, idays
!
! 
!
!=====================================================================
!  	... Subroutine
!======================================================================
!============================================================================
    contains
!============================================================================
      subroutine z_integ (nzd, dzd, x, integx)
         integer :: i, nzd
         double precision, dimension(nzd) :: dzd, x, integx
!
         integx(1) = x(1) * dzd(1) * 0.5d0
!
         do i = 2, nzd
           integx(i) = integx(i-1) + x(i) * dzd(i)
         end do
!
      end subroutine z_integ
!-----------------------------------------------------------------------
      subroutine z_differen (nzd, dzd, x, diffx)
         integer :: i, nzd
         double precision, dimension(nzd) :: dzd, x, diffx
!
         diffx(1) = (x(2) - x(1)) / dzd(1)
         diffx(nzd) = (x(nzd) - x(nzd-1)) / dzd(nzd)
!
         do i = 2, nzd-1
           diffx(i) = (x(i+1) - x(i-1)) / dzd(i) /2.d0
         end do
!
      end subroutine z_differen
!-----------------------------------------------------------------------
      subroutine read_temp (nzd, zdm, tmx, nmxt)
         integer :: i, nzd, nd0, nmxt
         double precision, allocatable, dimension(:) :: z0, t0
         double precision, dimension(nzd) :: zdm, temp
         double precision, dimension(nzd, 366) :: tmx  ! daily temp profiles
!
!  allow to read temperature as free vertical interval data
!   by using interpolate subroutine
!
         open (100, file = 'temp.d')
         read(100,*) nmxt, nd0
!
         allocate( z0(nd0), t0(nd0) )
!
	     do i = 1, nmxt
		  read(100,*) z0(:), t0(:)
!         interpolate
          call interpl(nd0, z0, t0, nzd, zdm, temp)
          tmx(:,i) = temp(:)
	     end do
!
         close(100)
!
!       if not given full-data, extrapolate with the same data
         if(nmxt .lt. 366) then
          do i = nmxt+1, 366
           tmx(:,i) = temp(:)
          end do
         end if
!
      end subroutine read_temp
!-----------------------------------------------------------------------
      subroutine read_dedy (nzd, zdm, dmx, nmxd)
         integer :: i, nzd, nd0, nmxd
         double precision, allocatable, dimension(:) :: z0, t0
         double precision, dimension(nzd) :: zdm, dedy
         double precision, dimension(nzd, 366) :: dmx  ! daily temp profiles
!
!  allow to read eddy diffusivity as free vertical interval data
!   by using interpolate subroutine
!
         open (110, file = 'dedy.d')
         read(110,*) nmxd, nd0
!
         allocate( z0(nd0), t0(nd0) )
!
	     do i = 1, nmxd
		  read(110,*) z0(:), t0(:)
!         interpolate
          call interpl(nd0, z0, t0, nzd, zdm, dedy)
          dmx(:,i) = dedy(:)
	     end do
!
         close(110)
!
!       if not given full-data, extrapolate with the same data
         if(nmxd .lt. 366) then
          do i = nmxd+1, 366
           dmx(:,i) = dedy(:)
          end do
         end if
!
      end subroutine read_dedy
!-----------------------------------------------------------------------
      subroutine read_w (nzd, zdm, wmx, nmxw)
         integer :: i, nzd, nd0, nmxw
         double precision, allocatable, dimension(:) :: z0, t0
         double precision, dimension(nzd) :: zdm, w
         double precision, dimension(nzd, 366) :: wmx  ! daily temp profiles
!
!  allow to read vertical wind as free vertical interval data
!   by using interpolate subroutine
!
         open (120, file = 'w.d')
         read(120,*) nmxw, nd0
!
         allocate( z0(nd0), t0(nd0) )
!
	     do i = 1, nmxw
		  read(120,*) z0(:), t0(:)
!         interpolate
          call interpl(nd0, z0, t0, nzd, zdm, w)
          wmx(:,i) = w(:)
	     end do
!
         close(120)
!
!       if not given full-data, extrapolate with the same data
         if(nmxw .lt. 366) then
          do i = nmxw+1, 366
           wmx(:,i) = w(:)
          end do
         end if
!
      end subroutine read_w
!-----------------------------------------------------------------------
      subroutine read_init (nzd, c)
         integer :: i, nzd
         double precision, dimension(nzd) :: c
         real :: dum
!
         open (140, file = 'cinit.d')
         read(140,'()')      ! skip heading, same format with result.d
	     do i = 1, nzd
		   read(140,*) dum, c(i)
	     end do
         close(140)
!
      end subroutine read_init
!-----------------------------------------------------------------------
    subroutine setatmc (idy, myy, mdd, atc, mnm)
!
      double precision :: yday = 365.24236d0
      double precision :: conv
      integer :: myy(300), mnm(0:300), mdd(300, 366)
      double precision :: atc(300, 366)
      double precision :: att, acc, dd
      integer :: nut, nac, my, my1, j, k, i, md, idy
!
!
!  time format
!  idy : number of years
!  myy(n) : n-th year (ex. myy(1)=2001 A.D.)
!  mnm(n) : number of data in n-th year (ex. mnm(1)=250, 250 data in 2001 A.D.)
!  mdd(n,m) : number of days at m-th data in n-th year (ex. mdd(1,55)=121, 121 days at 55th data in 2001 A.D.)
!  atc(n,m) : conc. at m-th data in n-th year
!
! read atmospheric conc.
!
! *** About data format
! *** data must cover entire period of target.
! *** start and stop year is set in main program (IYS and IYP).
! *** integration start at first day of starting year,
! *** and stop at end of final year
! *** so you should add some data before and after target period
! *** because daily data is interpolated from this sparse data.
!
!
      open (130, file='atmc.d')
      read (130, *) nut
      if ( nut == 6 ) then
        conv = 1.d-6
      else if ( nut == 0 ) then
        conv = 1.d0
      else if ( nut == 9 ) then
        conv = 1.d-9
      else if ( nut == 12 ) then
        conv = 1.d-12
      end if
!
      read (130, *) nac
      my1 = 1
      j = 0
      k = 0
      do i = 1, nac
        read (130, *) att, acc
        my = idint(att)
        if (my == my1) then
          j = j + 1
        else
          mnm(k) = j
          k = k + 1
          j = 1
        end if
        dd = att - dint(att)
        md = idnint(dd * yday) + 1
        myy(k) = my
        atc(k, j) = acc * conv
        mdd(k, j) = md
        my1 = my
      end do
      mnm(k) = j
      idy = k
      close (130)
!
    end subroutine setatmc
!-----------------------------------------------------------------------
    subroutine dayatc (idm, xac, myy, mdd, mnm, atc)
!
!   this subroutine will be called at starting every year, and set daily data
!
      integer :: i, j, idm
      integer :: myy(300), mnm(0:300), mdd(300, 366)
      integer :: ndy, ndy1, ix1
      double precision :: xac(366)   ! daily data in corresponding year
      double precision :: atc(300, 366)
      double precision :: y1, dum
!
      integer :: iyy, iur1, iur2, iur3  ! working, for leap year
!
      iyy = myy(idm)
      iur1 = mod(iyy, 4)
      iur2 = mod(iyy, 100)
      iur3 = mod(iyy, 400)
!
      if (iur1 == 0 .and. iur2 /= 0) then
        ndy = 366
      else if (iur3==0) then
        ndy = 366
      else
        ndy = 365
      end if
!
      iyy = myy(idm-1)
      iur1 = mod(iyy, 4)
      iur2 = mod(iyy, 100)
      iur3 = mod(iyy, 400)
      if (iur1==0 .and. iur2/=0) then
        ndy1 = 366
      else if (iur3==0) then
        ndy1 = 366
      else
        ndy1 = 365
      end if
!
! backward interpolation at atc starting point
      ix1 = ndy1 - mdd(idm-1, mnm(idm-1)) + mdd(idm, 1)
      y1 = (atc(idm,1)-atc(idm-1,mnm(idm-1)))/dble(ix1)
      dum = atc(idm-1, mnm(idm-1))
      do i = mdd(idm-1, mnm(idm-1)) + 1, ndy1
        dum = dum + y1
      end do
      do i = 1, mdd(idm, 1)
        dum = dum + y1
        xac(i) = dum
      end do
!
! forward interpolation at atc ending point
      ix1 = ndy - mdd(idm, mnm(idm)) + mdd(idm+1, 1)
      y1 = (atc(idm+1,1)-atc(idm,mnm(idm)))/dble(ix1)
      dum = atc(idm, mnm(idm))
      do i = mdd(idm, mnm(idm)) + 1, ndy
        dum = dum + y1
        xac(i) = dum
      end do
!
! inside interpolation
      if (mnm(idm)>1) then
        do j = 1, mnm(idm) - 1
          ix1 = mdd(idm, j+1) - mdd(idm, j)
          y1 = (atc(idm,j+1)-atc(idm,j))/dble(ix1)
          dum = atc(idm, j)
          do i = mdd(idm, j) + 1, mdd(idm, j+1)
            dum = dum + y1
            xac(i) = dum
          end do
        end do
      end if
!
    end subroutine dayatc
!-----------------------------------------------------------------------
    subroutine leapyear (iyy, ndy)
!
      integer :: ndy
      integer :: iyy, iur1, iur2, iur3  ! working, for leap year
!
      iur1 = mod(iyy, 4)
      iur2 = mod(iyy, 100)
      iur3 = mod(iyy, 400)
!
      if (iur1 == 0 .and. iur2 /= 0) then
        ndy = 366
      else if (iur3 == 0) then
        ndy = 366
      else
        ndy = 365
      end if
!
    end subroutine leapyear
!-----------------------------------------------------------------------
!  =====================================================================
      SUBROUTINE dgtsv( N, NRHS, DL, D, DU, B, LDB, INFO )
!
!  -- LAPACK driver routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   FACT, TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGTSV ', -INFO )
         RETURN
      END IF
!
      IF( N.EQ.0 )   RETURN
!
      IF( NRHS.EQ.1 ) THEN
         DO 10 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
!
!              No row interchange required
!
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
   10    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      ELSE
         DO 40 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
!
!              No row interchange required
!
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 20 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   20             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               DO 30 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   30          CONTINUE
            END IF
   40    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 50 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   50             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               DO 60 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   60          CONTINUE
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      END IF
!
!     Back solve with the matrix U from the factorization.
!
      IF( NRHS.LE.2 ) THEN
         J = 1
   70    CONTINUE
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 ) &
            B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 80 I = N - 2, 1, -1
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*  &
                        B( I+2, J ) ) / D( I )
   80    CONTINUE
         IF( J.LT.NRHS ) THEN
            J = J + 1
            GO TO 70
         END IF
      ELSE
         DO 100 J = 1, NRHS
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 )  &
               B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /  &
                             D( N-1 )
            DO 90 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*  &
                           B( I+2, J ) ) / D( I )
   90       CONTINUE
  100    CONTINUE
      END IF
!
      RETURN
!
!     End of dgtsv
!
      END subroutine dgtsv
!-----------------------------------------------------------------------
!  =====================================================================
      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
!     ..
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ',  &
            'an illegal value' )
!
!     End of XERBLA
!
      END subroutine XERBLA
!-----------------------------------------------------------------------
    subroutine interpl(nd0, z0, fz0, nzd, zdm, fzdm)
      implicit none
      integer :: i, j, nd0, nzd
      double precision :: z0(nd0), fz0(nd0)
      double precision :: zdm(nzd), fzdm(nzd)
!
!  interpolate at model height, zdm
      do i = 1, nzd
        if (zdm(i)<z0(1)) then
          fzdm(i) = fz0(1)
        else if (zdm(i)>=z0(nd0)) then
          fzdm(i) = fz0(nd0)
        else
          do j = 1, nd0
            if (zdm(i)>=z0(j) .and. zdm(i)<z0(j+1)) then
              fzdm(i) = fz0(j) + (fz0(j+1)-fz0(j))  &
                       /(z0(j+1)-z0(j))*(zdm(i)-z0(j))
              Go To 100
            end if
          end do
        end if
100   end do

    end subroutine 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
    end program moldif

