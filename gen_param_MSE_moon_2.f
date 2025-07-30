      implicit real*8 (a-h,m,o-z)

      real*8 fp(4)
      real*8 moon_xfe_1(4)
      real*8 Theia_fraction(4)
      real*8 tmp

      integer un, n, k, istat
      integer, allocatable :: seed(:)

C --------- set up random seed portably -----------------------
      call random_seed(size = n)
      allocate(seed(n))
      open(newunit=un, file="/dev/urandom", access="stream",
     & form="unformatted", action="read", status="old", iostat=istat)
      read(un) seed
      close(un)
      call random_seed(put=seed)

C --------- generate random numbers ---------------------------
      call random_number(fp)
      call random_number(moon_xfe_1)
      call random_number(Theia_fraction)

C --------- scale and order values ----------------------------
      do k=1,4
         fp(k) = 1.0 * fp(k) + 0.0
         moon_xfe_1(k) = 1.0 * moon_xfe_1(k) + 0.0
         Theia_fraction(k) = 1.0 * Theia_fraction(k) + 0.0
      end do

C --------- write results to file -----------------------------
      open(unit=1, file='PARAM_MSE.txt')
123   format(4(1f10.6,1x),1x,a12)
      write(1,123) (fp(k), k=1,4), 'p_factor'
      write(1,123) (moon_xfe_1(k), k=1,4),
     & 'moon_xfe_1'
      write(1,123) (Theia_fraction(k), k=1,4),
     & 'Theia_frac'

      write(1,*) 'PARAM_MSE.txt generated from gen_param_MSE_moon_2.f'
      close(1)

      end
