      implicit real*8 (a-h,m,o-z)

      real*8 fp(6)
      real*8 moon_xfe_1(6), moon_xfe_2(6)
      real*8 moon_dist_1(6), moon_dist_2(6)
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
      call random_number(moon_xfe_2)
      call random_number(moon_dist_1)
      call random_number(moon_dist_2)

C --------- scale and order values ----------------------------
      do k=1,6
         fp(k) = 0.2 * fp(k) + 0.55
         moon_xfe_1(k) = 1.0 * moon_xfe_1(k) + 0.0
         moon_xfe_2(k) = 1.0 * moon_xfe_2(k) + 0.0
         moon_dist_1(k) = 5.0 * moon_dist_1(k) + 2.9
         moon_dist_2(k) = 5.0 * moon_dist_2(k) + 2.9
         if (moon_dist_1(k) .gt. moon_dist_2(k)) then
            tmp = moon_dist_2(k)
            moon_dist_2(k) = moon_dist_1(k)
            moon_dist_1(k) = tmp
         endif
      end do

C --------- write results to file -----------------------------
      open(unit=1, file='PARAM_MSE.txt')
123   format(6(1f10.6,1x),1x,a12)
      write(1,123) (fp(k), k=1,6), 'p_factor'
      write(1,123) (moon_xfe_1(k), k=1,6),
     & 'moon_xfe_1'
      write(1,123) (moon_xfe_2(k), k=1,6),
     & 'moon_xfe_2'
      write(1,123) (moon_dist_1(k), k=1,6),
     & 'moon_dist_1'
      write(1,123) (moon_dist_2(k), k=1,6),
     & 'moon_dist_2'

      write(1,*) 'PARAM_MSE.txt generated from gen_param_MSE_moon.f'
      close(1)

      end
