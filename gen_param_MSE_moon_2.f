      implicit real*8 (a-h,m,o-z)

      real*8 fp(4)
      real*8 X_TM(4)
      real*8 X_TC(4)
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
      call random_number(X_TM)
      call random_number(X_TC)

C --------- scale and order values ----------------------------
      do k=1,4
         fp(k) = 1.0 * fp(k) + 0.0
         X_TM(k) = 0.5 * X_TM(k) + 0.0
         X_TC(k) = 0.5 * X_TC(k) + 0.0
      end do

C --------- write results to file -----------------------------
      open(unit=1, file='PARAM_MSE.txt')
123   format(4(1f10.6,1x),1x,a12)
      write(1,123) (fp(k), k=1,4), 'p_factor'
      write(1,123) (X_TM(k), k=1,4),
     & 'Theia mantle'
      write(1,123) (X_TC(k), k=1,4),
     & 'Theia core'

      write(1,*) 'PARAM_MSE.txt generated from gen_param_MSE_moon_2.f'
      close(1)

      end
