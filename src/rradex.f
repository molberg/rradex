      subroutine rdxinput(mol, f, tk, bg, ncol, dv)
      implicit none
      integer length
      character*80 mol
      real*8 f(2), tk, bg, ncol, dv

      include 'radex.inc'

      call rdxdefaults
      molfile = mol(1:length(mol))
      fmin = f(1)
      fmax = f(2)
      tkin = tk
      tbg = bg
      cdmol = ncol
      deltav = dv*1.0e5
      return
      end

      subroutine rdxdensity(partner, d)
      implicit none
      include 'radex.inc'
      character*6 partner
      real*8 d
      integer id

      if ((d.lt.1.0e-3).or.(d.gt.1.0e13)) then
         call rwarn('densities need to be between 1e-3 and 1e13')
         return
      endif

      id = 0
      if ((partner(1:1).eq.'h').or.(partner(1:1).eq.'H')) id=5
      if ((partner(1:2).eq.'h2').or.(partner(1:2).eq.'H2')) id=1
      if ((partner(1:1).eq.'p').or.(partner(1:1).eq.'P')) id=2
      if ((partner(1:1).eq.'o').or.(partner(1:1).eq.'O')) id=3
      if ((partner(1:1).eq.'e').or.(partner(1:1).eq.'E')) id=4
      if ((partner(1:2).eq.'he').or.(partner(1:2).eq.'He')) id=6
      if ((partner(1:2).eq.'h+').or.(partner(1:2).eq.'H+')) id=7
      if (id.eq.0) then
         call rwarn('Unknown species.'
     $            //' Choose from: H2 p-H2 o-H2 e H He H+')
      else
         density(id) = d
         if ((density(2).gt.0.0).or.(density(3).gt.0.0))
     $        density(1) = density(2)+density(3)
      endif
      return
      end

      subroutine rdxoutput(up,low,nrg,ghz,tx,tau,ant,pu,pl,flx,erg)
      implicit none
      include 'radex.inc'
      
      integer up(maxline), low(maxline)
      real*8 nrg(maxline), ghz(maxline), tx(maxline), tau(maxline)
      real*8 ant(maxline), pu(maxline), pl(maxline)
      real*8 flx(maxline), erg(maxline)

      integer iline    ! to loop over lines
      integer m,n      ! upper & lower level of the line

      real*8 xt        ! frequency cubed
      real*8 hnu       ! photon energy
      real*8 bnutex    ! line source function
      real*8 ftau      ! exp(-tau)
      real*8 toti      ! background intensity
      real*8 tbl       ! black body temperature
      real*8 wh        ! Planck correction
      real*8 tback     ! background temperature
      real*8 ta        ! line antenna temperature
      real*8 tr        ! line radiation temperature
      real*8 beta,escprob ! escape probability
      external escprob
      real*8 bnu       ! Planck function
      real*8 kkms      ! line integrated intensity (K km/s)
      real*8 ergs      ! line flux (erg / s / cm^2)
      real*8 wavel     ! line wavelength (micron)

      do iline=1,nline
         m  = iupp(iline)
         n  = ilow(iline)
         xt = xnu(iline)**3.
c     Calculate source function
         hnu = fk*xnu(iline)/tex(iline)
         if(hnu.ge.160.0d0) then
            bnutex = 0.0d0
         else
            bnutex = thc*xt/(dexp(fk*xnu(iline)/tex(iline))-1.d0)
         endif
c     Calculate line brightness in excess of background
         ftau = 0.0d0
         if(abs(taul(iline)).le.3.d2) ftau = dexp(-taul(iline))
         toti = backi(iline)*ftau+bnutex*(1.d0-ftau)
         if(toti.eq.0.0d0) then
            tbl = 0.0d0
         else
            wh = thc*xt/toti+1.d0
            if(wh.le.0.d0) then
               tbl = toti/(thc*xnu(iline)*xnu(iline)/fk)
            else
               tbl = fk*xnu(iline)/dlog(wh)
            endif
         endif
         if(backi(iline).eq.0.0d0) then
            tback = 0.0d0
         else
            tback = fk*xnu(iline)/dlog(thc*xt/backi(iline)+1.d0)
         endif
c     Calculate antenna temperature
         tbl = tbl-tback
         hnu = fk*xnu(iline)
         if(abs(tback/hnu).le.0.02) then
            ta = toti
         else
            ta = toti-backi(iline)
         endif
         ta = ta/(thc*xnu(iline)*xnu(iline)/fk)
c     Calculate radiation temperature
         beta = escprob(taul(iline))
         bnu  = totalb(iline)*beta+(1.d0-beta)*bnutex
         if(bnu.eq.0.0d0) then
            tr = totalb(iline)
         else
            wh = thc*xt/bnu+1.0
            if(wh.le.0.0) then
               tr = bnu/(thc*xnu(iline)*xnu(iline)/fk)
            else
               tr = fk*xnu(iline)/dlog(wh)
            endif
         endif

c     Check if line within output freq range
         up(iline) = m
         low(iline) = n
         nrg(iline) = eup(iline)
         ghz(iline) = spfreq(iline)
         tx(iline) = tex(iline)
         tau(iline) = taul(iline)
         pu(iline) = xpop(m)
         pl(iline) = xpop(n)
         if (spfreq(iline).lt.fmax.and.spfreq(iline).gt.fmin) then
            wavel = clight / spfreq(iline) / 1.0d5 ! unit =  micron
            kkms  = 1.0645*deltav*ta
            ergs  = fgaus*kboltz*deltav*ta*(xnu(iline)**3.)
c     Line flux in K*cm/s and in erg/s/cm^2 
            ant(iline) = ta
            flx(iline) = kkms/1.0e5
            erg(iline) = ergs
c            if (dabs((tex(iline))).lt.1000.0) then
c            write(*,113) qnum(m),qnum(n),eup(iline),spfreq(iline),wavel,
c     $         tex(iline),taul(iline),ta,xpop(m),xpop(n),kkms/1.0d5,ergs
c            else
c            write(*,114) qnum(m),qnum(n),eup(iline),spfreq(iline),wavel,
c     $         tex(iline),taul(iline),ta,xpop(m),xpop(n),kkms/1.0d5,ergs
c            endif
c 113        format(a,' -- ',a,f8.1,2(2x,f10.4),1x,f8.3,6(1x,1pe10.3))
c 114        format(a,' -- ',a,f8.1,2(2x,f10.4),7(1x,1pe10.3))
         endif
      enddo

      return
      end

      subroutine rdxrun(niter)
      implicit none

      include 'radex.inc'

      integer niter
      logical conv

c     Read data file
      call backrad

      niter = 0
      conv  = .false.

c     Set up rate matrix, splitting it in radiative and collisional parts
c     Invert rate matrix to get `thin' starting condition
      call matrix(niter,conv)

c     Start iterating
      do 22 niter=1,maxiter
c     Invert rate matrix using escape probability for line photons
         call matrix(niter,conv)
         if (conv) goto 23
 22   continue

 23   return
      end
      
      subroutine rdxmoldata(level, line, coll, part)
      implicit none
      integer level, line, coll, part
      include 'radex.inc'

      level = nlev
      line = nline
      coll = ncoll
      part = npart
      return
      end

      function length(str)
c     Returns the length of a string
      integer length,maxl,i
c      parameter(maxl=120)
      character*(*) str

      do i=1,len(str)
         if (str(i:i).eq.' ') then
            length=i-1
            return
         endif
      enddo
      call rexit('string length problem')
      end

      subroutine rdxdefaults
      implicit none
      include 'radex.inc'

c     Set physical parameters to default values

      integer ipart  ! to loop over collision partners

      tkin   = 30.0
      tbg    = 2.73
      cdmol  = 1.0e13
      deltav = 1.0

      density(1) = 1.0e5
      do ipart=2,maxpart
         density(ipart) = 0.0
      enddo

      return
      end
