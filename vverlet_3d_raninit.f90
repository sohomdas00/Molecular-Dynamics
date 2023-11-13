!---------------------------3D LJ Mixture at RHO=0.3,T=0.5---------------------------------
!-----------------------------implementing NVE MD-----------------------------------
!-------------------------particles are intially randomly placed--------------------------
	program molecular_dynamics
	implicit none
	integer*4,parameter :: n = 2000	! no. of particles
	integer*4 :: npart,i,cntmax,idum
	integer*4 :: cnt
	real*8 :: tmax,rho,en,rcut,en2,en1,en_var,en_sdev,e_drift
	real*8 :: l,sigma,ecut,fcorr,dt,t,etot
	real*8 :: fx(n),fy(n),fz(n),x(n),y(n),z(n),vx(n),vy(n),vz(n),temp
	character(len=100) :: name1,name2
	
	common /block1/ en2,en1,en_var,en_sdev,e_drift
	common /block2/ l,sigma,ecut,fcorr,dt,t,etot
	common /block3/ cnt
    
	name1 = 'energies_time_5e-03'
	open(15,file=name1,status='unknown')
        open(16,file='dt_energy-variance',status='unknown')
	idum = 45678

	rho = 0.3d0
	temp = 0.5d0
	dt = 0.005d0
	
	l = int((dble(n)/rho)**(1d0/3d0))	! length of cubical box
	sigma = 1.0d0	! diameter of particles
	
	rcut = 2.5d0*sigma
	t=0.0d0
	ecut = 4.0d0*( 1.0d0/rcut**12d0 - 1.0d0/rcut**6d0 )
	fcorr = 48.0d0*( 0.5d0/rcut**7d0 - 1.0d0/rcut**13d0 )
	
	cntmax = 3000

	call init(n,rho,x,y,z,vx,vy,vz,idum,temp)
	
	tmax = 50.0d0
	!cntmax = int(tmax/dt)		! energy variance is calculated in all the runs for equal time
	!print *, cntmax
	en2 = 0.d0
	en1 = 0.d0
	en_var = 0.d0
	en_sdev = 0.d0

	call force(n,fx,fy,fz,x,y,z,en)	! calculates forces on the initial state point
	cnt = 1					! number of time steps progressed
	do while (cnt.le.cntmax)
	!  if (cnt.eq.cntmax) then
	!    name1='snap_'
	!    call addnumtostring(name1,cnt)
	!    open (51,file=name1,status='unknown')
	!    do i=1,n
	!   	write(51,'(6(f13.10,2X))') x(i),y(i),z(i),vx(i),vy(i),vz(i)
	!    enddo
	!  endif
	call integrate(n,fx,fy,fz,x,y,z,vx,vy,vz,en)
	call sample(etot,cnt,cntmax,dt)	    
	cnt = cnt + 1
	t = t+dt  	
	print*, cnt
	enddo
	!dt = dt + 1d-03
	!print *, dt
	!if (dt .le. 0.010001) goto 15
	end program molecular_dynamics
!----------------------------------------------------------------------------------------------------- 
!*****************************************************************************************************
      subroutine force(n,fx,fy,fz,x,y,z,en)
      implicit none
      integer*4 :: n
      integer*4 :: i,j,cnt
      real*8 :: fx(n),fy(n),fz(n),x(n),y(n),z(n),en,sigma
      real*8 :: rcut,xr,yr,zr,xl,yl,zl,r2,ecut,fcorr,r2i,r6i,ff,r1,l,dt,t,etot
      common /block2/ l,sigma,ecut,fcorr,dt,t,etot
      common /block3/ cnt

        rcut = 2.5d0*sigma
	xl=l
	yl=l
	zl=l
      en=0.d0	
      do i=1,n
        fx(i)=0.d0
        fy(i)=0.d0
        fz(i)=0.d0
      enddo

     do i=1,n-1
       do j=i+1,n
        xr = x(i) - x(j)
        yr = y(i) - y(j)
	zr = z(i) - z(j)
!_______minimum image convention__________
        xr = xr - xl*dnint(xr/xl)
        yr = yr - yl*dnint(yr/yl)
	zr = zr - zl*dnint(zr/zl)
!_________________________________________
        r2=xr**2d0 + yr**2d0 + zr**2d0    
        if (r2.le.rcut**2d0) then
	  r1=dsqrt(r2)
	  !print *, r1
          r2i=1.0d0/r2
	  r6i=r2i**3d0
          ff=48.0d0*r2i*r6i*(r6i-0.5d0)
          fx(i) = fx(i) + (ff*xr + fcorr*xr/r1)
          fx(j) = fx(j) - (ff*xr + fcorr*xr/r1)
          fy(i) = fy(i) + (ff*yr + fcorr*yr/r1)
          fy(j) = fy(j) - (ff*yr + fcorr*yr/r1)
	  fz(i) = fz(i) + (ff*zr + fcorr*zr/r1)
	  fz(j) = fz(j) - (ff*zr + fcorr*zr/r1) 
	  en = en + (4.0d0*r6i*(r6i-1.0d0)) - ecut - fcorr*(r1-rcut)
        endif
       enddo
      enddo
      return
      end     
!*****************************************************************************************************
!-----------------------------------------------------------------------------------------------------
	subroutine integrate(n,fx,fy,fz,x,y,z,vx,vy,vz,en)
	implicit none
	integer*4 :: n 
	integer*4 :: i,cnt
	real*8 :: sumvx,sumvx2,sumvy,sumvy2,sumvz,sumvz2,dt,en,vxn,vyn,vzn,temp1,etot,ketot,entot,t,sumv
	real*8 :: x(n),y(n),z(n),fx(n),fy(n),fz(n),vx(n),vy(n),vz(n),sigma,l
        real*8 :: fxi(n),fyi(n),fzi(n),ecut,fcorr,dx,dy,dz,sumv2
	character(len=100) :: name1
	common /block2/ l,sigma,ecut,fcorr,dt,t,etot
	common /block3/ cnt

	sumvx = 0.d0
	sumvy = 0.d0
	sumvz = 0.d0
	sumvx2 = 0.d0
	sumvy2 = 0.d0
	sumvz2 = 0.d0

	do i=1,n
	  fxi(i)=fx(i)
	  fyi(i)=fy(i)	! stores forces at the ith time step
	  fzi(i)=fz(i)
	enddo
	do i=1,n	! position half update
	  x(i) = x(i) + vx(i)*dt + 0.5d0*fx(i)*dt**2d0
	  y(i) = y(i) + vy(i)*dt + 0.5d0*fy(i)*dt**2d0
	  z(i) = z(i) + vz(i)*dt + 0.5d0*fz(i)*dt**2d0
	end do

	do i=1,n
!_______periodic boundary conditions___________________
	  if (x(i) .gt. l) x(i) = x(i) - l
	  if (x(i) .lt. 0.d0)x(i) = x(i) + l
	  if (y(i) .gt. l) y(i) = y(i) - l
	  if (y(i) .lt. 0.d0) y(i) = y(i) + l
	  if (z(i) .gt. l) z(i) = z(i) - l
	  if (z(i) .lt. 0.d0) z(i) = z(i) + l
!______________________________________________________
	enddo	
	call force(n,fx,fy,fz,x,y,z,en)	! stores forces at (i+1)th time step
	do i=1,n
	  vx(i) = vx(i) + 0.5d0*(fx(i)+fxi(i))*dt
	  vy(i) = vy(i) + 0.5d0*(fy(i)+fyi(i))*dt	
	  vz(i) = vz(i) + 0.5d0*(fz(i)+fzi(i))*dt
	  sumvx = sumvx + vx(i)
          sumvy = sumvy + vy(i)
	  sumvz = sumvz + vz(i)
          sumvx2 = sumvx2 + vx(i)**2d0
          sumvy2 = sumvy2 + vy(i)**2d0
	  sumvz2 = sumvz2 + vz(i)**2d0
	enddo
        sumv2 = sumvx2 + sumvy2 + sumvz2
	temp1 = sumv2/(3.d0*dble(n))	! instantaneous temperature
	ketot = 0.5d0*sumv2 		    ! total ke 
	etot = (en + ketot)/dble(n)	    ! energy per particle
	entot = en/dble(n)	            ! pe per particle
        ketot = ketot/dble(n)           ! ke per particle
	write(15,'(i5,X,5(f13.10,2X))') cnt,cnt*dt,ketot,entot,etot,temp1
	return
	end
!---------------------------------------------------------------------------------------------
!*********************************************************************************************
	subroutine sample(etot,cnt,cntmax,dt)
	implicit none
	real*8 :: etot,en2,en1,en_var,en_sdev,dt,e_drift,e0
	integer*4 :: cnt,cntmax
	common /block1/ en2,en1,en_var,en_sdev,e_drift
	!common /block2/ l,sigma,ecut,fcorr,dt,t,etot
	!common /block3/ cnt
	
	if (cnt.eq.1)e0 = etot
	e_drift = e_drift + dabs((e0-etot)/e0)
	en2 = en2 + etot**2d0
	en1 = en1 + etot
	!print *, en2,en1
	if (cnt.eq.cntmax) then
	  en2 = en2/dble(cntmax)
	  en1 = en1/dble(cntmax)
	  e_drift = e_drift/(cntmax)
	  print *, en2,en1
	  en_var = en2 - en1*en1
	  en_sdev = dsqrt(en_var)
	  write(16,'(4(f13.10X))') dt,en_var,en_sdev,e_drift
	endif
	return
	end
!---------------------------------------------------------------------------------------------
!*********************************************************************************************
	  subroutine init(n,rho,x,y,z,vx,vy,vz,idum,temp)
	  implicit none
	  integer*4 :: n	! number of particles
         real*8 :: sigma, dist, npart, rho, x1, y1, x2, y2, rsq, ct, x(n), y(n),z(n), rij_sq, ran3
         real*8 :: vx(n),vy(n),vz(n),fs,sumv2,sumvx,sumvy,sumvz,sumvx2,sumvy2,sumvz2,temp, fsx,fsy,fsz
         integer*4, dimension(:,:,:), allocatable :: cell
         integer*4 :: idum, i, j, k, i1, i2, p, q, r, l, p1, q1
         real*8 :: x00,y00,z00,x0,y0,z0,dx,dy,dz
	  npart = dble(n)      !number of particles

 
	  !temp = 2.0  	  !reduced temperature

        l = int((npart/rho)**(1d0/3d0))   ! dimension of square box
        allocate (cell(l,l,l))
	    !allocate (part(l,l))  
        	  
        sigma = 1.d0	! radius of particle
	    i1=1	! dummy index stores the number of particles in box
        dist = sigma**2d0	! square of the minimum permissible distance b/w two particles
	  
	  open(13,file='rho_0.3',status='unknown')
!------ randomly generating 1st particle -----------------------------------------------------
	x00 = l*ran3(idum)
	y00 = l*ran3(idum)
	z00 = l*ran3(idum)
	x(i1) = x00
	y(i1) = y00
	z(i1) = z00
 !------generating subsequent particles -----------------------------------------------------
15	x0 = l*ran3(idum)
	y0 = l*ran3(idum)
	z0 = l*ran3(idum)
	do i=1,i1
	  dx = x0-x(i)
	  dy = y0-y(i)
	  dz = z0-z(i)
	  dx = dx-l*nint(dx/l)
	  dy = dy-l*nint(dy/l)
	  dz = dz-l*nint(dz/l)
	  rsq = dx**2 + dy**2 + dz**2
	  rsq = sqrt(rsq)
	  if (rsq.lt.sigma) goto 15
	enddo
 	!p = int(l*ran3(idum)) + 1
	!q = int(l*ran3(idum)) + 1
	!r = int(l*ran3(idum)) + 1
	!if (p.eq.l+1)p=1
	!if (q.eq.l+1)q=1
	!if (r.eq.l+1)r=1
	!if (cell(p,q,r).eq.1) goto 15
	i1 = i1 + 1
	x(i1) = x0
	y(i1) = y0
	z(i1) = z0
	!cell(p,q,r) = 1
	!print *,p,q,r,cell(p,q,r)
	if (i1.lt.n) goto 15

!-------write particle co-ordinates in file ------------------------------------------------
	do i=1,n
	  write(13,*) x(i),y(i),z(i)
	enddo	
!-------------------------------------------------------------------------------------------
        !checking pair distances and printing them
	!do i=1,n-1
       ! do j=i+1,n
	  !  rij_sq = (x(i)-x(j))**2d0 + (y(i)-y(j))**2d0 + (z(i)-z(j))**2d0
	   ! write (*,*) i,j,dsqrt(rij_sq)
	  !enddo
	!enddo
        sumvx=0.d0
        sumvy=0.d0
	sumvz=0.d0
        sumvx2=0.d0
        sumvy2=0.d0
	  sumvz2=0.d0
        do i=1,n
          vx(i)=ran3(idum)-0.5d0
          vy(i)=ran3(idum)-0.5d0
	  vz(i)=ran3(idum)-0.5d0
	    !print*, vx(i),vy(i),vz(i)
          sumvx=sumvx+vx(i)
          sumvy=sumvy+vy(i)
	  sumvz=sumvz+vz(i)
          sumvx2=sumvx2+vx(i)**2d0
          sumvy2=sumvy2+vy(i)**2d0
	  sumvz2=sumvz2+vz(i)**2d0
        enddo
        sumvx=sumvx/dble(n)
        sumvy=sumvy/dble(n)
	sumvz=sumvz/dble(n)
        sumvx2=sumvx2/dble(n)
        sumvy2=sumvy2/dble(n)
	sumvz2=sumvz2/dble(n)
        sumv2=sumvx2+sumvy2+sumvz2
        fsx=dsqrt(temp/sumvx2)
        fsy=dsqrt(temp/sumvy2)
        fsz=dsqrt(temp/sumvz2)
	do i=1,n
	  vx(i)=(vx(i)-sumvx)*fsx
	  vy(i)=(vy(i)-sumvy)*fsy
	  vz(i)=(vz(i)-sumvz)*fsz
	  !print*, vx(i),vy(i),vz(i)
	enddo
 	return
	end
!*********************************************************************************************
!---------------------------------------------------------------------------------------------      	  
!---------------------------------------------------------------------------------------------
!-----Random Number Generator taken from Num. Rec. -------------------------------------------    

     FUNCTION ran3(idum)
      INTEGER*4 idum
      INTEGER*4 MBIG,MSEED,MZ
!      REAL MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER*4 i,iff,ii,inext,inextp,k
      INTEGER*4 mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END

!__________________________________________________________________________________________________
! Subr: Generating files
       SUBROUTINE addnumtostring(string,number)
       implicit none
       integer*4 i,strlen,number,nodig,num,snum
       character*(*) string
       snum=number
       do i=len(string),1,-1
       if (string(i:i) .ne. ' ') goto 10
       enddo
   10  strlen=i
       nodig=int(log10(1.0*snum+0.1))+1
       do i=nodig,1,-1
       num=snum/10**(i-1)
       string(strlen+1:strlen+1)=char(48+num)
       strlen=strlen+1
       snum=snum-num*10**(i-1)
       enddo
       return
       end
       integer*4 function strlen(string)
       character*(*) string
       integer*4 i
       do i=len(string),1,-1
       if (string(i:i).ne. ' ') goto 100
       enddo ! i
 100   strlen=i
       end
	
	
