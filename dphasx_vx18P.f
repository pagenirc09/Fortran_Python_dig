c  x-z plane
c  larger number of particles, az is remmoved for thin jet 
c        September 7, 2002
c  Figure is rotated 90 degree cw with trimmed edges
c  contour: Ne: dense, Ni: densi, Jy: ajy, Eke: denqe, Eki: denqi
c           Btot: btot, Ey: aey, Etot: etot
c           GAMMA: vgam
c  arrows:  Bzx, emyx, emyz, Ezx; elyx, elyz
c           Vezx: fvlex, fvley, Vizx; fvlix, fvliy
c
c Contours with color for number, energy, flux, velocity, density
c with magnetic and electric fields
c  
c     parameter (Nproc=10000)
c     parameter (Nproc=3136)

      parameter (Nproc=4356)
c     parameter (Npx=80,Npy=2,Npz=2)
c     parameter (mxs=2005,mys=1005,mzs=1005)
      parameter (mxs=1285,mys=797,mzs=797)

      parameter (nptl=2500000)
      parameter (mb=nptl,mj=1900000)

c     parameter (nptl=6000000)
c     parameter (mb=nptl,mj=9000000)

c     parameter (mg = mb*Nproc,mjg = mj*Nproc)
c     integer*8,parameter:: mg=mb*Nproc
c     integer*8,parameter:: mjg=mj*Nproc
c Ken attention
      parameter (mdiag=3000)
c     parameter (mb=mdiag,mj=mdiag)
      integer*4,parameter:: mg=mdiag*Nproc
      integer*4,parameter:: mjg=mdiag*Nproc
c     integer*4,parameter:: mjg=mdiag*Nproc

c     integer*8 iong, lecg, ionjg, lecjg, ih
c     integer*8 iong, lecg, ih
      integer*4 iong, lecg, ih
      integer*4 ionjg, lecjg, ih4
c     integer*4 ionjg, lecjg, ih
      integer is
      parameter (mlpar=620000)
      parameter (mlpars=620000)
c for 1d
      parameter (mx1=mxs-30,nl=3)
c     parameter (mx1=400,nl=3)

      parameter (imax=32767)
c  attention 90 dgree rotation
      parameter (mxx=mxs-12,myy=mys-20,mzz=mzs-20)
c attention
c     parameter (mxr1=52,   myr1=27,   mzr1=27)
c     parameter (mxr=53,  myr=28,   mzr=28)
c     parameter (mxq=28,   myq=53,   mzq=53)
c     parameter (mxr1=50,   myr1=25,   mzr1=25)
c     parameter (mxr=50,  myr=25,   mzr=25)
c     parameter (mxq=50,   myq=25,   mzq=25)
      parameter (mxr1=54,   myr1=8,   mzr1=8)
      parameter (mxr=54,  myr=8,   mzr=8)
      parameter (mxq=54,  myq=8,   mzq=8)
c     real*8 me,mi,mi2

c     dimension ixx(mb),iyy(mb),izz(mb)
c  ambient ions and electrons positions and velocities arrays
c     real*4 xi(mb),yi(mb),zi(mb),ui(mb),vi(mb),wi(mb)
c     real*4 xe(mb),ye(mb),ze(mb),ue(mb),ve(mb),we(mb)

c  jet ions and electrons positions and velocities arrays
c     real*4 xij(mj),yij(mj),zij(mj),uij(mj),vij(mj),wij(mj)
c     real*4 xej(mj),yej(mj),zej(mj),uej(mj),vej(mj),wej(mj)

c  global field arrays
c     real*4 bxx,byy,bzz,exx,eyy,ezz
c     real*4 bxx(mxs,mys,mzs),byy(mxs,mys,mzs),bzz(mxs,mys,mzs)
c     real*4 exx(mxs,mys,mzs),eyy(mxs,mys,mzs),ezz(mxs,mys,mzs)
c     dimension bxx(mx,my,mz),byy(mx,my,mz),bzz(mx,my,mz)
c     dimension exx(mx,my,mz),eyy(mx,my,mz),ezz(mx,my,mz)

c  global particle arrays
c     real*4 xig,yig,zig,uig,vig,wig
c     real*4 xeg,yeg,zeg,ueg,veg,weg
c     real*4 xijg,yijg,zijg,uijg,vijg,wijg
c     real*4 xejg,yejg,zejg,uejg,vejg,wejg
      real*4 :: xig(mg),yig(mg),zig(mg),uig(mg),vig(mg),wig(mg)
      real*4 :: xeg(mg),yeg(mg),zeg(mg),ueg(mg),veg(mg),weg(mg)
      real*4 :: xijg(mjg),yijg(mjg),zijg(mjg)
      real*4 :: uijg(mjg),vijg(mjg),wijg(mjg)
      real*4 :: xejg(mjg),yejg(mjg),zejg(mjg)
      real*4 :: uejg(mjg),vejg(mjg),wejg(mjg)
c     dimension xig(mg),yig(mg),zig(mg),uig(mg),vig(mg),wig(mg)
c     dimension xeg(mg),yeg(mg),zeg(mg),ueg(mg),veg(mg),weg(mg)
c     dimension xijg(mjg),yijg(mjg),zijg(mjg)
c     dimension uijg(mjg),vijg(mjg),wijg(mjg)
c     dimension xejg(mjg),yejg(mjg),zejg(mjg)
c     dimension uejg(mjg),vejg(mjg),wejg(mjg)

      real*4 ::
     &   c4,DT4,qi4,qe4,qmi4,qme4,vithml4,vethml4,vijet4,vejet4,
     &   vithmj4,vethmj4,refli4,refle4,rselect4,xj04,yj04,zj04,
     &   b0x4,b0y4,e0z4,
     &         dlxj4,dlyj4,dlzj4

c
      real*4 :: xece(mlpar), yece(mlpar), zece(mlpar),
     1          vxce(mlpar), vyce(mlpar), vzce(mlpar)
     1    ,vegm(mlpar)
      real*4 :: xeci(mlpar), yeci(mlpar), zeci(mlpar),
     1          vxci(mlpar), vyci(mlpar), vzci(mlpar)
     1    ,vigm(mlpar)
      real*4 :: zsp(mxx),zjb(mxx),zjg(mxx),zjt(mxx)
      real*4 :: xp(mlpars),yp(mlpars)

      dimension gridr(4),gridm(4),gride(4),bnum(2),bnu1(2) 
      dimension xd(3),yd(3)                                             
c Display size
c     parameter (mxr=70   , myr= 30   , mzr= 30   )
c     dimension ex(mxs,mys,mzs),ey(mxs,mys,mzs),ez(mxs,mys,mzs),
c    &bx(mxs,mys,mzs),by(mxs,mys,mzs),bz(mxs,mys,mzs)
      dimension sm(27),ms(27),table(64)
c
      dimension delzx(mxx,myy), delzy(mxx,myy),
     &          demzx(mxx,myy), demzy(mxx,myy),
     &          delyx(mxx,mzz), delyz(mxx,mzz),
     &          demyx(mxx,mzz), demyz(mxx,mzz),
     &          delxy(myy,mzz), delxz(myy,mzz),
     &          demxy(myy,mzz), demxz(myy,mzz)
c
      dimension elzx(mxr1,myr1), elzy(mxr1,myr1),
     &          emzx(mxr1,myr1), emzy(mxr1,myr1),
     &          elyx(mxr1,mzr1), elyz(mxr1,mzr1),
     &          emyx(mxr1,mzr1), emyz(mxr1,mzr1),
     &          elxy(myq,mzq), elxz(myq,mzq),
     &          emxy(myq,mzq), emxz(myq,mzq)
c for 1d
      dimension xsp(mx1),del(nl,mx1),din(nl,mx1),cha(nl,mx1),
     1          ex1(nl,mx1),ey1(nl,mx1),ez1(nl,mx1),
     1          bx1(nl,mx1),by1(nl,mx1),bz1(nl,mx1),
     1          cx1(nl,mx1),cy1(nl,mx1),cz1(nl,mx1),
     1          fxx(mx1)
      
c
c     data gride/0.10,0.75,0.35,0.593/
c     data gridr/0.275,0.70,0.15,0.95/
c     data gridr/0.10,0.75,0.300,0.460/
      data gridr/0.10,0.75,0.300,0.625/
      data gride/0.10,0.75,0.25,0.7/,bnum/0.030,0.030/
      data gridm/0.25,0.95,0.25,0.5478/,bnu1/0.030,0.030/
      data xd/11.,12.,15./, yd/1.0,2.0,3.0/                             
c                                                                       
c     dimension x(nptl),y(nptl),z(nptl),
c    1          u(nptl),v(nptl),w(nptl)
c     dimension     mcont(mbeam),kcont(mback),
c    1              mcon2(mbeam),kcon2(mback),
c    1              mcone(mlpar),kconi(mlpar)
      COMMON /VEC1/   ASH        ,EXT        ,ICTRFG     ,ILAB       ,
     +                IOFFD      ,IOFFM      ,ISX        ,ISY        ,
     +                RMN        ,RMX        ,SIDE       ,SIZE       ,
     +                XLT        ,YBT        ,ZMN        ,ZMX
      common /srfip1/ ifr    ,istp   ,irots  ,idrx ,idry,               
     1                idrz   ,iupper ,iskirt ,ncla ,theta,              
     2                hskirt ,chi    ,clo    ,cinc ,ispval              
c------------------------------------------------------------
      character num*3,num1,num2,num3
      character step*6,step1,step2,step3,step4,step5
      character st01,st02,st03,st0*3,st*4,hyph
      character dchar1*6,dchar2*4,file1*10,file2*8
      character strpd*6
      character message2*6

      character message*29                                              
      character message1*22
      character*1 ittmp(21)
      open(6,file='dwriteot',status='unknown',form='formatted')

      call opngks
      call gsclip (0)

c specify "nst" - which data need to be processed
c      open(unit=2, file='nstn', status='old')

c      read(2,*) nst
c      do ix=11,8,-4
       do ix=10,10
c      do ix=7,7
c      do ix=6,2,-1

       nst=ix
      
c      nst = 5
      
       c = 1.0

c      strpd = 'partd_'
       hyph='_'

       st01 = char(int(nst/100.)+48)
       st02 = char(int((nst-int(nst/100.)*100)/10.)+48)
       st03 = char(int(nst-int(nst/10.)*10)+48)
       st0 = st01//st02//st03
       st = hyph//st0

       write(6,*) 'before reading'  

        dchar2='ions'
        file2= dchar2//st
cheung  
        open(17,file=file2,form='unformatted')
        read(17) iong
        read(17) xig(1:iong)
        read(17) yig(1:iong)
        read(17) zig(1:iong)
        read(17) uig(1:iong)
        read(17) vig(1:iong)
        read(17) wig(1:iong)

        close(17)
c      write(6,*) 'after reading 17,iong=',iong
c      print *, 'after reading 17,iong=',iong

        dchar2='lecs'
        file2= dchar2//st

cheung
        open(18,file=file2,form='unformatted')
        read(18) lecg
        read(18) xeg(1:lecg)
        read(18) yeg(1:lecg)
        read(18) zeg(1:lecg)
        read(18) ueg(1:lecg)
        read(18) veg(1:lecg)
        read(18) weg(1:lecg)

        close(18)
       write(6,*) 'after reading 18,lecg=',lecg
       write(6,*) 'xeg(50),ueg(50)=',xeg(50),ueg(50)      

        dchar2='ionj'
        file2= dchar2//st

cheung
        open(19,file=file2,form='unformatted')
        read(19) ionjg
        read(19) xijg(1:ionjg)
        read(19) yijg(1:ionjg)
        read(19) zijg(1:ionjg)
        read(19) uijg(1:ionjg)
        read(19) vijg(1:ionjg)
        read(19) wijg(1:ionjg)
        
        close(19) 
c      write(6,*) 'after reading 19,ionjg=',ionjg
c      print *, 'after reading 19,ionjg=',ionjg
        
        dchar2='lecj'
        file2= dchar2//st
        
cheung  
        open(20,file=file2,form='unformatted')
        read(20) lecjg
        read(20) xejg(1:lecjg)
        read(20) yejg(1:lecjg)
        read(20) zejg(1:lecjg)
        read(20) uejg(1:lecjg)
        read(20) vejg(1:lecjg)
        read(20) wejg(1:lecjg)

        close(20)
c      write(6,*) 'after reading 20,lecjg=',lecjg
c      write(6,*) 'xejg(50),uejg(50)=',xejg(50),uejg(50)      
c      print *, 'after reading 20,lecjg=',lecjg


c       message2='kenout'

c       dchar1='kenout'
c       file1= dchar1//st
c       open(21,file=file1,form='unformatted')

c       open(21,file=message2,form='unformatted')
c        read(21) nstep
c    &  ,c,DT,qi,qe,qmi,qme,vithml,vethml,vijet,vejet,
c    &         vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c new Jacek
c    &          b0y,e0z
c       read(21) nstep,
c    &   c4,DT4,qi4,qe4,qmi4,qme4,vithml4,vethml4,vijet4,vejet4,
c    &   vithmj4,vethmj4,refli4,refle4,rselect4,xj04,yj04,zj04,
c    &   b0x4,b0y4,e0z4

c       read(21) nstep
c    &  ,c,DT,qi,qe,qmi,qme,vithml,vethml,vijet,vejet,
c    & vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,
c    &  vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c    &  vithmj,vethmj


c,DT,qi,qe,mi,me,qmi,qme,vithml,vethml,vijet,vejet,
c    &         vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c new Jacek
c    &          b0y,e0z
c       close(21)
c      write(6,*) 'after reading 21,nstep=',nstep


c     write(6,*) ions, lecs,maxptl,maxhlf
      je=47
c     ke=48
c     ie=70
      in=43
      jn=43
      kn=43
      imid=in
      jmid=jn
      kmid=kn
      t = nstep
      xmax=1183.
c
c attention        x-z plane Ex,z, Bx,z
c center of jet
       yjc = 381.
       zjc = 381.
       rjt = 100.
       is = -5
      rjt2 = rjt**2

      sran = 0.03*mlpar/lecjg
c     sran = 0.5
c for positrons (ions)
      i=0
      do ih4=1,ionjg
      aran =ran1(is)
      if(aran.gt.sran) go to 215
      i=i+1
      xece(i)=xijg(ih4)
      yece(i)=yijg(ih4)
      zece(i)=zijg(ih4)
      vxce(i)=uijg(ih4)/c
      vyce(i)=vijg(ih4)/c
      vzce(i)=wijg(ih4)/c
      vegm(i)=1.0/sqrt(1.0
     1    -(vxce(i)**2+vyce(i)**2+vzce(i)**2))
  215 continue
      enddo
       write(6,*) 'xece(50),vxce(50)=',xece(50),vxce(50)     
      mmae=i

       write(6,*) 'after 215,mmae=,vxce(i)=',i,vxce(i)
      sran = 0.03*mlpar/float(lecg)
c     sran = 0.01
      is = -2
      i=0
c for positrons (ions)
      do ih=1,iong 
c       rjet = (yeg(ih)-yjc)**2 +(zeg(ih)-zjc)**2
c      if(rjet.ge.rjt2) go to 816
      if(xig(ih).gt.xmax) go to 816
      aran =ran1(is)
      if(aran.gt.sran) go to 816
      i=i+1
      xeci(i)=xig(ih)
      yeci(i)=yig(ih)
      zeci(i)=zig(ih)
      vxci(i)=uig(ih)/c
      vyci(i)=vig(ih)/c
      vzci(i)=wig(ih)/c
      vigm(i)=1.0/sqrt(1.0
     1    -(vxci(i)**2+vyci(i)**2+vzci(i)**2))
  816 continue
      enddo
       write(6,*) 'xeci(50),vxci(50)=',xeci(50),vxci(50)      
      kmai=i

       write(6,*) 'after 816,kmai=',kmai
      xmin=3.0
c     xmax=2003.
c     xmax=1283.
      xmax=1183.
c     ymin=-0.05
c     ymax=0.4
c     ymax=0.2
c     ymax=0.25
      ymin=-10.0
c     ymin=-20.0
c     ymax=100.00
      ymax=50.00
c     ymax=30.00
      isc=1
      call dscale(gridr,bnu1,xmin,xmax,ymin,ymax,
     1            0.04,-0.03,'Y$','Z$',message,isc)

       write(6,*) 'after dscale,kmai=',kmai

      ratn = 0.01
c     ratn = 0.3
      ratn = 1.00
      mr = 0
      gmax= 50.
      do m = 1, mmae
c     arat = rand()
c     arat4 = arat
c     if(arat4.gt.ratn) go to 11
      mr = mr +1
      xp(mr) = xece(m)
c     yp(mr) = vxce(m)**2+vyce(m)**2
c     yp(mr) = sqrt(vxce(m)**2+vyce(m)**2)
      gam    = 1/sqrt(1.0-(vxce(m)**2+vyce(m)**2
     1          +vzce(m)**2))
c     gam    = 1/sqrt(1.0- vzce(m)**2))
      yp(mr) = vxce(m)*gam
      if(yp(mr).ge.ymax) yp(mr)=ymax
      if(yp(mr).le.ymin) yp(mr)=ymin
c  11 continue
      enddo
       write(6,*) 'after 11,mmae=',mmae
       write(6,*) 'gam, yp(mr)=', gam, yp(mr)
      CALL DFCLRS
c 15: red 6:green  need to be checked
        CALL GSPMCI (6)
        CALL GSTXCI (6)
        CALL GSPLCI (6)
      mam = mr
c       CALL DFCLRS
c       CALL GSPMCI (0)
c       CALL GSTXCI (1)
c       CALL GSPLCI (0)
      call points(xp,yp,mam,-1,0)
c     call points(xp,yp,mam,0,0)

       write(6,*) 'after points,mam=',mam
c ambient electrons
      ratn = 1.00
      mr = 0
      gmax= 50.
      do m = 1, kmai
c     arat = rand()
c     arat4 = arat
c     if(arat4.gt.ratn) go to 11
      mr = mr +1
      xp(mr) = xeci(m)
c     yp(mr) = vxce(m)**2+vyci(m)**2
c     yp(mr) = sqrt(vxce(m)**2+vyce(m)**2)
      gam    = 1/sqrt(1.0-(vxci(m)**2+vyci(m)**2
     1          +vzci(m)**2))
c     gam    = 1/sqrt(1.0- vzce(m)**2))
      yp(mr) = vxci(m)*gam
      if(yp(mr).ge.ymax) yp(mr)=ymax
      if(yp(mr).le.ymin) yp(mr)=ymin
c  11 continue
      enddo
      CALL DFCLRS
        CALL GSPMCI (2)
        CALL GSTXCI (2)
        CALL GSPLCI (2)
      mam = mr
c       CALL DFCLRS
c       CALL GSPMCI (0)
c       CALL GSTXCI (1)
c       CALL GSPLCI (0)
      call points(xp,yp,mam,-1,0)
c     call points(xp,yp,mam,0,0)

       write(6,*) 'after pointsg,mam=',mam

        CALL CPCLRS
c       CALL GSPMCI (0)
        CALL GSTXCI (2)
        CALL GSPLCI (2)
      write(message,27) t
c       CALL LABTOP (message,.018)
      xpob=0.425
      ypob=0.25
      xpol=0.05
      ypol=0.475
      size=0.018
c   for z - gamma
      ichs=1 
c   for x - z
c     ichs=2
      call eplch(xpob,ypob,xpol,ypol,size,ichs)
      
      call frame

      end do
       write(6,*) 'after frame,mam=',mam
  137 format('ELECTRON DEN   T=',f6.1)

  110 format('ELE DENSITY    T=',f8.1)
  111 format('BF ENERGY      T=',f8.1)
  112 format('Z-CURRNET DEN  T=',f6.1)
  113 format('Y-MAGNE FIELD  T=',f8.1)
c 113 format('Z-ELECT FIELD  T=',f6.1)
  114 format('BTOT    (MAG)  T=',f6.1)
  115 format('ETOT    (ELE)  T=',f6.1)
  116 format('GAMMAE  (EFLU) T=',f6.1)
  117 format('EY      (ELE)  T=',f6.1)
  118 format('EKIN    (EFLU) T=',f6.1)
  119 format('IKIN    (IFLU) T=',f6.1)
  214 format('CUR DEN (MAG)  T=',f6.1)
  275 format('CUR DEN (FLUX) T=',f6.1)
  216 format('CUR DEN (CUR)  T=',f6.1)
  223 format('ELE ENE (MAG)  T=',f6.1)
  224 format('ION ENE (MAG)  T=',f6.1)
  323 format('ELE ENE (FLUX) T=',f6.1)
  324 format('ION ENE (FLUX) T=',f6.1)
   25 format('X-Y SPACE (ION) T=',F6.0,'$')
   26 format('X-Y SPACE (ELE) T=',F8.0,'$')
   27 format('X-VX  SPACE     T=',F8.0,'$')

      call clsgks
      stop
      end

      subroutine Vol_Weighting(mxs,mys,mzs,lecg,rhe,eke,ueg,veg,weg,
     &                         xeg,yeg,zeg,fex,fey,fez)
      implicit none

c     integer mxs,mys,mzs,mg,lecg
      integer mxs,mys,mzs,mg
      integer*8 lecg,n0
      real*4 xeg(lecg),yeg(lecg),zeg(lecg)
      real*4 ueg(lecg),veg(lecg),weg(lecg)
      real rhe(mxs,mys,mzs),eke(mxs,mys,mzs)
      real fex(mxs,mys,mzs),fey(mxs,mys,mzs),fez(mxs,mys,mzs)
      real*4 dx,dy,dz,cx,cy,cz
      real*4 sl,sr,sn,s     
      integer i,j,k,l,m,n

C  Volume weighting:

c     mhlec = maxhlf +lecs
      do 3 n0=1, lecg
c     if(n0.gt.mhlec.and.n0.lt.lecm) go to 3

      i=xeg(n0)
      dx=xeg(n0)-i
      cx=1.-dx
      j=yeg(n0)
      dy=yeg(n0)-j
      cy=1.-dy
      k=zeg(n0)
      dz=zeg(n0)-k
      cz=1.-dz
C  Smoothing with the (.25,.5,.25) profile in each dimension:
      sl=.5
      do 187 l=-1,1
      sl=.75-sl
      sr=.5
      do 187 m=-1,1
      sr=.75-sr
      sn=.5
      do 187 n=-1,1
      sn=.75-sn
      s=sl*sr*sn
      rhe(i+l  ,j+m  ,k+n  )=rhe(i+l  ,j+m  ,k+n  )+s*cx*cy*cz
      rhe(i+l+1,j+m  ,k+n  )=rhe(i+l+1,j+m  ,k+n  )+s*dx*cy*cz
      rhe(i+l  ,j+m+1,k+n  )=rhe(i+l  ,j+m+1,k+n  )+s*cx*dy*cz
      rhe(i+l  ,j+m  ,k+n+1)=rhe(i+l  ,j+m  ,k+n+1)+s*cx*cy*dz
      rhe(i+l  ,j+m+1,k+n+1)=rhe(i+l  ,j+m+1,k+n+1)+s*cx*dy*dz
      rhe(i+l+1,j+m  ,k+n+1)=rhe(i+l+1,j+m  ,k+n+1)+s*dx*cy*dz
      rhe(i+l+1,j+m+1,k+n  )=rhe(i+l+1,j+m+1,k+n  )+s*dx*dy*cz
      rhe(i+l+1,j+m+1,k+n+1)=rhe(i+l+1,j+m+1,k+n+1)+s*dx*dy*dz

c
      eke(i+l  ,j+m  ,k+n  )=eke(i+l  ,j+m  ,k+n  )
     1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*cx*cy*cz
      eke(i+l+1,j+m  ,k+n  )=eke(i+l+1,j+m  ,k+n  )
     1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*dx*cy*cz
      eke(i+l  ,j+m+1,k+n  )=eke(i+l  ,j+m+1,k+n  )
     1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*cx*dy*cz
      eke(i+l  ,j+m  ,k+n+1)=eke(i+l  ,j+m  ,k+n+1)
     1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*cx*cy*dz
      eke(i+l  ,j+m+1,k+n+1)=eke(i+l  ,j+m+1,k+n+1)
     1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*cx*dy*dz
      eke(i+l+1,j+m  ,k+n+1)=eke(i+l+1,j+m  ,k+n+1)
     1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*dx*cy*dz
      eke(i+l+1,j+m+1,k+n  )=eke(i+l+1,j+m+1,k+n  )
     1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*dx*dy*cz
      eke(i+l+1,j+m+1,k+n+1)=eke(i+l+1,j+m+1,k+n+1)
     1 +(ueg(n0)*ueg(n0)+veg(n0)*veg(n0)+weg(n0)*weg(n0))*s*dx*dy*dz
c
      fex(i+l  ,j+m  ,k+n  )=fex(i+l  ,j+m  ,k+n  )
     1        +ueg(n0)*s*cx*cy*cz
      fex(i+l+1,j+m  ,k+n  )=fex(i+l+1,j+m  ,k+n  )
     1        +ueg(n0)*s*dx*cy*cz
      fex(i+l  ,j+m+1,k+n  )=fex(i+l  ,j+m+1,k+n  )
     1        +ueg(n0)*s*cx*dy*cz
      fex(i+l  ,j+m  ,k+n+1)=fex(i+l  ,j+m  ,k+n+1)
     1        +ueg(n0)*s*cx*cy*dz
      fex(i+l  ,j+m+1,k+n+1)=fex(i+l  ,j+m+1,k+n+1)
     1        +ueg(n0)*s*cx*dy*dz
      fex(i+l+1,j+m  ,k+n+1)=fex(i+l+1,j+m  ,k+n+1)
     1        +ueg(n0)*s*dx*cy*dz
      fex(i+l+1,j+m+1,k+n  )=fex(i+l+1,j+m+1,k+n  )
     1        +ueg(n0)*s*dx*dy*cz
      fex(i+l+1,j+m+1,k+n+1)=fex(i+l+1,j+m+1,k+n+1)
     1        +ueg(n0)*s*dx*dy*dz
c
      fey(i+l  ,j+m  ,k+n  )=fey(i+l  ,j+m  ,k+n  )
     1        +veg(n0)*s*cx*cy*cz
      fey(i+l+1,j+m  ,k+n  )=fey(i+l+1,j+m  ,k+n  )
     1        +veg(n0)*s*dx*cy*cz
      fey(i+l  ,j+m+1,k+n  )=fey(i+l  ,j+m+1,k+n  )
     1        +veg(n0)*s*cx*dy*cz
      fey(i+l  ,j+m  ,k+n+1)=fey(i+l  ,j+m  ,k+n+1)
     1        +veg(n0)*s*cx*cy*dz
      fey(i+l  ,j+m+1,k+n+1)=fey(i+l  ,j+m+1,k+n+1)
     1        +veg(n0)*s*cx*dy*dz
      fey(i+l+1,j+m  ,k+n+1)=fey(i+l+1,j+m  ,k+n+1)
     1        +veg(n0)*s*dx*cy*dz
      fey(i+l+1,j+m+1,k+n  )=fey(i+l+1,j+m+1,k+n  )
     1        +veg(n0)*s*dx*dy*cz
      fey(i+l+1,j+m+1,k+n+1)=fey(i+l+1,j+m+1,k+n+1)
     1        +veg(n0)*s*dx*dy*dz
c
      fez(i+l  ,j+m  ,k+n  )=fez(i+l  ,j+m  ,k+n  )
     1        +weg(n0)*s*cx*cy*cz
      fez(i+l+1,j+m  ,k+n  )=fez(i+l+1,j+m  ,k+n  )
     1        +weg(n0)*s*dx*cy*cz
      fez(i+l  ,j+m+1,k+n  )=fez(i+l  ,j+m+1,k+n  )
     1        +weg(n0)*s*cx*dy*cz
      fez(i+l  ,j+m  ,k+n+1)=fez(i+l  ,j+m  ,k+n+1)
     1        +weg(n0)*s*cx*cy*dz
      fez(i+l  ,j+m+1,k+n+1)=fez(i+l  ,j+m+1,k+n+1)
     1        +weg(n0)*s*cx*dy*dz
      fez(i+l+1,j+m  ,k+n+1)=fez(i+l+1,j+m  ,k+n+1)
     1        +weg(n0)*s*dx*cy*dz
      fez(i+l+1,j+m+1,k+n  )=fez(i+l+1,j+m+1,k+n  )
     1        +weg(n0)*s*dx*dy*cz
      fez(i+l+1,j+m+1,k+n+1)=fez(i+l+1,j+m+1,k+n+1)
     1        +weg(n0)*s*dx*dy*dz
  187 continue
    3 continue
      return
      end subroutine Vol_Weighting

c
c  Subroutine for color contour plot
c
c       subroutine colorcn(ZDAT,nx,ny,ZMIN,ZMAX,grid,message,ici,ilb)
        subroutine colorcn(ZDAT,nx,ny,ZMIN,ZMAX,grid,message,ici)
C
C Declare required data arrays and workspace arrays.
C
c       parameter (nx=43,ny=19)
c       parameter (nspc=1000,nwsp=12000)
c       parameter (nspc=6000,nwsp=82000)
c       parameter (nspc=36000,nwsp=482000)
        parameter (nspc=66000,nwsp=882000)
        dimension grid(4)
        DIMENSION ZDAT(nx,ny),RWRK(nspc),IWRK(nspc),IAMA(nwsp) 
        DIMENSION IASF(13)
        DIMENSION XCRA(nspc),YCRA(nspc)
        DIMENSION IAIA(10),IGIA(10)
C
C Declare arrays to hold the list of indices and the list of labels
C required by the label-bar routine.
C
        DIMENSION LIND(14)
        CHARACTER*10 LLBS(15)
        character*(*) message
C
C Declare the routine which will color the areas.
C
        EXTERNAL COLRAM
C
C Define an array for GKS aspect source flags.
C
        DATA IASF / 13*1 /
C
C Define the list of indices required by the label-bar routine.
C
        DATA LIND / 2,3,4,5,6,7,8,9,10,11,12,13,14,15 /
C
C Open GKS.
C
cm      CALL OPNGKS
C
C Turn off the clipping indicator.
C
cm      CALL GSCLIP (0)
c 
c Put a label at the top of theframe
c
cp      call set (0.,1.,0.,1.,0.,1.,0.,1.,1)
cp      call wtstr (0.5,.9765,message,2.0,0)
C
C Set all aspect source flags to "individual".
C
        CALL GSASF (IASF)
C
C Force solid fill.
C
        CALL GSFAIS (1)
C
C Define color indices.
C
c choosing color indices
c
        if(ici.eq.1) CALL DFCLRS
        if(ici.eq.2) CALL CPCLRS
C
C Generate an array of test data.
C
cx      CALL GENDAT (ZDAT,23,23,14,20,20,-136.148,451.834)
C
C Get the current elapsed time, in seconds.
C
c       TIME=SECOND(DUMI)
C
C Force the plot into the upper portion of the frame.
c       CALL CPSETR ('VPB - VIEWPORT BOTTOM',.25)
C
C Force the plot into the left side of the frame.
C
c       CALL CPSETR ('VPR - VIEWPORT RIGHT',.80)
C
       xvl=grid(1)
       xvr=grid(2)
       yvb=grid(3)
       yvt=grid(4)
c
        CALL CPSETR ('VPL - VIEWPORT LEFT',xvl)
        CALL CPSETR ('VPR - VIEWPORT RIGHT',xvr)
        CALL CPSETR ('VPB - VIEWPORT BOTTOM',yvb)
        CALL CPSETR ('VPT - VIEWPORT TOP',yvt)
        CALL CPSETR ('VPS - VIEWPORT SHAPE',0.)
C
C Disallow the trimming of trailing zeroes.
C
        CALL CPSETI ('NOF - NUMERIC OMISSION FLAGS',0)
C
C Tell CONPACK to use 13 contour levels, splitting the range into 14
C equal bands, one for each of the 14 colors available.
C
        CALL CPSETI ('CLS - CONTOUR LEVEL SELECTOR',-13)
C
C Initialize the drawing of the contour plot.
C
c       CALL CPRECT (ZDAT,23,23,14,RWRK,1000,IWRK,1000)
        CALL CPRECT (ZDAT,nx,nx,ny,RWRK,nspc,IWRK,nspc)
C
C Initialize the area map and put the contour lines into it.
C
        CALL ARINAM (IAMA,nwsp)  
        CALL CPCLAM (ZDAT,RWRK,IWRK,IAMA)
C
C Color the map.
C
        CALL ARSCAM (IAMA,XCRA,YCRA,nspc,IAIA,IGIA,10,COLRAM)
C
C Put black contour lines over the colored map.
C
c       CALL GSPLCI (0)
c       CALL CPCLDR (ZDAT,RWRK,IWRK)
c       CALL GSPLCI (1)
C
C Draw a label bar for the plot, relating colors to values.
C
c       CALL CPGETR ('ZMN',ZMIN)
c       CALL CPGETR ('ZMX',ZMAX)
C
        DO 102 I=1,15
          CALL CPSETR ('ZDV - Z DATA VALUE',
     +                 ZMIN+REAL(I-1)*(ZMAX-ZMIN)/14.)
          CALL CPGETC ('ZDV - Z DATA VALUE',LLBS(I))
  102   CONTINUE
C
c       CALL LBSETI ('CBL - COLOR OF BOX LINES',0)
c       CALL LBLBAR (0,.05,.95,.15,.25,14,1.,.5,LIND,0,LLBS,15,1)
c
c   Drow a color bar for the plot.
c
        CALL GSTXCI (0)
        CALL LBSETI ('CBL - COLOR OF BOX LINES',2)
        call pcsetr ('CS - CONSTANT SPACING FLAG',1.25)
c       if(grid(4).gt.0.63) then
c       CALL LBLBAR (1,.79,.89,.12,.72,14,.333,1.,LIND,0,LLBS,15,1)
c       CALL LBLBAR (1,.79,.89,.12,.72,14,.333,1.,LIND,0,LLBS,15,1)
c       else
c       CALL LBLBAR (1,.79,.89,.20,.80,14,.333,1.,LIND,0,LLBS,15,1)
        CALL LBLBAR (0,.10,.95,.23,.30,14,1.,0.5,LIND,0,LLBS,15,1)
c       endif
        CALL GSTXCI (1)
C
C Compute and print statistics for the plot, label it, and put a
C boundary line at the edge of the plotter frame.
C
c       CALL CAPSAP ('EXAMPLE 7',TIME,IAMA,12000)
c       CALL CAPSAP ('EXAMPLE 7',TIME,IAMA,nwsp) 
c       CALL LABTOP ('EXAMPLE 7',.017)
c       CALL GSTXCI (0)
c      if(ilb.eq.1)
c    .  CALL LABTOP (message,.018)
        CALL LABTOP (message,.018)
c      if(ilb.eq.2)
c    .  CALL LABLFT (message,0.018)
c       CALL GSTXCI (1)
c       CALL BNDARY
C
C Advance the frame.
C
c       CALL FRAME
C
C Close GKS.
C
cm      CALL CLSGKS
C
C Done.
C
c       STOP
        return
C
       END
       SUBROUTINE COLRAM (XCRA,YCRA,NCRA,IAIA,IGIA,NAIA)
C
        DIMENSION XCRA(*),YCRA(*),IAIA(*),IGIA(*)
C
C The arrays XCRA and YCRA, for indices 1 to NCRA, contain the X and Y
C coordinates of points defining a polygon.  The area identifiers in
C the array IAIA, each with an associated group identifier in the array
C IGIA, tell us whether the polygon is to be color-filled or not.
C
C
C Assume the polygon will be filled until we find otherwise.
C
        IFLL=1
C
C If any of the area identifiers is negative, don't fill the polygon.
C
        DO 101 I=1,NAIA
          IF (IAIA(I).LT.0) IFLL=0
  101   CONTINUE
C
C Otherwise, fill the polygon in the color implied by its area
C identifier relative to edge group 3 (the contour-line group).
C
        IF (IFLL.NE.0) THEN
          IFLL=0
          DO 102 I=1,NAIA
            IF (IGIA(I).EQ.3) IFLL=IAIA(I)
  102     CONTINUE
          IF (IFLL.GT.0.AND.IFLL.LT.15) THEN
            CALL GSFACI (IFLL+1)
            CALL GFA (NCRA-1,XCRA,YCRA)
          END IF
        END IF
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE DFCLRS
C
C Define a set of RGB color triples for colors 1 through 15.
C
        DIMENSION RGBV(3,15)
C
C Define the RGB color triples needed below.
C
c  Rainbow color bar
        DATA RGBV / 1.00 , 1.00 , 1.00 ,
c    +              0.70 , 0.70 , 0.70 ,
c    +              0.75 , 0.50 , 1.00 ,
c    +              0.50 , 0.00 , 1.00 ,
c   added
     +              0.00 , 0.00 , 1.00 ,
     +              0.00 , 0.50 , 1.00 ,
     +              0.00 , 1.00 , 1.00 ,
     +              0.00 , 1.00 , 0.60 ,
     +              0.00 , 1.00 , 0.00 ,
c   added
     +              0.20 , 1.00 , 0.00 ,
     +              0.50 , 1.00 , 0.00 ,
     +              0.75 , 1.00 , 0.00 ,
     +              1.00 , 1.00 , 0.00 ,
c    +              1.00 , 0.75 , 0.00 ,
c    +              1.00 , 0.38 , 0.38 ,
     +              1.00 , 0.80 , 0.00 ,
     +              1.00 , 0.55 , 0.25 ,
     +              1.00 , 0.38 , 0.00 ,
     +              1.00 , 0.20 , 0.00 ,
     +              1.00 , 0.00 , 0.00 /
C
C
C
C
C Define 16 different color indices, for indices 0 through 15.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
        CALL GSCR (1,0,0.,0.,0.)
C
        DO 101 I=1,15
          CALL GSCR (1,I,RGBV(1,I),RGBV(2,I),RGBV(3,I))
  101   CONTINUE
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE BFCLRS
C
C Define a set of RGB color triples for colors 1 through 15.
C
        DIMENSION RGBV(3,15)
C
C Define the RGB color triples needed below.
C
c  Rainbow color bar
        DATA RGBV / 1.00 , 0.00 , 0.00 ,
c    +              0.70 , 0.70 , 0.70 ,
c    +              0.75 , 0.50 , 1.00 ,
c    +              0.50 , 0.00 , 1.00 ,
c   added
     +              0.00 , 0.00 , 1.00 ,
     +              0.00 , 0.50 , 1.00 ,
     +              0.00 , 1.00 , 1.00 ,
     +              0.00 , 1.00 , 0.60 ,
     +              0.00 , 1.00 , 0.00 ,
c   added
     +              0.20 , 1.00 , 0.00 ,
     +              0.50 , 1.00 , 0.00 ,
     +              0.75 , 1.00 , 0.00 ,
     +              1.00 , 1.00 , 0.00 ,
c    +              1.00 , 0.75 , 0.00 ,
c    +              1.00 , 0.38 , 0.38 ,
     +              1.00 , 0.80 , 0.00 ,
     +              1.00 , 0.55 , 0.25 ,
     +              1.00 , 0.38 , 0.00 ,
     +              1.00 , 0.20 , 0.00 ,
     +              1.00 , 0.00 , 0.00 /
C
C
C
C
C Define 16 different color indices, for indices 0 through 15.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
        CALL GSCR (1,0,0.,0.,0.)
C
        DO 101 I=1,15
          CALL GSCR (1,I,RGBV(1,I),RGBV(2,I),RGBV(3,I))
  101   CONTINUE
C
C Done.
C
        RETURN
C
      END
      SUBROUTINE CPCLRS
C
C Define a set of RGB color triples for colors 1 through 15.
C
        DIMENSION RGBV(3,15)
C
C Define the RGB color triples needed below.
C
c   Black to Red   (through White)
c
        DATA RGBV / 1.00 , 1.00 , 1.00 ,
     +              0.00 , 0.00 , 0.00 ,
     +              0.00 , 0.50 , 0.63 ,
     +              0.00 , 0.75 , 1.00 ,
     +              0.25 , 1.00 , 1.00 ,
     +              0.50 , 1.00 , 1.00 ,
     +              0.70 , 1.00 , 1.00 ,
     +              0.85 , 1.00 , 1.00 ,
c    +              1.00 , 1.00 , 1.00 ,
c    +              1.00 , 1.00 , 1.00 ,
     +              1.00 , 1.00 , 0.85 ,
     +              1.00 , 0.90 , 0.55 ,
     +              1.00 , 0.83 , 0.25 ,
     +              1.00 , 0.50 , 0.12 ,
     +              0.75 , 0.25 , 0.25 ,
     +              0.75 , 0.12 , 0.38 ,
     +              1.00 , 0.00 , 0.00 /
C
C Define 16 different color indices, for indices 0 through 15.  The
C color corresponding to index 0 is black and the color corresponding
C to index 1 is white.
C
        CALL GSCR (1,0,0.,0.,0.)
C
        DO 101 I=1,15
          CALL GSCR (1,I,RGBV(1,I),RGBV(2,I),RGBV(3,I))
  101   CONTINUE
C
C Done.
C
        RETURN
C
      END
c
c=================  LABTOP ======================
c
      SUBROUTINE LABTOP(LABL,SIZE)
c
      CHARACTER*(*) LABL
c
c  Put a label just above the top of the plot The SET call is re-done
c  to allow for the use of fractional coordinates, and the text extent
c  capabilities of the package PLOTCHAR are used to determine the label
c  position.
c
c  To put a slight margin from the top of the frame
      TMAG=0.05
c     
        CALL GETSET (XVPL,XVPR,YVPB,YVPT,XWDL,XWDR,YWDB,YWDT,LNLG)
        SZFS=SIZE*(XVPR-XVPL)
        CALL    SET (0.,1.,0.,1.,0.,1.,0.,1.,1)
        CALL PCGETI ('QU - QUALITY FLAG',IQUA)
        CALL PCSETI ('QU - QUALITY FLAG',0)
        CALL PCSETI ('TE - TEXT EXTENT COMPUTATION FLAG',1)
c       CALL PLCHHQ (.5,.5,LABL,SZFS,360.,0.)
        CALL PCGETR ('DB - DISTANCE TO BOTTOM OF STRING',DBOS)
        CALL GQPLCI (IERR,ICIR)
        CALL GSPLCI (0)
        CALL PLCHHQ (.5*(XVPL+XVPR),YVPT+SZFS+DBOS+TMAG,
     1             LABL,SZFS,0.,0.)
        CALL GSPLCI (1)
        CALL PCSETI ('QU - QUALITY FLAG',IQUA)
        CALL    SET (XVPL,XVPR,YVPB,YVPT,XWDL,XWDR,YWDB,YWDT,LNLG)
c
c Done
c
      return
c
      end
c
c=================  LABLFT ======================
c
      SUBROUTINE LABLFT(LABL,SIZE)
c
      CHARACTER*(*) LABL
c
c  Put a label just above the top of the plot The SET call is re-done
c  to allow for the use of fractional coordinates, and the text extent
c  capabilities of the package PLOTCHAR are used to determine the label
c  position.
c
c  To put a slight margin from the top of the frame
      TMAG=0.05
c     
      CALL GETSET (XVPL,XVPR,YVPB,YVPT,XWDL,XWDR,YWDB,YWDT,LNLG) 
      SZFS=SIZE*(XVPR-XVPL)
      CALL    SET (0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL PCGETI ('QU - QUALITY FLAG',IQUA)
      CALL PCSETI ('QU - QUALITY FLAG',0)    
      CALL PCGETI ('TE - TEXT EXTENT COMPUTATION FLAG',1)
c     CALL PLCHHQ (.5,.5,LABL,SZFS,360.,0.)
      CALL PCGETR ('DB - DIATANCE TO BOTTOM OF STRING',DBOS)
c  This treatment does not work
c       CALL GSTXCI (0)
        CALL GQPLCI (IERR,ICIR)
      write(6,*) 'ierr= ',IERR
        CALL GSPLCI (0)
c     CALL PLCHHQ (.5*(XVPL+XVPR),YVPT+SZFS+DBOS+TMAG,
      CALL PLCHHQ (0.1,0.55,
     1             LABL,SZFS,90.,0.)
c       CALL GSTXCI (1)
        CALL GSPLCI (1)
      CALL PCSETI ('QU - QUALITY FLAG',IQUA) 
      CALL    SET (XVPL,XVPR,YVPB,YVPT,XWDL,XWDR,YWDB,YWDT,LNLG) 
c
c Done
c
      return
c
      end
c
      subroutine dscale(grid,bum,xmin,xmax,ymin,ymax,
     1                 alnch,alofx,msxlb,msylb,message,isc)
      dimension grid(4),bum(2)
      character*(*) msxlb,msylb,message
      dimension xd(3),yd(3)                                             
      data xd/11.,12.,15./, yd/1.0,2.0,3.0/                             
       if(isc.eq.1) ysc=10.0
       if(isc.eq.2) ysc=1.0
       if(isc.eq.3) ysc=0.1
      call agsetp('GRID WINDOW.',grid,4)                               
      call agsetf('X/MIN.',xmin)                                        
      call agsetf('X/MAX.',xmax)                                        
      call agseti('X/NICE.',0)                                          
      call agseti('BOTTOM/MAJOR/TYPE.',1)
      call agsetf('BOTTOM/MAJOR/BASE.',500.)
      call agseti('BOTTOM/MINOR/SPACING.',4)
      call agsetf('Y/MIN.',ymin)
      call agsetf('Y/MAX.',ymax)
      call agseti('Y/NICE.',0)
      call agseti('LEFT/MAJOR/TYPE.',1)
c     call agsetf('LEFT/MAJOR/BASE.',20.)
      call agsetf('LEFT/MAJOR/BASE.',ysc)
      call agseti('LEFT/MINOR/SPACING.',4)

c     if(xmax.gt.200.) then
c     call agseti('LEFT/MAJOR/TYPE.',1)                                
c     call agsetf('LEFT/MAJOR/BASE.',20.)                             
c     call agseti('LEFT/MINOR/SPACING.',1)                             
c     else
      call agseti('Y/NICE.',0)                                          
c     endif
c          
c      call agseti('BOTTOM/MAJOR/TYPE.',1)
c      call agsetf('BOTTOM/MAJOR/BASE.',20.)
c      call agseti('BOTTOM/MINOR/SPACING.',1)
c LEFT LABEL
      call agsetc('LABEL/NAME.','L')
      call agseti('LINE/NUMBER.',100)
      call agsetc('LINE/TEXT.',' ')
      call agsetf('LINE/CH.',0.040)
      call agsetf('LAB/OF/X.',-0.01)
c     call agsetc('LINE/TEXT.',msylb)
c     call agsetf('LINE/CH.',alnch)
c     call agsetf('LAB/OF/X.',alofx)
      call agsetp('LEFT/WIDTH.',bum,2)
c TOP LABEL
      call agsetc('LABEL/NAME.','T')
      call agseti('LINE/NUMBER.',100)
      call agsetc('LINE/TEXT.',message)
c     call agsetf('LINE/CH.',0.040)
c     call agsetf('LAB/OF/Y.',0.08)
c BOTTOM LABEL
      call agsetc('LABEL/NAME.','B')
      call agseti('LINE/NUMBER.',-100)
      call agsetf('LAB/OF/Y.',-0.04)
      call agsetp('BOTTOM/WIDTH.',bum,2)
c     call agsetc('LINE/TEXT.',msxlb)
      call agsetc('LINE/TEXT.',' ')
      call agsetf('LINE/CH.',0.040)
c
      call agstup(xd(1),1,0,3,1,yd(1),1,0,3,1)
      call agback
      return
      end
      subroutine aminmax(kxr,kyr,deni,dene,
     1             gmmaxi,gmmini,gmmaxe,gmmine)
      dimension deni(kxr,kyr),dene(kxr,kyr)
      gmmaxi=-10.
      gmmini=1000.
      gmmaxe=-10.
      gmmine=1000.
      do 767 m=1,kxr
      do 767 n=1,kyr
      if(deni(m,n).gt.gmmaxi) gmmaxi=deni(m,n)
      if(deni(m,n).lt.gmmini) gmmini=deni(m,n)
      if(dene(m,n).gt.gmmaxe) gmmaxe=dene(m,n)
      if(dene(m,n).lt.gmmine) gmmine=dene(m,n)
  767 continue
      return
      end
      subroutine arwmax(kxr,kyr,fvex,fvey,fvix,fviy,
     1            flo1,hi1,flo2,hi2)
      dimension fvex(kxr,kyr),fvey(kxr,kyr),
     1          fvix(kxr,kyr),fviy(kxr,kyr)
      flo1=1.0e10
      hi1=-1.0
      flo2=1.0e10
      hi2=-1.0
c
      do 31 i = 1, kxr
      do 31 j = 1, kyr
      exy=sqrt(fvex(i,j)**2 +fvey(i,j)**2)
      bxy=sqrt(fvix(i,j)**2 +fviy(i,j)**2)
      if(exy.le.flo1) flo1=exy
      if(exy.ge.hi1)  hi1=exy
      if(bxy.le.flo2) flo2=bxy
      if(bxy.ge.hi2)  hi2=bxy
   31 continue
c
      return
      end
      subroutine eplch(xpob,ypob,xpol,ypol,size,ichs)
C
C ---------------------------------------------------------------------
C
C Define a couple of temporary variables of type CHARACTER.
C
      CHARACTER*16 CTMP
C
C Do a call to SET which allows us to use fractional coordinates.
C
      CALL SET (0.,1.,0.,1.,0.,1.,0.,1.,1)
c
C Change to the appropriate character set.
C
        CALL PCSETI ('CD',0)
C     
C     
C Force constant spacing of the characters used for the column and row
C labels, so that they will line up better with each other.
C     
          CALL PCSETR ('CS',1.25)
C     
      if(ichs.ne.3) then
              CTMP=':IRU:X/:IGU:D:'
      else
              CTMP=':IRU:Y/:IGU:D:'
      endif
              CALL PLCHHQ (xpob,ypob,CTMP,size,0.,0.)
c
      if(ichs.eq.1) then 
c             CTMP=':IRU:Y/:IGU:D:'
              CTMP=':IRU:F:'
      else
              CTMP=':IRU:Z/:IGU:D:'
      endif
              CALL PLCHHQ (xpol,ypol,CTMP,size,90.,0.)
C
C Done.
C
      return
C
C ---------------------------------------------------------------------
C
      END
c    random number generator
c ** this is ran2 routine from Numerical Recepies (name changed to ran1) **
c     real*8 FUNCTION ran1(idum)
      real*4 FUNCTION ran1(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
c      REAL ran2,AM,EPS,RNMX
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     &            IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     &      IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
       enddo 
       iy=iv(1)
      endif
      
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2 
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
c      ran2=min(AM*iy,RNMX)
      ran1=min(AM*iy,RNMX)
      return
      END

