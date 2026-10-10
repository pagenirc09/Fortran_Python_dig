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
c     parameter (Nproc=3136)

      parameter (Nproc=4356)
c     parameter (Npx=16,Npy=2,Npz=2)
c     parameter (mxs=645,mys=761,mzs=761)
      parameter (mxs=1285,mys=797,mzs=797)

      parameter (nptl=6000000)
      parameter (mb=nptl,mj=9000000)

c     parameter (mg = mb*Nproc,mjg = mj*Nproc)
c     integer*8,parameter:: mg=mb*Nproc
c     integer*8,parameter:: mjg=mj*Nproc

      integer*8 iong, lecg
      integer*4 ionjg, lecjg

      parameter (imax=32767)
c  attention 90 dgree rotation
      parameter (mxx=200,myy=200,mzz=200)
c attention
c     parameter (mxr1=52,   myr1=27,   mzr1=27)
c     parameter (mxr=53,  myr=28,   mzr=28)
c     parameter (mxq=28,   myq=53,   mzq=53)
      parameter (mxr1=40,   myr1=40,   mzr1=40)
      parameter (mxr=40,  myr=40,   mzr=40)
      parameter (mxq=40,   myq=40,   mzq=40)
c     parameter (mxr1=50,   myr1=12,   mzr1=12)
c     parameter (mxr=50,  myr=12,   mzr=12)
c     parameter (mxq=50,  myq=12,   mzq=12)
c     real*8 me,mi,mi2
      dimension rhi(mxs,mys,mzs),rhe(mxs,mys,mzs)
      dimension eki(mxs,mys,mzs),eke(mxs,mys,mzs)
      dimension fex(mxs,mys,mzs),fey(mxs,mys,mzs),fez(mxs,mys,mzs)
      dimension fix(mxs,mys,mzs),fiy(mxs,mys,mzs),fiz(mxs,mys,mzs)
      dimension rhij(mxs,mys,mzs),rhej(mxs,mys,mzs)
      dimension ekij(mxs,mys,mzs),ekej(mxs,mys,mzs)
      dimension fexj(mxs,mys,mzs),feyj(mxs,mys,mzs),fezj(mxs,mys,mzs)
      dimension fixj(mxs,mys,mzs),fiyj(mxs,mys,mzs),fizj(mxs,mys,mzs)

c     dimension ixx(mb),iyy(mb),izz(mb)
c  ambient ions and electrons positions and velocities arrays
c     real*4 xi(mb),yi(mb),zi(mb),ui(mb),vi(mb),wi(mb)
c     real*4 xe(mb),ye(mb),ze(mb),ue(mb),ve(mb),we(mb)

c  jet ions and electrons positions and velocities arrays
c     real*4 xij(mj),yij(mj),zij(mj),uij(mj),vij(mj),wij(mj)
c     real*4 xej(mj),yej(mj),zej(mj),uej(mj),vej(mj),wej(mj)

c  global field arrays
c     real bxx,byy,bzz,exx,eyy,ezz
      real*4 bxx(mxs,mys,mzs),byy(mxs,mys,mzs),bzz(mxs,mys,mzs)
      real*4 exx(mxs,mys,mzs),eyy(mxs,mys,mzs),ezz(mxs,mys,mzs)
c  global particle arrays
c     real xig,yig,zig,uig,vig,wig
c     real xeg,yeg,zeg,ueg,veg,weg
c     real xijg,yijg,zijg,uijg,vijg,wijg
c     real xejg,yejg,zejg,uejg,vejg,wejg
c     real*4 xig(mg),yig(mg),zig(mg),uig(mg),vig(mg),wig(mg)
c     real*4 xeg(mg),yeg(mg),zeg(mg),ueg(mg),veg(mg),weg(mg)
c     real*4 xijg(mjg),yijg(mjg),zijg(mjg)
c     real*4 uijg(mjg),vijg(mjg),wijg(mjg)
c     real*4 xejg(mjg),yejg(mjg),zejg(mjg)
c     real*4 uejg(mjg),vejg(mjg),wejg(mjg)

c
      dimension dense(myy,mzz),densi(myy,mzz),ajy(myy,mzz),
     1          btot(myy,mzz),etot(myy,mzz),vgam(myy,mzz),
     1          denqe(myy,mzz),denqi(myy,mzz),aey(myy,mzz)
      dimension fvlex(myr,mzr),fvley(myr,mzr),fvlix(myr,mzr),
     1   fvliy(myr,mzr),fvqex(myq,mzq),fvqey(myq,mzq),
     1   fvqix(myq,mzq),fvqiy(myq,mzq)
c     dimension argx(mxs),argy(mys)
c     dimension xeb(mbeam), yeb(mbeam),
c    1          zeb(mbeam), vzb(mbeam),
c    2          xeg0(mback), yeg0(mback),
c    3          zeg0(mback), vzg(mback)
c
      dimension rho(mxs,mys),rho1(mxs,mzs)
      dimension gridr(4),gridm(4),gride(4),bnum(2),bnu1(2) 
      dimension xd(3),yd(3)                                             
c Display size
c     parameter (mxr=70   , myr= 30   , mzr= 30   )
c     dimension ex(mxs,mys,mzs),ey(mxs,mys,mzs),ez(mxs,mys,mzs),
c    &bx(mxs,mys,mzs),by(mxs,mys,mzs),bz(mxs,mys,mzs)
      dimension sm(27),ms(27),table(64)
c
      dimension delzx(myy,myy), delzy(myy,mzz),
     &          demzx(myy,myy), demzy(myy,mzz),
     &          delyx(mxx,mzz), delyz(mxx,mzz),
     &          demyx(mxx,mzz), demyz(mxx,mzz),
     &          delxy(myy,mzz), delxz(myy,mzz),
     &          demxy(myy,mzz), demxz(myy,mzz)
c
      dimension elzx(myr1,mzr1), elzy(myr1,mzr1),
     &          emzx(myr1,mzr1), emzy(myr1,mzr1),
     &          elyx(mxr1,mzr1), elyz(mxr1,mzr1),
     &          emyx(mxr1,mzr1), emyz(mxr1,mzr1),
     &          elxy(myq,mzq), elxz(myq,mzq),
     &          emxy(myq,mzq), emxz(myq,mzq)
c
c     data gride/0.10,0.75,0.35,0.593/
c     data gridr/0.275,0.70,0.15,0.95/
c     data gridr/0.10,0.75,0.300,0.460/
      data gridr/0.10,0.75,0.10,0.75/
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

      character message*28                                              
      character message1*22
      character*1 ittmp(21)
      open(6,file='dwriteot',status='unknown',form='formatted')

      call opngks
      call gsclip (0)

c specify "nst" - which data need to be processed
c     open(unit=2, file='nstn', status='old')

c      read(2,*) nst
      
c      nst = 5
      
       c = 1.0

       do ix=10,10
c      do ix=7,7

c      do ix=8,8

c      do ix=6,6
c      do ix=9,1,-1
c      do ix=7,26,2

       nst = ix

       strpd = 'partd_'
       hyph='_'

       st01 = char(int(nst/100.)+48)
       st02 = char(int((nst-int(nst/100.)*100)/10.)+48)
       st03 = char(int(nst-int(nst/10.)*10)+48)
       st0 = st01//st02//st03
       st = hyph//st0

       write(6,*) 'before reading'  

c  read from file
        dchar1='bfield'
        file1= dchar1//st

        open(15,file=file1,form='unformatted')
c       open(15,file=bfield_007,form='unformatted')
        read(15) bxx,byy,bzz
        close(15)
       write(6,*) 'after reading 15'  
       write(6,*)'bxx(300,35,35),byy(300,35,35),bzz(300,35,35)=',
     1 bxx(300,35,35),byy(300,35,35),bzz(300,35,35)

        dchar1='efield'
        file1= dchar1//st

        open(16,file=file1,form='unformatted')
c       open(16,file=efield_007,form='unformatted')
        read(16) exx,eyy,ezz
        close(16)
       write(6,*)'exx(300,35,35),eyy(300,35,35),ezz(300,35,35)=',
     1 exx(300,35,35),eyy(300,35,35),ezz(300,35,35)

       write(6,*) 'after reading 16'  
        dchar1='ionj3D'
        file1= dchar1//st

        open(18,file=file1,form='unformatted')
        read(18) rhij,fixj,fiyj,fizj,ekij
        close(18)

       write(6,*) 'after reading 17'  
        dchar1='lecs3D'
        file1= dchar1//st

        open(19,file=file1,form='unformatted')
        read(19) rhe,fex,fey,fez,eke
        close(19)

       write(6,*) 'after reading 16'  
        dchar1='lecj3D'
        file1= dchar1//st

        open(20,file=file1,form='unformatted')
        read(20) rhej,fexj,feyj,fezj,ekej
        close(20)


c       message2='kenout'
c       open(21,file=message2,form='unformatted')
c       dchar1='kenout'
c       file1= dchar1//st
c       open(21,file=file1,form='unformatted')
c       read(21) nstep
c    &   ,c,DT,qi,qe,mi,me,qmi,qme,vithml,vethml,vijet,vejet,
c    &         vithmj,vethmj
c    &         vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c new Jacek
c    &          b0y,e0z
c       close(21)
       write(6,*) 'after reading 21,nstep=',nstep




      je=47
c     ke=48
c     ie=70
c     in=130
c     in=160
      in=180
      in=280
c     in=380
      in=480
      in=101
c     in=580
      jn=381
      kn=381
c x_start
      istart= 419
c y_start
      jstart= 280
c z_start
      kstart= 280
c  for nst=6
      if(ix.eq.1) then
c     in=500
      in=430
      endif
      if(ix.eq.2) then
      in=450
      endif
      if(ix.eq.3) then
      in=460
      endif
      if(ix.eq.4) then
      in=470
      endif
      if(ix.eq.5) then
      in=480
      endif
      if(ix.eq.6) then
      in=485
      endif
      if(ix.eq.7) then
      in=490
      endif
      if(ix.eq.8) then
      in=500
      endif
      if(ix.eq.9) then
      in=900
      endif
      if(ix.eq.10) then
      in=520
      endif
      imid=in
      jmid=jn
      kmid=kn
c
c     do 81 i=1,mxs
c     do 81 j=1,mys
c     rho(i,j)=ex(i,j,kn)-ex(i-1,j,kn)
c    &        +ey(i,j,kn)-ey(i,j-1,kn)
c    &        +ez(i,j,kn)-ez(i,j,kn-1)
c  81 continue
c
c
c     do 82 i=1,mxs
c     do 82 k=1,mzs
c     rho1(i,k)=ex(i,jn,k)-ex(i-1,jn,k)
c    &         +ey(i,jn,k)-ey(i,jn-1,k)
c    &         +ez(i,jn,k)-ez(i,jn,k-1)
c  82 continue
c
c attention  x - y -->    y - z  plane  z-->x  x --> y   y ---> z
c     do 91 i=1, myy
c     do 91 j=1, mzz
c     i2=i+5
c     i3=i+5
c     j2=j+5
c     j3=j+5
c     delzx(i,j)=0.5*(eyy(imid,i2,j2) +eyy(imid,i2-1,j2))
c     delzy(i,j)=0.5*(ezz(imid,i2,j2) +ezz(imid,i2,j2-1))
c     demzx(i,j)=0.5*(byy(imid,i2,j2) +byy(imid-1,i2,j2-1))
c     demzy(i,j)=0.5*(bzz(imid,i2,j2) +bzz(imid-1,i2-1,j2))
c  91 continue
c
c attention        x-z plane Ex,z, Bx,z
c     b0z=b0z*0.992
c     do 92 i=1, mxx
c     do 92 k=1, mzz
c     k2=k+6
c     k3=k+6
c     i2=i+6
c     i3=i+6
c     delyx(i,k)=0.5*(exx(i3,jmid,k2) +exx(i3-1,jmid,k2))
c     delyz(i,k)=0.5*(ezz(i2,jmid,k2) +ezz(i2,jmid,k2-1))
c     demyx(i,k)=0.5*(bxx(i2,jmid,k2) +bxx(i2,jmid-1,k2-1))
c     demyz(i,k)=0.5*(bzz(i3,jmid,k2) +bzz(i3-1,jmid-1,k2))
c    1    -b0z
c  92 continue
c
      do 93 j=1, myy
      do 93 k=1, mzz
      j2=j+jstart
      j3=j+jstart
      k2=k+kstart
      k3=k+kstart
      delxy(j,k)=0.5*(eyy(imid,j3,k2) +eyy(imid,j3-1,k2))
      delxz(j,k)=0.5*(ezz(imid,j2,k3) +ezz(imid,j2,k3-1))
      demxy(j,k)=0.25*(byy(imid,j2,k3) +byy(imid-1,j2,k3-1))
     1          +0.25*(byy(imid-1,j2,k3) +byy(imid,j2,k3-1))
      demxz(j,k)=0.25*(bzz(imid,j3,k2) +bzz(imid-1,j3-1,k2))
     1          +0.25*(bzz(imid-1,j3,k2) +bzz(imid,j3-1,k2))
c    1    -b0z
   93 continue
c
c Averaging for matching with density plot
c
c     do 387 ii=1,myr1
c         i=18*ii -12
c     do 387 jj=1,mzr1
c         j=18*jj -12
c     elzx(ii,jj)=0.25*(delzx(i,j)+delzx(i+1,j)+delzx(i,j+1)
c    1                                         +delzx(i+1,j+1))
c     elzy(ii,jj)=0.25*(delzy(i,j)+delzy(i+1,j)+delzy(i,j+1)
c    1                                         +delzy(i+1,j+1))
c     emzx(ii,jj)=0.25*(demzx(i,j)+demzx(i+1,j)+demzx(i,j+1)
c    1                                         +demzx(i+1,j+1))
c     emzy(ii,jj)=0.25*(demzy(i,j)+demzy(i+1,j)+demzy(i,j+1)
c    1                                         +demzy(i+1,j+1))
c387  continue

c     write(6,*) 'emzx, elzx= ',emzx(5,3) ,elzx(5,3)
c
      imin=0 
      do 587 ii=1,myq  
      i = 5*(ii-1) +2
      do 587 jj=1,mzq 
      j = 5*(jj-1) +2
c  attention
      elxy(ii,jj)=delxy(i,j)
c     elxy(ii,jj)=0.25*(delxy(i,j)+delxy(i+1,j)+delxy(i,j+1)
c    1                                         +delxy(i+1,j+1))
      elxz(ii,jj)=delxz(i,j)
c     elxz(ii,jj)=0.25*(delxz(i,j)+delxz(i+1,j)+delxz(i,j+1)
c    1                                         +delxz(i+1,j+1))
c for y - z
      emxy(ii,jj)=demxy(i,j)
c     emxy(ii,jj)=0.25*(demxy(i,j)+demxy(i+1,j)+demxy(i,j+1)
c    1                                         +demxy(i+1,j+1))
c for y - z
      emxz(ii,jj)=demxz(i,j)
c     emxz(ii,jj)=0.25*(demxz(i,j)+demxz(i+1,j)+demxz(i,j+1)
c    1                                         +demxz(i+1,j+1))
c for x - z (attention)
c     emxy(ii,jj)=demyx(i+imin+1,j)
c     emxz(ii,jj)=demyz(i+imin+1,j)
c for x - y (attention)
c     emxy(ii,jj)=demzx(i+imin+1,j)
c     emxz(ii,jj)=demzy(i+imin+1,j)
 587  continue
c
      write(6,*) 'emxz, elxz= ',emxz(5,3) ,elxz(5,3)

      amp=0.20
c     call arwmax(myr1,mzr1,elzx,elzy,emzx,emzy,
c    1            ezlo,ezhi,azlo,azhi)
c     call arwmax(mxr1,mzr1,elyx,elyz,emyx,emyz,
c    1            eylo,eyhi,aylo,ayhi)
      call arwmax(myq,mzq,elxy,elxz,emxy,emxz,
     1            exlo,exhi,axlo,axhi)
c     ezhi=ezhi*amp
c     azhi=azhi*amp
c     eyhi=sqrt(eyhi)*amp
c     ayhi=sqrt(ayhi)*amp
      exhi=exhi*amp
      axhi=axhi*amp
c     write(6,*)'ezhi,azhi,eyhi,ayhi,exhi,axhi=',
c    1  ezhi,azhi,eyhi,ayhi,exhi,axhi
c
c     do 377 i = 1, myr1
c     do 377 j = 1, mzr1
c     el=sqrt(elyx(i,j)**2+elyz(i,j)**2)
c     em=sqrt(emyx(i,j)**2+emyz(i,j)**2)
c     elyx(i,j)=elyx(i,j)/(el+eyhi)
c     elyz(i,j)=elyz(i,j)/(el+eyhi)
c     emyx(i,j)=emyx(i,j)/(em+ayhi)
c     emyz(i,j)=emyz(i,j)/(em+ayhi)
c     elyx(i,j)=elyx(i,j)/(abs(elyx(i,j))+eyhi)
c     elyz(i,j)=elyz(i,j)/(abs(elyz(i,j))+eyhi)
c     emyx(i,j)=emyx(i,j)/(abs(emyx(i,j))+ayhi)
c     emyz(i,j)=emyz(i,j)/(abs(emyz(i,j))+ayhi)
c 377 continue
      do 372 i = 1, myq
      do 372 j = 1, mzq
      el=sqrt(elxy(i,j)**2+elxz(i,j)**2)
      em=sqrt(emxy(i,j)**2+emxz(i,j)**2)
      elxy(i,j)=elxy(i,j)/(el+exhi)
      elxz(i,j)=elxz(i,j)/(el+exhi)
c attention
      emxy(i,j)=emxy(i,j)/(em+axhi)
      emxz(i,j)=emxz(i,j)/(em+axhi)
  372 continue
c
c
c     call arwmax(mxr1,myr1,elzx,elzy,emzx,emzy,
c    1            ezlo,ezhi,azlo,azhi)
c     call arwmax(mxr1,mzr1,elyx,elyz,emyx,emyz,
c    1            eylo,eyhi,aylo,ayhi)
      call arwmax(myq,mzq,elxy,elxz,emxy,emxz,
     1            exlo,exhi,axlo,axhi)
c
c
       write(6,*) 'after field calculation'  
c
      ymin=float(jn)-2.0
      ymax=float(jn)+2.0
      zmin=float(kn)-2.0
      zmax=float(kn)+2.0
c     xliml=150.
c     xlimr=290.
c     xliml=0.
c     xlimr=150.
c
c  attention the sustem is original size  
c
c attention 
      do 11 j=1,mzz
         jj=j+kstart
      do 11 i=1,myy
         ii=i+jstart
      densi(i,j)= rhij(in,ii,jj)
      dense(i,j)=rhe(in,ii,jj)+rhej(in,ii,jj)
      btot(i,j)=sqrt(
     1  (0.5*(bxx(in,ii,jj) +bxx(in,ii-1,jj-1)))**2
     1 +(0.5*(byy(in,ii,jj) +byy(in-1,ii,jj-1)))**2
     1 +(0.5*(bzz(in,ii,jj) +bzz(in-1,ii-1,jj)))**2)
      etot(i,j)=sqrt(
     1  (0.5*(exx(in,ii,jj) +exx(in-1,ii,jj)))**2
     1 +(0.5*(eyy(in,ii,jj) +eyy(in,ii-1,jj)))**2
     1 +(0.5*(ezz(in,ii,jj) +ezz(in,ii,jj-1)))**2)
c attention   ajy = Jx
      ajy(i,j)=0.5*(exx(in,ii,jj) +exx(in-1,ii,jj))
c     ajy(i,j)=0.5*(fix(in,ii,jj)-fex(in,ii,jj))
c    1         +0.25*(fix(in,ii,jj-1)-fex(in,ii,jj-1))
c    1         +0.25*(fix(in,ii,jj+1)-fex(in,ii,jj+1))
c    1        +0.5*(fixj(in,ii,jj)-fexj(in,ii,jj))
c    1         +0.25*(fixj(in,ii,jj-1)-fexj(in,ii,jj-1))
c    1         +0.25*(fixj(in,ii,jj+1)-fexj(in,ii,jj+1))
      vtemp= 1.0 -(fex(in,ii,jj)**2+fey(in,ii,jj)**2
     1      +fez(in,ii,jj)**2)/(c*rhe(in,ii,jj))**2
     1           -(fexj(in,ii,jj)**2+feyj(in,ii,jj)**2
     1      +fezj(in,ii,jj)**2)/(c*rhej(in,ii,jj))**2
      if(vtemp.le.0.0) vtemp=0.0001
      vgam(i,j)=1.0/sqrt(vtemp)
   11 continue

c      write(6,*)'dense(300,35),densi(300,35),ajy(300,35)=',
c    1  dense(300,35),densi(300,35),ajy(300,35)

      do 311 j=1,mzr
      jj=5*(j-1)+kstart +3 
      do 311 i=1,myr
      ii=5*(i-1)+jstart+3 
      fvlix(i,j)=0.5*fiy(in,ii,jj)
     1        +0.25*(fiy(in,ii,jj-1)+fiy(in,ii,jj+1))
     1          +0.5*fiyj(in,ii,jj)
     1        +0.25*(fiyj(in,ii,jj-1)+fiyj(in,ii,jj+1))
      fvliy(i,j)=0.5*fiz(in,ii,jj)
     1        +0.25*(fiz(in,ii,jj-1)+fiz(in,ii,jj+1))
     1          +0.5*fizj(in,ii,jj)
     1        +0.25*(fizj(in,ii,jj-1)+fizj(in,ii,jj+1))
      fvlex(i,j)=0.5*fey(in,ii,jj)
     1        +0.25*(fey(in,ii,jj-1)+fey(in,ii,jj+1))
     1          +0.5*feyj(in,ii,jj)
     1        +0.25*(feyj(in,ii,jj-1)+feyj(in,ii,jj+1))
      fvley(i,j)=0.5*fez(in,ii,jj)
     1        +0.25*(fez(in,ii,jj-1)+fez(in,ii,jj+1))
     1          +0.5*fezj(in,ii,jj)
     1        +0.25*(fezj(in,ii,jj-1)+fezj(in,ii,jj+1))
  311 continue
c
      call aminmax(myy,mzz,densi,dense,
     1             gmmaxi,gmmini,gmmaxe,gmmine)
      call aminmax(myy,mzz,btot,etot,    
     1             gmmaxbt,gmminbt,gmmaxet,gmminet)
      call aminmax(myy,mzz,ajy,vgam,
     1             gmmaxjy,gmminjy,gmmaxvg,gmminvg)

      write(6,*)'gmmine,gmmaxe=',gmmine,gmmaxe
c
c  attention setting the minimum and maximum with logarismic scale 
c     dense(1,1)=7.0
c     dense(mxx,mzz)=160.

c     do j=1,mzz
c     do i=1,mxx
c     densi(i,j)=alog(densi(i,j))
c     dense(i,j)=alog(dense(i,j))
c     enddo
c     enddo
c     call aminmax(myy,mzz,densi,dense,
c    1             gmmaxi,gmmini,gmmaxe,gmmine)
c
c
 7139 continue
c
c
c     go to 3927
      amp=0.01
      call arwmax(myr,mzr,fvlex,fvley,fvlix,fvliy,
     1            flo1,hi1,flo2,hi2)
      hi1=hi1*amp
      hi2=hi2*amp
c
      do 37 i = 1, myr
      do 37 j = 1, mzr
      fve=sqrt(fvlex(i,j)**2+fvley(i,j)**2)
      fvi=sqrt(fvlix(i,j)**2+fvliy(i,j)**2)
      fvlex(i,j)=fvlex(i,j)/(fve+hi1)
      fvley(i,j)=fvley(i,j)/(fve+hi1)
      fvlix(i,j)=fvlix(i,j)/(fvi+hi2)
      fvliy(i,j)=fvliy(i,j)/(fvi+hi2)
   37 continue
c
 3927 continue
c
      call arwmax(myr,mzr,fvlex,fvley,fvlix,fvliy,
     1            flo1,hi1,flo2,hi2)

c attention   x  --> y    y --> z
      do 215 k=1,mzq
      kk=4*k -2
      do 215 i=1,myq
      ii=4*i -2
c attention for notation x, y
      fvqix(i,k)=fiz(in,ii,kk)
      fvqiy(i,k)=fix(in,ii,kk)
      fvqex(i,k)=fez(in,ii,kk)
      fvqey(i,k)=fex(in,ii,kk)
  215 continue
c
      do 205 j=1,mzz
         ii=j+7
      do 205 i=1,myy
         jj=i+7
      denqi(i,j)=eki(in,jj,ii)
      denqe(i,j)=eke(in,jj,ii)
c attention for notation x, y
      aey(i,j)  =0.5*(eyy(in,jj,ii) +eyy(in,jj-1,ii))
  205 continue
c
      call aminmax(myy,mzz,denqi,denqe,
     1             gmmaxki,gmminki,gmmaxke,gmminke)
      call aminmax(myy,mzz,aey,aey,   
     1             gmmaxey,gmminey,gmmaxey,gmminey)
      call arwmax(mxq,mzq,fvqex,fvqey,fvqix,fvqiy,
     1            flo1,hi1,flo2,hi2)

      t=nstep-1
c
      nset=-1
      length=0
      ispv=0
      spv=1.0
c attention
      xmin=istart +1
c     xmax=152.0
c     xmax=148.0
      xmax=istart+201
      ymin=jstart +1
c     ymax=77.0
c     ymax=73.0
      ymax=jstart+201
      zmin=kstart +1
c     zmax=77.0 
      zmax=kstart+201
c
c     if(gmmaxey.gt.abs(gmminey)) then
c     gmminey=-gmmaxey
c     aey(1,1)=-gmmaxey
c     else
c     gmmaxey=-gmminey
c     aey(1,1)=-gmminey
c     endif
      if(gmmaxjy.gt.abs(gmminjy)) then
      gmminjy=-gmmaxjy
      ajy(1,1)=-gmmaxjy
      else
      gmmaxjy=-gmminjy
      ajy(1,1)=-gmminjy
      endif


      ICHOS=315


c     ICHOS=114
c
      if(ichos.eq.110)
     .write(message,110) t
      if(ichos.eq.111)
     .write(message,111) t
      if(ichos.eq.112)
     .write(message,112) t
      if(ichos.eq.113)
     .write(message,113) t
      if(ichos.eq.114)
     .write(message,114) t
      if(ichos.eq.115)
     .write(message,115) t
      if(ichos.eq.116)
     .write(message,116) t
      if(ichos.eq.117)
     .write(message,117) t
      if(ichos.eq.118)
     .write(message,118) t
      if(ichos.eq.119)
     .write(message,119) t
      if(ichos.eq.315)
     .write(message,315) in


c     write(message,113) t                                               

        CALL DFCLRS
        CALL GSTXCI (1)
        CALL GSPLCI (1)

      if(ichos.eq.110.or.ichos.eq.315)
     .call colorcn(ajy,myy,mzz,gmminjy,gmmaxjy,gridr,message,2)

      if(ichos.eq.111)
     .call colorcn(dense,myy,mzz,gmmine,gmmaxe,gridr,message,1)

      if(ichos.eq.112)
     .call colorcn(dense,myy,mzz,gmmine,gmmaxe,gridr,message,1)

      if(ichos.eq.113)
     .call colorcn(densi,myy,mzz,gmmini,gmmaxi,gridr,message,1)

      if(ichos.eq.114)
     .call colorcn(btot,myy,mzz,gmminbt,gmmaxbt,gridr,message,1)

      if(ichos.eq.115)
     .call colorcn(etot,myy,mzz,gmminet,gmmaxet,gridr,message,1)

      if(ichos.eq.116)
     .call colorcn(vgam,myy,mzz,gmminvg,gmmaxvg,gridr,message,1)

      if(ichos.eq.117)
     .call colorcn(aey,myy,mzz,gmminey,gmmaxey,gridr,message,1)

      if(ichos.eq.118)
     .call colorcn(denqe,myy,mzz,gmminke,gmmaxke,gridr,message,1)

      if(ichos.eq.119)
     .call colorcn(denqi,myy,mzz,gmminki,gmmaxki,gridr,message,1)


       xvpl=gridr(1)+0.01
       xvpr=gridr(2)-0.01
       yvpb=gridr(3)+0.01
       yvpt=gridr(4)-0.01

       CALL BFCLRS
c red
c       CALL GSTXCI (1)
c       CALL GSPLCI (1)
c black arrows
        CALL GSTXCI (0)
        CALL GSPLCI (0)

      call set(xvpl,xvpr,yvpb,yvpt,0.0,645.0,0.0,85.0,1)

      IOFFM = 1

      if(ichos.eq.111.or.ichos.eq.315)
     .call velvct(emxy,myr1,emxz,myr1,myr1,mzr1,axlo,axhi,nset,length,
     &            ispv,spv)

      if(ichos.eq.115.or.ichos.eq.117.or.ichos.eq.316)
     .call velvct(elxy,myr1,elxz,myr1,myr1,mzr1,exlo,exhi,nset,length,
     &            ispv,spv)

      if(ichos.eq.110.or.ichos.eq.112.or.ichos.eq.116.or.ichos.eq.118)
     .call velvct(fvlex,myr1,fvley,myr1,myr1,mzr1,flo1,hi1,nset,length,
     &            ispv,spv)

      if(ichos.eq.115.or.ichos.eq.119)
     .call velvct(fvlix,myr1,fvliy,myr1,myr1,mzr1,flo2,hi2,nset,length,
     &            ispv,spv)


       CALL DFCLRS
        CALL GSTXCI (0)
        CALL GSPLCI (0)

      call dscale(gridr,bnu1,xmin,xmax,zmin,zmax, 
     1            0.04,-0.03,'Y$','Z$',message)

      xpob=0.425
      ypob=0.025
c     xpol=0.205
c     xpol=0.105
      xpol=0.020
      ypol=0.425

      size=0.018
c  for x - y
c     ichs=1
c  for x - z
c     ichs=2
c  for y - z
      ichs=3
      call eplch(xpob,ypob,xpol,ypol,size,ichs)
       CALL DFCLRS
        CALL GSTXCI (1)
        CALL GSPLCI (1)
      call frame
 8173 continue
      enddo

c 110 format('JY      (ELE)  T=',f8.1)
  110 format('JX      (EFLU) T=',f8.1)
  111 format('ELE DEN (MAG)  T=',f8.1)
  112 format('ELE DEN (EFLU) T=',f8.1)
c 113 format('ION DEN (IFLU) T=',f8.1)
  113 format('IODE JT (ELE)  T=',f8.1)
  114 format('BTOT    (MAG)  T=',f8.1)
  115 format('ETOT    (ELE)  T=',f8.1)
  116 format('GAMMAE  (EFLU) T=',f8.1)
  117 format('EY      (ELE)  T=',f8.1)
  118 format('EKIN    (EFLU) T=',f8.1)
  119 format('IKIN    (IFLU) T=',f8.1)
  214 format('CUR DEN (MAG)  T=',f8.1)
c 315 format('CUR DEN (MAG)  T=',f8.1)
  315 format('Ex      (MAG)  X=',i8)
  275 format('CUR DEN (FLUX) T=',f8.1)
  216 format('CUR DEN (CUR)  T=',f8.1)
  223 format('ELE ENE (MAG)  T=',f8.1)
  224 format('ION ENE (MAG)  T=',f8.1)
  323 format('ELE ENE (FLUX) T=',f8.1)
  324 format('ION ENE (FLUX) T=',f8.1)
   25 format('X-Y SPACE (ION) T=',F8.0,'$')                            
   26 format('X-Y SPACE (ELE) T=',F8.0,'$')                            
   27 format('Z-VZ  SPACE     T=',F8.0,'$')                            
c                                                                       
c     go to 4
      call clsgks
      stop
      end


      subroutine Vol_Weighting(mxs,mys,mzs,lecg,rhe,eke,ueg,veg,weg,
     &                         xeg,yeg,zeg,fex,fey,fez)
      implicit none

      integer mxs,mys,mzs,mg,lecg
      real*4 xeg(lecg),yeg(lecg),zeg(lecg)
      real*4 ueg(lecg),veg(lecg),weg(lecg)
      real rhe(mxs,mys,mzs),eke(mxs,mys,mzs)
      real fex(mxs,mys,mzs),fey(mxs,mys,mzs),fez(mxs,mys,mzs)
      real*4 dx,dy,dz,cx,cy,cz
      real*4 sl,sr,sn,s     
      integer i,j,k,l,m,n,n0

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
        parameter (nspc=366000,nwsp=5882000)
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
        CALL LBLBAR (1,.79,.89,.12,.72,14,.333,1.,LIND,0,LLBS,15,1)
c       else
c       CALL LBLBAR (1,.79,.89,.20,.80,14,.333,1.,LIND,0,LLBS,15,1)
c       CALL LBLBAR (0,.10,.95,.23,.30,14,1.,0.5,LIND,0,LLBS,15,1)
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
     1                 alnch,alofx,msxlb,msylb,message)
      dimension grid(4),bum(2)
      character*(*) msxlb,msylb,message
      dimension xd(3),yd(3)                                             
      data xd/11.,12.,15./, yd/1.0,2.0,3.0/                             
      call agsetp('GRID WINDOW.',grid,4)                               
      call agsetf('X/MIN.',xmin)                                        
      call agsetf('X/MAX.',xmax)                                        
      call agseti('X/NICE.',0)                                          
      call agseti('BOTTOM/MAJOR/TYPE.',1)
      call agsetf('BOTTOM/MAJOR/BASE.',100.)
      call agseti('BOTTOM/MINOR/SPACING.',4)
      call agsetf('Y/MIN.',ymin)
      call agsetf('Y/MAX.',ymax)
      call agseti('Y/NICE.',0)
      call agseti('LEFT/MAJOR/TYPE.',1)
      call agsetf('LEFT/MAJOR/BASE.',100.)
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
c     call agsetc('LINE/TEXT.',message)                                 
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
              CTMP=':IRU:Y/:IGU:D:'
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
