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
      parameter (Nproc=3136)
c     parameter (Npx=1,Npy=9,Npz=9)
c     parameter (mxs=2005,mys=1005,mzs=1005)
      parameter (mxs=1285,mys=789,mzs=789)

c     parameter (nptl=4600000)
c     parameter (mb=nptl,mj=3500000)

      parameter (nptl=6000000)
      parameter (mb=nptl,mj=9000000)

c     integer*8,parameter:: mg=mb*Nproc
c     integer*8,parameter:: mjg=mj*Nproc
      integer*8, parameter :: mg = int(mb,8) * int(Nproc,8)
      integer*8, parameter :: mjg = int(mj,8) * int(Nproc,8)


      integer*4 iong, lecg
      integer*4 ionjg, lecjg, n0


      parameter (imax=32767)
c  attention 90 dgree rotation
c     parameter (mxx=mxs-12,myy=mys-12,mzz=mzs-12)
c New ken
      parameter (mxx=mxs-25,myy=mys-15,mzz=mzs-15)
c     parameter (mxx=mxs-5,myy=mys-5,mzz=mzs-5)
c     parameter (mxk=630,myk=387,mzk=387)
      parameter (mxk=510,myk=125,mzk=125)
c     parameter (mxx=mxs-95,myy=mys-12,mzz=mzs-12)
c attention
c     parameter (mxr1=52,   myr1=27,   mzr1=27)
c     parameter (mxr=53,  myr=28,   mzr=28)
c     parameter (mxq=28,   myq=53,   mzq=53)
c     parameter (mxr1=50,   myr1=25,   mzr1=25)
c     parameter (mxr=50,  myr=25,   mzr=25)
c     parameter (mxq=50,   myq=25,   mzq=25)
c     parameter (mxr1=41,   myr1=25,   mzr1=25)
c     parameter (mxr=41,  myr=25,   mzr=25)
c     parameter (mxq=41,  myq=25,   mzq=25)
      parameter (mxr1=41,   myr1=12,   mzr1=12)
      parameter (mxr=40,  myr=12,   mzr=12)
      parameter (mxq=40,  myq=12,   mzq=12)
c     real*8 me,mi,mi2
      real*4 :: rhi(mxs,mys,mzs),rhe(mxs,mys,mzs)
      real*4 :: eki(mxs,mys,mzs),eke(mxs,mys,mzs)
      real*4 :: fex(mxs,mys,mzs),fey(mxs,mys,mzs),fez(mxs,mys,mzs)
      real*4 :: fix(mxs,mys,mzs),fiy(mxs,mys,mzs),fiz(mxs,mys,mzs)
      real*4 :: rhij(mxs,mys,mzs),rhej(mxs,mys,mzs)
      real*4 :: ekij(mxs,mys,mzs),ekej(mxs,mys,mzs)
      real*4 :: fexj(mxs,mys,mzs),feyj(mxs,mys,mzs),fezj(mxs,mys,mzs)
      real*4 :: fixj(mxs,mys,mzs),fiyj(mxs,mys,mzs),fizj(mxs,mys,mzs)

c     dimension ixx(mb),iyy(mb),izz(mb)
c  ambient ions and electrons positions and velocities arrays
c     real*4 xi(mb),yi(mb),zi(mb),ui(mb),vi(mb),wi(mb)
c     real*4 xe(mb),ye(mb),ze(mb),ue(mb),ve(mb),we(mb)

c  jet ions and electrons positions and velocities arrays
c     real*4 xij(mj),yij(mj),zij(mj),uij(mj),vij(mj),wij(mj)
c     real*4 xej(mj),yej(mj),zej(mj),uej(mj),vej(mj),wej(mj)

c  global field arrays
c     real bxx,byy,bzz,exx,eyy,ezz
      real*4 :: bxx(mxs,mys,mzs),byy(mxs,mys,mzs),bzz(mxs,mys,mzs)
      real*4 :: exx(mxs,mys,mzs),eyy(mxs,mys,mzs),ezz(mxs,mys,mzs)
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

      real*4 
     &   c4,DT4,qi4,qe4,qmi4,qme4,vithml4,vethml4,vijet4,vejet4,
     &   vithmj4,vethmj4,refli4,refle4,rselect4,xj04,yj04,zj04,
     &   b0x4,b0y4,e0z4
       real*4 gam
c
      dimension dense(mxk,mzk),densi(mxk,mzk),ajy(mxk,mzk),
     1          btot(mxk,mzk),etot(mxk,mzk),vgam(mxk,mzk),
     1          denqe(mxk,mzk),denqi(mxk,mzk),aey(mxk,mzk)
      dimension fvlex(mxr,mzr),fvley(mxr,mzr),fvlix(mxr,mzr),
     1   fvliy(mxr,mzr),fvqex(mxq,mzq),fvqey(mxq,mzq),
     1   fvqix(mxq,mzq),fvqiy(mxq,mzq)
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
c
c     data gride/0.10,0.75,0.35,0.593/
c     data gridr/0.275,0.70,0.15,0.95/
c     data gridr/0.10,0.75,0.300,0.460/
c     data gridr/0.10,0.95,0.380,0.871/
c     data gride/0.10,0.75,0.25,0.6191/
      data gridr/0.10,0.95,0.380,0.588/
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
      
c      do ix=1,9
c      do ix=10,10
       do ix=9,9
c      do ix=8,8
       nst=ix

c      nst = 5
      
       c = 1.0

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
c      write(6,*) 'after reading 15'  
c      write(6,*)'bxx(300,35,35),byy(300,35,35),bzz(300,35,35)=',
c    1 bxx(300,35,35),byy(300,35,35),bzz(300,35,35)

        dchar1='efield'
        file1= dchar1//st

        open(16,file=file1,form='unformatted')
c       open(16,file=efield_007,form='unformatted')
        read(16) exx,eyy,ezz
        close(16)
c      write(6,*)'exx(300,35,35),eyy(300,35,35),ezz(300,35,35)=',
c    1 exx(300,35,35),eyy(300,35,35),ezz(300,35,35)

       write(6,*) 'after reading 16'  
        dchar1='ions3D'
        file1= dchar1//st

cheung
        open(17,file=file1,form='unformatted')
        read(17) rhi,fix,fiy,fiz,eki
        close(17)

        dchar1='ionj3D'
        file1= dchar1//st

        open(18,file=file1,form='unformatted')
        read(18) rhij,fixj,fiyj,fizj,ekij
        close(18)

       dchar1='lecs3D'
        file1= dchar1//st

        open(19,file=file1,form='unformatted')
        read(19) rhe,fex,fey,fez,eke
        close(19)

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
c    &  ,c,DT,qi,qe,qmi,qme,vithml,vethml,vijet,vejet,
c    & vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,
c    &  vithmj,vethmj,refli,refle,rselect,xj0,yj0,zj0,b0x,
c    &  vithmj,vethmj
c new Jacek
c    &  b0y,e0z
c       close(21)
c      write(6,*) 'after reading 21,nstep=',nstep
c      print *, 'after reading 21,nstep=',nstep



c     write(6,*) ions, lecs,maxptl,maxhlf
c     write(6,*) 'ionm,ionj,lecm,lecj= ',ionm,ionj,lecm,lecj
c      ionm = 42  ! added by Johnny
c      ionj = 42  ! added by Johnny
c      lecm = 42  ! added by Johnny
c      lecj = 42  ! added by Johnny
c      print *, 'ionm,ionj,lecm,lecj= ',ionm,ionj,lecm,lecj
      je=47
c     ke=48
c     ie=70
      in=1000
      jn=381
      kn=381
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
c  x - y
c     do 91 i=1, mxx
c     do 91 j=1, myy
c     i2=i+4
c     i3=i+4
c     j2=j+4
c     j3=j+4
c     delzx(i,j)=0.5*(ex(i3,j2,kmid) +ex(i3-1,j2,kmid))
c     delzy(i,j)=0.5*(ey(i2,j3,kmid) +ez(i2,j3-1,kmid))
c     demzx(i,j)=0.5*(bx(i2,j3,kmid) +bx(i2,j3-1,kmid-1))
c     demzy(i,j)=0.5*(by(i3,j2,kmid) +by(i3-1,j2,kmid-1))
c  91 continue
c
c attention        x-z plane Ex,z, Bx,z
c     b0z = 42    ! added by Johnny
c     b0z=b0z*0.992
c
c   shifting the location
c
c      ishift = vejet*DT*nstep -350
       ishift = 0
c      ishift = 1800 -8
c     if(ishift.lt.0) ishift = 0
c     if((ishift+408).ge.4000) ishift = 4000-408
      write(6,*)'DT,vejet,nstep,ishift=',DT,vejet,nstep,ishift

      do 92 i=1, mxx
      do 92 k=1, mzz
      k2=k+7
      k3=k+7
      i2=i+12+ishift
      i3=i+12+ishift
      delyx(i,k)=0.5*(exx(i3,jmid,k2) +exx(i3-1,jmid,k2))
      delyz(i,k)=0.5*(ezz(i2,jmid,k2) +ezz(i2,jmid,k2-1))
      demyx(i,k)=0.25*(bxx(i2,jmid,k2) +bxx(i2,jmid-1,k2-1))
     1          +0.25*(bxx(i2,jmid-1,k2) +bxx(i2,jmid,k2-1))
      demyz(i,k)=0.25*(bzz(i3,jmid,k2) +bzz(i3-1,jmid-1,k2))
     1          +0.25*(bzz(i3-1,jmid,k2) +bzz(i3,jmid-1,k2))
c    1    -b0z
   92 continue
c
c     do 93 j=1, myy
c     do 93 k=1, mzz
c     j2=j+2
c     j3=j+2
c     k2=k+2
c     k3=k+2
c     delxy(j,k)=0.5*(ey(imid,j3,k2) +ey(imid,j3-1,k2))
c     delxz(j,k)=0.5*(ez(imid,j2,k3) +ez(imid,j2,k3-1))
c     demxy(j,k)=0.5*(by(imid,j2,k3) +by(imid-1,j2,k3-1))
c     demxz(j,k)=0.5*(bz(imid,j3,k2) +bz(imid-1,j3-1,k2))
c    1    -b0z
c  93 continue
c
c Averaging for matching with density plot
c
      do 387 ii=1,mxr1
c         i=31*ii -20 
          i=24*ii +96
      do 387 jj=1,mzr1
c         j=31*jj -15
          j=24*jj +229
c     elzx(ii,jj)=delzx(i+1,j+1)
c     elzy(ii,jj)=delzy(i+1,j+1)
c     elzy(ii,jj)=0.25*(delzy(i,j)+delzy(i+1,j)+delzy(i,j+1)
c    1                                         +delzy(i+1,j+1))
      elyx(ii,jj)=delyx(i+1,j+1)
c     elyx(ii,jj)=0.25*(delyx(i,j)+delyx(i+1,j)+delyx(i,j+1)
c    1                                         +delyx(i+1,j+1))
      elyz(ii,jj)=delyz(i+1,j+1)
c     elyz(ii,jj)=0.25*(delyz(i,j)+delyz(i+1,j)+delyz(i,j+1)
c    1                                         +delyz(i+1,j+1))
c     emzx(ii,jj)=demzx(i+1,j+1)
c     emzx(ii,jj)=0.25*(demzx(i,j)+demzx(i+1,j)+demzx(i,j+1)
c    1                                         +demzx(i+1,j+1))
c     emzy(ii,jj)=demzy(i+1,j+1)
c     emzy(ii,jj)=0.25*(demzy(i,j)+demzy(i+1,j)+demzy(i,j+1)
c    1                                         +demzy(i+1,j+1))
      emyx(ii,jj)=demyx(i+1,j+1)
c     emyx(ii,jj)=0.25*(demyx(i,j)+demyx(i+1,j)+demyx(i,j+1)
c    1                                         +demyx(i+1,j+1))
      emyz(ii,jj)=demyz(i+1,j+1)
c     emyz(ii,jj)=0.25*(demyz(i,j)+demyz(i+1,j)+demyz(i,j+1)
c    1                                         +demyz(i+1,j+1))
 387  continue

      write(6,*) 'emyx, emyz= ',emyx(5,3) ,emyz(5,3)
      print *, 'emyx, emyz= ',emyx(5,3) ,emyz(5,3)
c
      imin=0 
c     do 587 i=4,82,3 
c     ii=(i-1)/3
c     do 587 j=4,157,3
c     jj=(j-1)/3
c  attention
c     elxy(ii,jj)=delyx(i+1,j+1)
c     elxy(ii,jj)=0.25*(delxy(i,j)+delxy(i+1,j)+delxy(i,j+1)
c    1                                         +delxy(i+1,j+1))
c     elxz(ii,jj)=delyz(i+1,j+1)
c     elxz(ii,jj)=0.25*(delxz(i,j)+delxz(i+1,j)+delxz(i,j+1)
c    1                                         +delxz(i+1,j+1))
c for y - z
c     emxy(ii,jj)=demxy(i+1,j+1)
c     emxy(ii,jj)=0.25*(demxy(i,j)+demxy(i+1,j)+demxy(i,j+1)
c    1                                         +demxy(i+1,j+1))
c for y - z
c     emxz(ii,jj)=demxz(i+1,j+1)
c     emxz(ii,jj)=0.25*(demxz(i,j)+demxz(i+1,j)+demxz(i,j+1)
c    1                                         +demxz(i+1,j+1))
c for x - z (attention)
c     emxy(ii,jj)=demyx(i+imin+1,j)
c     emxz(ii,jj)=demyz(i+imin+1,j)
c for x - y (attention)
c     emxy(ii,jj)=demzx(i+imin+1,j)
c     emxz(ii,jj)=demzy(i+imin+1,j)
c587  continue
c
      amp=0.20
c     call arwmax(mxr1,mzr1,elzx,elzy,emzx,emzy,
c    1            ezlo,ezhi,azlo,azhi)
      call arwmax(mxr1,mzr1,elyx,elyz,emyx,emyz,
     1            eylo,eyhi,aylo,ayhi)
c     call arwmax(myq,mzq,elxy,elxz,emxy,emxz,
c    1            exlo,exhi,axlo,axhi)
c     ezhi = 42   ! added by Johnny
c     ezhi=ezhi*amp
c     azhi = 42   ! added by Johnny
c     azhi=azhi*amp
      eyhi=sqrt(eyhi)*amp
      ayhi=sqrt(ayhi)*amp
c     exhi = 42   ! added by Johnny
c     exhi=exhi*amp*50.
c     axhi = 42   ! added by Johnny
c     axhi=axhi*amp
      write(6,*)'ezhi,azhi,eyhi,ayhi,exhi,axhi=',
     1  ezhi,azhi,eyhi,ayhi,exhi,axhi
c
      do 377 i = 1, mxr1
      do 377 j = 1, mzr1
      el=sqrt(elyx(i,j)**2+elyz(i,j)**2)
      em=sqrt(emyx(i,j)**2+emyz(i,j)**2)
      elyx(i,j)=elyx(i,j)/(el+eyhi)
      elyz(i,j)=elyz(i,j)/(el+eyhi)
      emyx(i,j)=emyx(i,j)/(em+ayhi)
      emyz(i,j)=emyz(i,j)/(em+ayhi)
c     elyx(i,j)=elyx(i,j)/(abs(elyx(i,j))+eyhi)
c     elyz(i,j)=elyz(i,j)/(abs(elyz(i,j))+eyhi)
c     emyx(i,j)=emyx(i,j)/(abs(emyx(i,j))+ayhi)
c     emyz(i,j)=emyz(i,j)/(abs(emyz(i,j))+ayhi)
  377 continue
c     do 372 i = 1, myq
c     do 372 j = 1, mzq
c     el=sqrt(elxy(i,j)**2+elxz(i,j)**2)
c     em=sqrt(emxy(i,j)**2+emxz(i,j)**2)
c     elxy(i,j)=elxy(i,j)/(el+exhi)
c     elxz(i,j)=elxz(i,j)/(el+exhi)
c attention
c     emxy(i,j)=emxy(i,j)/(em+axhi)
c     emxz(i,j)=emxz(i,j)/(em+axhi)
c 372 continue
c
c
c     call arwmax(mxr1,myr1,elzx,elzy,emzx,emzy,
c    1            ezlo,ezhi,azlo,azhi)
      call arwmax(mxr1,mzr1,elyx,elyz,emyx,emyz,
     1            eylo,eyhi,aylo,ayhi)
c     call arwmax(myq,mzq,elxy,elxz,emxy,emxz,
c    1            exlo,exhi,axlo,axhi)
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
c     do 8 k=1,mzs
c     do 8 j=1,mys
c     do 8 i=1,mxs
c     rhi(i,j,k)=0.
c     rhe(i,j,k)=0.
c     eke(i,j,k)=0.
c     fex(i,j,k)=0.
c     fey(i,j,k)=0.
c     fez(i,j,k)=0.
c     fix(i,j,k)=0.
c     fiy(i,j,k)=0.
c     fiz(i,j,k)=0.
c     rhij(i,j,k)=0.
c     rhej(i,j,k)=0.
c     fexj(i,j,k)=0.
c     feyj(i,j,k)=0.
c     fezj(i,j,k)=0.
c     fixj(i,j,k)=0.
c     fiyj(i,j,k)=0.
c     fizj(i,j,k)=0.
c   8 continue

C  Volume weighting:
c for ions
c     call Vol_Weighting(mxs,mys,mzs,iong,rhi,eki,uig,vig,wig,
c    &                         xig,yig,zig,fix,fiy,fiz)
c
c      write(6,*)'after background ions'

c for electrons
c     call Vol_Weighting(mxs,mys,mzs,lecg,rhe,eke,ueg,veg,weg,
c    &                         xeg,yeg,zeg,fex,fey,fez)
c      write(6,*)'after background electrons'
c       write(6,*)' non-zero rhe, zero rhej'
c     write(6,*)'rhe(350,25,25),rhej(350,25,25)=',
c    1 rhe(350,25,25),rhej(350,25,25)

C  Volume weighting:
c for jet ions
c     call Vol_Weighting(mxs,mys,mzs,ionjg,rhij,ekij,uijg,vijg,wijg,
c    &                         xijg,yijg,zijg,fixj,fiyj,fizj)
c     write(6,*)'after jet ions'

c
c for jet electrons
c     call Vol_Weighting(mxs,mys,mzs,lecjg,rhej,ekej,uejg,vejg,wejg,
c    &                         xejg,yejg,zejg,fexj,feyj,fezj)

cheung
c     write(6,*)'after jet electrons'

c       write(6,*)' non-zero rhe, non-zero rhej'
c     write(6,*)'rhe(350,25,25),rhej(350,25,25)=',
c    1 rhe(350,25,25),rhej(350,25,25)
c     s = 42    !  added by Johnny
c     cx = 42   !  added by Johnny
c     dx = 42   !  added by Johnny
c     cy = 42   !  added by Johnny
c     dy = 42   !  added by Johnny
c     cz = 42   !  added by Johnny
c     dz = 42   !  added by Johnny
c     write(6,*)"s cx dx cy dy cz dz", s,cx,dx,cy,dy,cz,dz

c
c calculating Lorentz factor in the jet
c     61 < y < 65 yjc =503
c
c for jet electron Lorentz factor
c      print *, 'before rhe eke  ,lecjg=',lecjg
c     do n0 = 1,lecjg

c since fluid quantiies are calculated the following lins are deleted!
c
c  attention the sustem is original size  
c
c attention z --->x  x --> z
      do 11 j=1,mzk
c        ii=2*(j-1)+7
         ii=2*(j-1)+257
      do 11 i=1,mxk
c        jj=2*(i-1)+12+ishift
         jj=2*(i-1)+101+ishift
c     densi(i,j)=rhe(jj,jn,ii)
c     dense(i,j)=rhe(jj,jn,ii)+rhej(jj,jn,ii)
c calculationg Lorentz factor
      if(rhij(jj,jn,ii).eq.0.0) then
      dense(i,j)=1.0
      else
c     dense(i,j)=eke(jj,jn,ii)/rhe(jj,jn,ii)
c     dense(i,j)=0.5*(fixj(jj,jn,ii)**2+fiyj(jj,jn,ii)**2
c    +  +fizj(jj,jn,ii)**2)/rhe(jj,jn,ii)**2
      dense(i,j)=ekij(jj,jn,ii)/rhij(jj,jn,ii)
       end if
c     btot(i,j)=sqrt(
c    1  (0.5*(bxx(jj,jn,ii) +bxx(jj,jn-1,ii-1)))**2
c    1 +(0.5*(byy(jj,jn,ii) +byy(jj-1,jn,ii-1)))**2
c    1 +(0.5*(bzz(jj,jn,ii) +bzz(jj-1,jn-1,ii)))**2)
c     etot(i,j)=sqrt(
c    1  (0.5*(exx(jj,jn,ii) +exx(jj-1,jn,ii)))**2
c    1 +(0.5*(eyy(jj,jn,ii) +eyy(jj,jn-1,ii)))**2
c    1 +(0.5*(ezz(jj,jn,ii) +ezz(jj,jn,ii-1)))**2)
c attention   ajy = Jx
c*    ajy(i,j)=bxx(jj,jn,ii)
c     ajy(i,j)=0.5*(fix(jj,jn,ii)-fex(jj,jn,ii))
c    1         +0.25*(fix(jj,jn,ii-1)-fex(jj,jn,ii-1))
c    1         +0.25*(fix(jj,jn,ii+1)-fex(jj,jn,ii+1))
c    1        +0.5*(fixj(jj,jn,ii)-fexj(jj,jn,ii))
c    1         +0.25*(fixj(jj,jn,ii-1)-fexj(jj,jn,ii-1))
c    1         +0.25*(fixj(jj,jn,ii+1)-fexj(jj,jn,ii+1))
c     vtemp= 1.0 -(fex(jj,jn,ii)**2+fey(jj,jn,ii)**2
c    1      +fez(jj,jn,ii)**2)/(c*rhe(jj,jn,ii))**2
c    1           -(fexj(jj,jn,ii)**2+feyj(jj,jn,ii)**2
c    1      +fezj(jj,jn,ii)**2)/(c*rhej(jj,jn,ii))**2
c     if(vtemp.le.0.0) vtemp=0.0001
c     vgam(i,j)=1.0/sqrt(vtemp)
   11 continue

c      print *,'dense(300,35),densi(300,35),ajy(300,35)=',
c    1  dense(300,35),densi(300,35),ajy(300,35)
c      write(6,*)'dense(300,35),densi(300,35),ajy(300,35)=',
c    1  dense(300,35),densi(300,35),ajy(300,35)

      do 311 j=1,mzr
c     jj=26*j -12
      jj=24*j +229
      do 311 i=1,mxr
c     ii=26*i -12 +ishift
      ii=24*i +67
      fvliy(i,j)=0.5*fiy(ii,jn,jj)
     1        +0.25*(fiy(ii,jn,jj-1)+fiy(ii,jn,jj+1))
     1          +0.5*fiyj(ii,jn,jj)
     1        +0.25*(fiyj(ii,jn,jj-1)+fiyj(ii,jn,jj+1))
      fvlix(i,j)=0.5*fix(ii,jn,jj)
     1        +0.25*(fix(ii,jn,jj-1)+fix(ii,jn,jj+1))
     1          +0.5*fixj(ii,jn,jj)
     1        +0.25*(fixj(ii,jn,jj-1)+fixj(ii,jn,jj+1))
      fvley(i,j)=0.5*fey(ii,jn,jj)
     1        +0.25*(fey(ii,jn,jj-1)+fey(ii,jn,jj+1))
     1          +0.5*feyj(ii,jn,jj)
     1        +0.25*(feyj(ii,jn,jj-1)+feyj(ii,jn,jj+1))
      fvlex(i,j)=0.5*fex(ii,jn,jj)
     1        +0.25*(fex(ii,jn,jj-1)+fex(ii,jn,jj+1))
     1          +0.5*fexj(ii,jn,jj)
     1        +0.25*(fexj(ii,jn,jj-1)+fexj(ii,jn,jj+1))
  311 continue
c
      call aminmax(mxk,mzk,densi,dense,
     1             gmmaxi,gmmini,gmmaxe,gmmine)
      call aminmax(mxk,mzk,btot,etot,    
     1             gmmaxbt,gmminbt,gmmaxet,gmminet)
      call aminmax(mxk,mzk,ajy,vgam,
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
      call aminmax(mxk,mzk,densi,dense,
     1             gmmaxi,gmmini,gmmaxe,gmmine)
c
c
 7139 continue
c
c
      go to 3927
      amp=0.01
      call arwmax(mxr,mzr,fvlex,fvley,fvlix,fvliy,
     1            flo1,hi1,flo2,hi2)
      hi1=hi1*amp
      hi2=hi2*amp
c
      do 37 i = 1, mxr
      do 37 j = 1, mzr
      fve=sqrt(fvlex(i,j)**2+fvley(i,j)**2)
      fvi=sqrt(fvlix(i,j)**2+fvliy(i,j)**2)
      fvlex(i,j)=fvlex(i,j)/(fve+hi1)
      fvley(i,j)=fvley(i,j)/(fve+hi1)
      fvlix(i,j)=fvlix(i,j)/(fvi+hi2)
      fvliy(i,j)=fvliy(i,j)/(fvi+hi2)
   37 continue
c
c3927 continue
c
      call arwmax(mxr,mzr,fvlex,fvley,fvlix,fvliy,
     1            flo1,hi1,flo2,hi2)

c     do 215 k=1,mzq
c     kk=3*k -2
c     do 215 i=1,mxq
c     ii=3*i -2
c attention for notation x, y
c     fvqix(i,k)=fiz(ii,jn,kk)
c     fvqiy(i,k)=fix(ii,jn,kk)
c     fvqex(i,k)=fez(ii,jn,kk)
c     fvqey(i,k)=fex(ii,jn,kk)
c 215 continue
c
c     do 205 j=1,mzk
c        ii=2*(j-1)+1
c     do 205 i=1,mxk
c        jj=2*(i-1)+1+ishift
c     denqi(i,j)=eki(jj,jn,ii)
c     denqe(i,j)=eke(jj,jn,ii)
c attention for notation x, y
c     aey(i,j)  =0.5*(eyy(jj,jn,ii) +eyy(jj,jn-1,ii))
c 205 continue
c
c     call aminmax(mxk,mzk,denqi,denqe,
c    1             gmmaxki,gmminki,gmmaxke,gmminke)
c     call aminmax(mxk,mzk,aey,aey,   
c    1             gmmaxey,gmminey,gmmaxey,gmminey)
c     call arwmax(mxq,mzq,fvqex,fvqey,fvqix,fvqiy,
c    1            flo1,hi1,flo2,hi2)
 3927 continue
      t=nstep-1
c
      nset=-1
      length=0
      ispv=0
      spv=1.0
c attention
c     xmin=16.0  +ishift
      xmin=100.0  +ishift
c     xmax=152.0
c     xmax=148.0
c     xmax=1276.0 +ishift
      xmax=1120.0 +ishift
c     ymin=7.0 
      ymin=256.0
c     ymax=77.0
c     ymax=73.0
c New Ken
c     xmin = 8.0 +float(ishift)
c     xmax = 408.0 +float(ishift)
c     ymax=782.0
      ymax=506.0
c     zmin=7.0
      zmin=256.0
c     zmax=77.0 
c     zmax=782.0 
      zmax=506.0
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

       print *, 'before ICHOS'
      ICHOS=111
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
      if(ichos.eq.214)
     .write(message,214) t


c     write(message,113) t                                               

        CALL DFCLRS
        CALL GSTXCI (1)
        CALL GSPLCI (1)

c     if(ichos.eq.110.or.ichos.eq.214)
c    .call colorcn(ajy,mxk,mzk,gmminjy,gmmaxjy,gridr,message,2)

      print *, 'before calling dense'
      print *, 'dense(500,125)=',dense(500,125)
      print *, 'mxk,mzk,gmmine,gmmaxe=',mxk,mzk,gmmine,gmmaxe
      if(ichos.eq.111)
     .call colorcn(dense,mxk,mzk,gmmine,gmmaxe,gridr,message,1)



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

      call set(xvpl,xvpr,yvpb,yvpt,0.0,1000.0,0.0,500.0,1)

      IOFFM = 1

      print *, 'before emyx'
      if(ichos.eq.111.or.ichos.eq.214)
     .call velvct(emyx,mxr1,emyz,mxr1,mxr1,mzr1,aylo,ayhi,nset,length,
     &            ispv,spv)


       CALL DFCLRS
        CALL GSTXCI (0)
        CALL GSPLCI (0)

      call dscale(gridr,bnu1,xmin,xmax,zmin,zmax, 
     1            0.04,-0.03,'Y$','Z$',message)
      xpob=0.525 
      ypob=0.340
c     xpol=0.205
c     xpol=0.105
      xpol=0.050
c     ypol=0.47
c     ypol=0.5925
      ypol=0.4925
      size=0.018
c  for x - y
c     ichs=1
c  for x - z
      ichs=2
c  for y - z
c     ichs=3
      call eplch(xpob,ypob,xpol,ypol,size,ichs)
       CALL DFCLRS
        CALL GSTXCI (1)
        CALL GSPLCI (1)
      call frame
 8173 continue
       enddo

c 110 format('JY      (ELE)  T=',f8.1)
  110 format('JX      (EFLU) T=',f8.1)
  111 format('LOR_Ele (MAG)  T=',f8.1)
c 112 format('ELE DEN (EFLU) T=',f8.1)
  112 format('ELE DEN (ELE)  T=',f8.1)
  113 format('ION DEN (IFLU) T=',f8.1)
  114 format('BTOT    (MAG)  T=',f8.1)
  115 format('ETOT    (ELE)  T=',f8.1)
  116 format('GAMMAE  (EFLU) T=',f8.1)
  117 format('EY      (ELE)  T=',f8.1)
  118 format('EKIN    (EFLU) T=',f8.1)
  119 format('IKIN    (IFLU) T=',f8.1)
  214 format('CUR DEN (MAG)  T=',f8.1)
  215 format('CUR DEN (ELE)  T=',f8.1)
  275 format('CUR DEN (FLUX) T=',f8.1)
  216 format('CUR DEN (CUR)  T=',f8.1)
  223 format('ELE ENE (MAG)  T=',f8.1)
  224 format('ION ENE (MAG)  T=',f8.1)
  323 format('ELE ENE (FLUX) T=',f8.1)
  324 format('ION ENE (FLUX) T=',f8.1)
   25 format('X-Y SPACE (ION) T=',F8.1,'$')                            
   26 format('X-Y SPACE (ELE) T=',F8.1,'$')                            
   27 format('Z-VZ  SPACE     T=',F8.1,'$')                            
c                                                                       
c     go to 4
      call clsgks
      stop
      end
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
c       parameter (nspc=266000,nwsp=3882000)
        parameter (nspc=566000,nwsp=9882000)
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

