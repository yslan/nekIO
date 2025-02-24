c     imode=1: onto uniform grid excluding elements boundary
c     imode=2: uniform grid, no element interface, yes bdry
c              deform boudary elements in usrdat2
c              set (nelx, nely, nelz)
c              call box_outfld_rescale_mesh
c     nrg, new lx1 (# grid in 1D) for output file
c-----------------------------------------------------------------------
      subroutine myoutpost_box(v1,v2,v3,vp,vt,name3,imode)

      include 'SIZE'
      include 'INPUT'

      real v1(1),v2(1),v3(1),vp(1),vt(1)
      character*3 name3

      itmp=0
      if (ifto) itmp=1
      call myoutpost2_box(v1,v2,v3,vp,vt,itmp,name3,imode)

      return
      end
c-----------------------------------------------------------------------
      subroutine myoutpost2_box(v1,v2,v3,vp,vt,nfldt,name3,imode)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'

      parameter(ltot1=lx1*ly1*lz1*lelt)
      parameter(ltot2=lx2*ly2*lz2*lelv)
      common /outtmp/  w1(ltot1),w2(ltot1),w3(ltot1),wp(ltot2)
     &                ,wt(ltot1,ldimt)
c
      real v1(1),v2(1),v3(1),vp(1),vt(ltot1,1)
      character*3 name3
      logical if_save(ldimt)
      integer nrg_s

      nrg_s = nrg ! used and modified
      call setup_interp(imode)
c
      ntot1  = lx1*ly1*lz1*nelt
      ntot1t = lx1*ly1*lz1*nelt
      ntot2  = lx2*ly2*lz2*nelt

      if(nfldt.gt.ldimt) then
        write(6,*) 'ABORT: outpost data too large (nfldt>ldimt)!'
        call exitt
      endif

c store solution
      call copy(w1,vx,ntot1)
      call copy(w2,vy,ntot1)
      call copy(w3,vz,ntot1)
      call copy(wp,pr,ntot2)
      do i = 1,nfldt
         call copy(wt(1,i),t(1,1,1,1,i),ntot1t)
      enddo

c swap with data to dump
      call copy(vx,v1,ntot1)
      call copy(vy,v2,ntot1)
      call copy(vz,v3,ntot1)
      call copy(pr,vp,ntot2)
      do i = 1,nfldt
         call copy(t(1,1,1,1,i),vt(1,i),ntot1t)
      enddo

c dump data
      if_save(1) = ifto
      ifto = .false.
      if(nfldt.gt.0) ifto = .true.
      do i = 1,ldimt-1
         if_save(i+1) = ifpsco(i)
         ifpsco(i) = .false.
         if(i+1.le.nfldt) ifpsco(i) = .true.
      enddo

      call myprepost_box(.true.,name3,imode)

      ifto = if_save(1)
      do i = 1,ldimt-1
         ifpsco(i) = if_save(i+1)
      enddo

c restore solution data
      call copy(vx,w1,ntot1)
      call copy(vy,w2,ntot1)
      call copy(vz,w3,ntot1)
      call copy(pr,wp,ntot2)
      do i = 1,nfldt
         call copy(t(1,1,1,1,i),wt(1,i),ntot1t)
      enddo

      nrg = nrs_s
      return
      end
c-----------------------------------------------------------------------
      subroutine myprepost_box(ifdoin,prefin,imode)

c     Store results for later postprocessing
c
c     Recent updates:
c
c     p65 now indicates the number of parallel i/o files; iff p66 >= 6
c
c     we now check whether we are going to checkpoint in set_outfld
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)

      character*3    prefin,prefix

      logical  ifdoin

      if (ioinfodmp.eq.-2) return

c#ifdef TIMER
      etime1=dnekclock_sync()
c#endif

      prefix = prefin
      if (prefix.eq.'his') prefix = '   '

      if (ifdoin) then
         icalld=icalld+1
         nprep=icalld

         call prepost_map(0) ! map pr and axisymm. arrays
         call myoutfld_box(prefix,imode)
         call prepost_map(1) ! map back axisymm. arrays

c#ifdef TIMER
         tprep=tprep+dnekclock_sync()-etime1
c#endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine myoutfld_box(prefix,imode)

c     output .fld file

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
C
C     Work arrays and temporary arrays
C
      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)
c
c     note, this usage of CTMP1 will be less than elsewhere if NELT ~> 3.
      parameter (lxyz=lx1*ly1*lz1)
      parameter (lpsc9=ldimt1+9)
      common /ctmp1/ tdump(lxyz,lpsc9)
      real*4         tdump
      real           tdmp(4)
      equivalence   (tdump,tdmp)

      real*4         test_pattern

      character prefix*(*)

      character*1    fhdfle1(132)
      character*132   fhdfle
      equivalence   (fhdfle,fhdfle1)

      character*1 excode(30)
      character*10 frmat

      common /nopenf/ nopen(99)

      common /rdump/ ntdump
      data ndumps / 0 /

      logical ifxyo_s

      if(nio.eq.0) then
        WRITE(6,1001) istep,time
 1001   FORMAT(/,i9,1pe12.4,' Write checkpoint')
      endif
      call nekgsync()

      p66 = param(66)
      if (abs(p66).eq.6) then
         call mymfo_outfld_box(prefix,imode)
         return
      endif
      call exitti('myoutfld_box only support p66=6 for now!$',1)

      return
      end
c-----------------------------------------------------------------------
      subroutine mymfo_outfld_box(prefix,imode)  ! muti-file output

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)  ! mapped pressure

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzo8
      character prefix*(*)
      logical ifxyo_s, ifinterp

      integer cnt, cntg

      common /SCRUZ/  ur1(lxo*lxo*lxo*lelt)
     &              , ur2(lxo*lxo*lxo*lelt)
     &              , ur3(lxo*lxo*lxo*lelt)

      tiostart=dnekclock_sync()

      call io_init

      ifxyo_s = ifxyo
      ifxyo_  = ifxyo
      nout = nelt ! ump all fields based on the t-mesh to avoid different topologies
      nxo  = lx1
      nyo  = ly1
      nzo  = lz1
      ifinterp = .false.
      if (imode.ne.0) ifinterp = .true.
      if (ifinterp) then ! do interpolation
         if (nrg.gt.lxo) then
            if (nid.eq.0) write(6,*)
     &         'WARNING: nrg too large, reset to lxo!'
            nrg = lxo
         endif
         nxo  = nrg
         nyo  = nrg
         nzo  = 1
         if(if3d) nzo = nrg
      endif

      cnt = 0
      do iel = 1,nelt
        if(out_mask(iel).ne.0) cnt = cnt + 1
      enddo
      cntg = iglsum(cnt, 1)

      ierr=0
      if (nid.eq.pid0) then
         call mymfo_open_files_box(prefix,imode,ierr)         ! open files on i/o nodes
      endif
      call err_chk(ierr,'Error opening file in mfo_open_files. $')
      call bcast(ifxyo_,lsize)
      ifxyo = ifxyo_

      call blank(rdcode1,10)
      i = 1
      IF (IFXYO) THEN
         rdcode1(i)='X'
         i = i + 1
      ENDIF
      IF (IFVO) THEN
         rdcode1(i)='U'
         i = i + 1
      ENDIF
      IF (IFPO) THEN
         rdcode1(i)='P'
         i = i + 1
      ENDIF
      IF (IFTO) THEN
         rdcode1(i)='T'
         i = i + 1
      ENDIF
      IF (LDIMT.GT.1) THEN
         NPSCALO = 0
         do k = 1,ldimt-1
           if(ifpsco(k)) NPSCALO = NPSCALO + 1
         enddo
         IF (NPSCALO.GT.0) THEN
            rdcode1(i) = 'S'
            WRITE(rdcode1(i+1),'(I1)') NPSCALO/10
            WRITE(rdcode1(i+2),'(I1)') NPSCALO-(NPSCALO/10)*10
         ENDIF
      ENDIF

      call mfo_write_hdr(rdcode1) ! including element mapping

      nxyzo8  = nxo*nyo*nzo

      ! only relevant for single shared file
      offs0 = iHeaderSize + 4 + isize*cntg
      stride  = cntg* nxyzo8*wdsizo
      strideB = nelB * nxyzo8*wdsizo

      ioflds = 0
      if (ifxyo) then
         offs = offs0 + ldim*strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifinterp) then
            call mymap2reg_box(ur1,nrg,xm1,nout,imode)
            call mymap2reg_box(ur2,nrg,ym1,nout,imode)
            if (if3d) call mymap2reg_box(ur3,nrg,zm1,nout,imode)
            call mfo_outv(ur1,ur2,ur3,nout,nxo,nyo,nzo)
         else
            call mfo_outv(xm1,ym1,zm1,nout,nxo,nyo,nzo)
         endif
         ioflds = ioflds + ldim
      endif
      if (ifvo ) then
         offs = offs0 + ioflds*stride + ldim*strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifinterp) then
             call mymap2reg_box(ur1,nrg,vx,nout,imode)
             call mymap2reg_box(ur2,nrg,vy,nout,imode)
             if (if3d) call mymap2reg_box(ur3,nrg,vz,nout,imode)
             call mfo_outv(ur1,ur2,ur3,nout,nxo,nyo,nzo)
         else
            call mfo_outv(vx,vy,vz,nout,nxo,nyo,nzo)  ! B-field handled thru outpost
         endif
         ioflds = ioflds + ldim
      endif
      if (ifpo ) then
         offs = offs0 + ioflds*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifinterp) then
            call mymap2reg_box(ur1,nrg,pm1,nout,imode)
            call mfo_outs(ur1,nout,nxo,nyo,nzo)
         else
            call mfo_outs(pm1,nout,nxo,nyo,nzo)
         endif
         ioflds = ioflds + 1
      endif
      if (ifto ) then
         offs = offs0 + ioflds*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifinterp) then
            call mymap2reg_box(ur1,nrg,t,nout,imode)
            call mfo_outs(ur1,nout,nxo,nyo,nzo)
         else
            call mfo_outs(t,nout,nxo,nyo,nzo)
         endif
         ioflds = ioflds + 1
      endif
      do k=1,ldimt-1
         if(ifpsco(k)) then
           offs = offs0 + ioflds*stride + strideB
           call byte_set_view(offs,ifh_mbyte)
           if (ifinterp) then
              call mymap2reg_box(ur1,nrg,t(1,1,1,1,k+1),nout,imode)
              call mfo_outs(ur1,nout,nxo,nyo,nzo)
           else
              call mfo_outs(t(1,1,1,1,k+1),nout,nxo,nyo,nzo)
           endif
           ioflds = ioflds + 1
         endif
      enddo
      dnbyte = 1.*ioflds*cnt*wdsizo*nxo*nyo*nzo

      ! add meta data to the end of the file
      if (if3d) then
         offs0 = offs0 + ioflds*stride
         stride  = cntg *2*4
         strideB = nelB *2*4
         ioflds  = 0
         if (ifxyo) then
            offs = offs0 + ioflds*stride + ldim*strideB
            call byte_set_view(offs,ifh_mbyte)
            call mfo_mdatav(xm1,ym1,zm1,nout)
            ioflds = ioflds + ldim
         endif
         if (ifvo ) then
            offs = offs0 + ioflds*stride + ldim*strideB
            call byte_set_view(offs,ifh_mbyte)
            call mfo_mdatav(vx,vy,vz,nout)
            ioflds = ioflds + ldim
         endif
         if (ifpo ) then
            offs = offs0 + ioflds*stride + strideB
            call byte_set_view(offs,ifh_mbyte)
            call mfo_mdatas(pm1,nout)
            ioflds = ioflds + 1
         endif
         if (ifto ) then
            offs = offs0 + ioflds*stride + strideB
            call byte_set_view(offs,ifh_mbyte)
            call mfo_mdatas(t,nout)
            ioflds = ioflds + 1
         endif
         do k=1,ldimt-1
            offs = offs0 + ioflds*stride + strideB
            call byte_set_view(offs,ifh_mbyte)
            if(ifpsco(k)) call mfo_mdatas(t(1,1,1,1,k+1),nout)
            ioflds = ioflds + 1
         enddo
         dnbyte = dnbyte + 2.*ioflds*cnt*wdsizo
      endif

      ierr = 0
      if (nid.eq.pid0) then
         if(ifmpiio) then
           call byte_close_mpi(ifh_mbyte,ierr)
         else
           call byte_close(ierr)
         endif
      endif
      call err_chk(ierr,'Error closing file in mfo_outfld. Abort. $')

      tio = dnekclock_sync()-tiostart
      if (tio.le.0) tio=1.

      dnbyte = glsum(dnbyte,1)
      dnbyte = dnbyte + iHeaderSize + 4. + isize*cntg
      dnbyte = dnbyte/1e9
      if(nio.eq.0) write(6,7) istep,time,dnbyte,dnbyte/tio,
     &             nfileo
    7 format(/,i9,1pe12.4,' done :: Write checkpoint',/,
     &       30X,'file size = ',3pG12.2,'GB',/,
     &       30X,'avg data-throughput = ',0pf7.1,'GB/s',/,
     &       30X,'io-nodes = ',i5,/)

      ifxyo = ifxyo_s ! restore old value

      return
      end
c-----------------------------------------------------------------------
      subroutine mymfo_open_files_box(prefix,imode,ierr) ! open files

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      character prefix*(*)

      character*3 prefx3
      data        prefx3 / "   " /

      character*132  name
      character*1    nam1(132)
      equivalence   (nam1,name)

      character*132  fname
      character*1    fnam1(132)
      equivalence   (fnam1,fname)

      character*6  six,str
      save         six
      data         six / "??????" /


      character*1 slash,dot
      save        slash,dot
      data        slash,dot  / '/' , '.' /

      integer nopen(1000,3)
      save    nopen
      data    nopen  / 3000*0 /

      call blank(fname,132)
      call blank(name,132)

      iprefix = i_find_prefix(prefix,1000)
      if (imode.eq.1) then
         nopen(iprefix,2) = nopen(iprefix,2)+1
         nfld             = nopen(iprefix,2)
      elseif (imode.eq.2) then
         nopen(iprefix,3) = nopen(iprefix,3)+1
         nfld             = nopen(iprefix,3)
      else
         nopen(iprefix,1) = nopen(iprefix,1)+1
         nfld             = nopen(iprefix,1)
      endif

      call chcopy(prefx3,prefix,min(len(prefix),3))        ! check for full-restart request
      if (prefx3.eq.'rst'.and.max_rst.gt.0) nfld = mod1(nfld,max_rst)

      call restart_nfld(nfld, prefix) ! Check for Restart option.
      if (prefx3.eq.'   '.and.nfld.eq.1) ifxyo_ = .true. ! 1st file

      if(ifmpiio) then
        rfileo = 1
      else
        rfileo = nfileo
      endif
      ndigit = log10(rfileo) + 1

      lenp = 0 !ltrunc(path,132)
      call chcopy(fnam1(1),path,lenp)
      k = 1 + lenp

      call blank(name,132)
      kk = 1

      if (ifdiro) then                                  !  Add directory
         call chcopy(fnam1(k),'A',1)
         k = k + 1
         call chcopy(fnam1(k),six,ndigit)  ! put ???? in string
         k = k + ndigit
         call chcopy(fnam1(k),slash,1)
         k = k + 1
      endif

      if (len(prefix) .gt. 0) then !  Add prefix
         if(prefx3 .ne. '   ') then
           call chcopy(fnam1(k),prefix,len(prefix))
           k = k + len(prefix)

           call chcopy(nam1(kk),prefix,len(prefix))
           kk = kk + len(prefix)
         endif
      endif

      ll=ltrunc(session,132)
      call chcopy(fnam1(k),session,ll)
      k = k+ll
      call chcopy(nam1(kk),session,ll)
      kk = kk + ll

      if (imode.eq.1) then
         ll=4
         call chcopy(fnam1(k),'_ubx',ll)
         k = k + ll
         call chcopy(nam1(kk),'_ubx',ll)
         kk = kk + ll
      elseif (imode.eq.2) then
         ll=4
         call chcopy(fnam1(k),'_bbx',ll)
         k = k + ll
         call chcopy(nam1(kk),'_bbx',ll)
         kk = kk + ll
      endif

      if(nid.eq.0) then
        call chcopy(nam1(kk),'.nek5000',8)
        open(unit=101,file=name)
        write(101,*) "filetemplate: ", trim(fname), "%01d.f%05d"
        write(101,*) "firsttimestep: 1"
        write(101,*) "numtimesteps: ", nfld
        close(101)
      endif

      call chcopy(fnam1(k),six,ndigit)                  !  Add file-id holder
      k = k + ndigit

      call chcopy(fnam1(k  ),dot,1)                     !  Add .f appendix
      call chcopy(fnam1(k+1),'f',1)
      k = k + 2

      write(str,4) nfld                                 !  Add nfld number
    4 format(i5.5)
      call chcopy(fnam1(k),str,5)
      k = k + 5

      call addfid(fname,fid0)

      if(ifmpiio) then
        if(nio.eq.0)    write(6,*) '      FILE:',fname
        call byte_open_mpi(fname,ifh_mbyte,.false.,ierr)
      else
        if(nid.eq.pid0) write(6,*) '      FILE:',fname
        call byte_open(fname,ierr)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mymap2reg_box(ur,n,u,nel,imode)
c
c     Map scalar field u() to regular n x n x n array ur
c
      implicit none
      include 'SIZE'
      integer n, nel, imode
      real ur(1),u(lx1*ly1*lz1,1)

      if (imode.eq.1) then
         call mymap2reg_box1(ur,n,u,nel)
      elseif (imode.eq.2) then
         call mymap2reg_box2(ur,n,u,nel)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mymap2reg_box1(ur,n,u,nel)
      implicit none
      include 'SIZE'

      real ur(1),u(lx1*ly1*lz1,1)

      integer n,nel,e,k,ldr,l
      parameter (l=50)

      real j, jt, w, z
      common /my_box_interp1/ j(l*l),jt(l*l),w(l*l),z(l)

      ldr = n**ldim

      k=1
      do e=1,nel
         if (ldim.eq.2) call apply_tensor2(ur(k),n,u(1,e),lx1,j,jt)
         if (ldim.eq.3) call apply_tensor3(ur(k),n,u(1,e),lx1,j,jt)
         k = k + ldr
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mymap2reg_box2(ur,n,u,nel)
      implicit none
      include 'SIZE'
      real ur(1),u(lx1*ly1*lz1,1)

      integer n,nel,e,k,ldr,l
      parameter (l=50)

      integer etype, j1,j2,j3
      real jx, jxt, we, ze
      common /my_box_interp2/ jx(l*l,3),jxt(l*l,3),we(l*l),ze(l)
      common /my_box_interp2i/ etype(3,lelt)

      ldr = n**ldim

      k=1
      do e=1,nel
         j1 = etype(1,e) ! different spacing in each direction
         j2 = etype(2,e)

         if (ldim.eq.2)
     $      call apply_tensor2b(ur(k),n,u(1,e),lx1
     $                         ,jx(1,j1),jxt(1,j2))
         if (ldim.eq.3) then
            j3 = etype(3,e)
            call apply_tensor3b(ur(k),n,u(1,e),lx1
     $                         ,jx(1,j1),jxt(1,j2),jxt(1,j3))
         endif
         k = k + ldr
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_tensor3(uf,n,uc,m,j,jt) ! uc -> uf
      implicit none
      integer m, n, mm, mn, nn, iv, iw, k
      real uf(n,n,n),uc(m,m,m),j(1),jt(1)

      integer l
      parameter (l=50)
      real v, w
      common /box_scratch/ v(l*l*l),w(l*l*l)

      mm = m*m
      mn = m*n
      nn = n*n

      call mxm(j,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxm(v(iv),n,jt,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxm(w,nn,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_tensor2(uf,n,uc,m,j,jt) ! uc -> uf
      implicit none
      integer m, n
      real uf(n,n),uc(m,m),j(1),jt(1)

      integer l
      parameter (l=50)
      real v, w
      common /box_scratch/ v(l*l*l),w(l*l*l)

      call mxm(j,n,uc,m,w ,m)
      call mxm(w,n,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_tensor3b(uf,n,uc,m,j1,j2t,j3t) ! uc -> uf
      implicit none
      integer m, n, mm, mn, nn, iv, iw, k
      real uf(n,n,n),uc(m,m,m),j1(1),j2t(1),j3t(1)

      integer l
      parameter (l=50)
      real v, w
      common /box_scratch/ v(l*l*l),w(l*l*l)

      mm = m*m
      mn = m*n
      nn = n*n

      call mxm(j1,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxm(v(iv),n,j2t,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxm(w,nn,j3t,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_tensor2b(uf,n,uc,m,j1,j2t) ! uc -> uf
      implicit none
      integer m, n
      real uf(n,n),uc(m,m),j1(1),j2t(1)

      integer l
      parameter (l=50)
      real v, w
      common /box_scratch/ v(l*l*l),w(l*l*l)

      call mxm(j1,n,uc,m,w ,m)
      call mxm(w,n,j2t,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_interp(imode)
      implicit none
      include 'SIZE'
      include 'RESTART' ! nrg

      integer imode, ierr

      logical yes_setup(2)
      data yes_setup / .true., .true. /
      save yes_setup

      integer nrg_prev(2)
      data nrg_prev /0, 0/

      if (imode.ne.1.AND.imode.ne.2)
     $   call exitti('box_interp invalid imode$',imode)

      if (nrg.eq.0) nrg = lx1
      if (nrg.ne.nrg_prev(imode)) yes_setup(imode) = .true.

      if (yes_setup(imode)) then

         call check_unif_box_mesh(imode, ierr)
         if (imode.eq.1) then
            call setup_interp_matrix1(lx1,nrg)
         elseif (imode.eq.2) then
            call setup_interp_matrix2(lx1,nrg)
         endif
         if (nio.eq.0) write(*,*)'setup_box_interp',imode,lx1,nrg

         nrg_prev(imode) = nrg
         yes_setup(imode) = .false.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_interp_matrix1(ni,no)
      implicit none

      integer ni, no, l
      parameter (l=50)

      real j, jt, w, z
      common /my_box_interp1/ j(l*l),jt(l*l),w(l*l),z(l)

      if (no.gt.l) call exitti('box_interp1, memory$',no)

      call zwgll(z,w,ni)
      call zuni_int(w,no)
      call gen_int_gz(j,jt,w,no,z,ni)

      return
      end
c-----------------------------------------------------------------------
      subroutine zuni_int(z,nz)
      implicit none
c             -1                         1
c     (-1,1)   o---x-----x ... x-----x---o
c              |dx2| dx1 |     | dx1 |dx2|
      real z(1), dz, dz2
      integer i, nz

      dz = 2./nz
      dz2= 0.5*dz
      z(1) = -1. + dz2
      do i = 2,(nz-1)
         z(i) = z(i-1) + dz
      enddo
      z(nz) = 1.0 - dz2

      return
      end
c-----------------------------------------------------------------------
      subroutine zuni_int_left(z,nz)
c             -1                           1
c     [-1,1)   x-----x-----x ... x-----x---o
c              | dx1 | dx1 |     | dx1 |dx2|
      implicit none
      real z(1), dz, dz2
      integer i, nz

      dz = 2.0 / (nz-0.5)
      dz2 = 0.5*dz
      z(1) = -1.0
      do i = 2,(nz-1)
         z(i) = z(i-1) + dz
      enddo
      z(nz) = 1.0 - dz2

      return
      end
c-----------------------------------------------------------------------
      subroutine zuni_int_right(z,nz)
c             -1                     1
c     (-1,1]   o---x-----x ... x-----x
c              |dx2| dx1 |     | dx1 |
      implicit none
      real z(1), dz, dz2
      integer i, nz

      dz = 2.0 / (nz-0.5)
      dz2 = 0.5*dz
      z(1) = -1.0 + dz2
      do i = 2,(nz-1)
         z(i) = z(i-1) + dz
      enddo
      z(nz) = 1.0

      return
      end
c-----------------------------------------------------------------------
      subroutine check_unif_box_mesh(imode, ierr)
      implicit none
      include 'SIZE'
      include 'PARALLEL' ! nelgt
      include 'GEOM' ! xm1
      include 'ZPER' ! nelx
      real tol
      integer imode, ierr

      tol = 1e-6 ! single precision rea will suffer this...

      ierr = 0
      if (nelx.eq.0) ierr = 1
      if (nely.eq.0) ierr = 2
      if (ldim.eq.3) then
         if (nelz.eq.0) ierr = 3
      else
         nelz = 1
      endif
      if (ierr.gt.0) then
         if (imode.eq.1) then
            if (nio.eq.0) write(*,*)'box-warn: nelx,nely,nelz not set'
            return
         elseif (imode.eq.2) then
            call exitti('box-Err nelx,nely,nelz not set$',ierr)
         endif
      endif

      ierr = 0
      if (nelx*nely*nelz.ne.nelgt) ierr = 1
      if (ierr.ne.0) then
         if (imode.eq.1) then
            if (nio.eq.0) write(*,*)'box-warn: prod(nelx,y,z) ne nelgt'
     $                             ,nelx*nely*nelz
            return
         elseif (imode.eq.2) then
            call exitti('box-Err prod(nelx,y,z) ne nelgt$'
     $                 ,nelx*nely*nelz)
         endif
      endif

      ierr = 0
      if (imode.eq.1) then 
         call check_box1_length(xm1, nelx, tol, 'x', ierr)
         call check_box1_length(ym1, nely, tol, 'y', ierr)
         if (ldim.eq.3)
     $      call check_box1_length(zm1, nelz, tol, 'z', ierr)
         if (ierr.ne.0.AND.nio.eq.0)
     $      write(*,*)'box-warn: not an uniform boxmesh!!!'
      elseif (imode.eq.2) then
         call check_box2_length(tol, ierr)
         if (ierr.ne.0)
     $      call exitti('mesh distribusion check fails, abort$',ierr)
      endif

c      call check_box_orientation() ! single box from genbox will pass
      return
      end
c-----------------------------------------------------------------------
      subroutine check_box1_length(xx, nelx, tol, s1, ierr)
      implicit none
      include 'SIZE'
      real xx(lx1*ly1*lz1,lelt), tol
      integer nelx, e, nxyz, ierr
      character*1 s1

      real vlmin,vlmax,glmin,glmax
      real dx,dxmin,dxmax,dxmid,dxref

      nxyz = lx1*ly1*lz1

      dxmin = 1e30
      dxmax = -1e30
      do e=1,nelt
         dx = vlmax(xx(1,e), nxyz) - vlmin(xx(1,e), nxyz)
         dxmin = min(dxmin, dx)
         dxmax = max(dxmax, dx)
      enddo

      dxmin = glmin(dxmin,1)
      dxmax = glmax(dxmax,1)
      dxmid = 0.5 * (dxmax + dxmin)
      if (dxmax-dxmin.gt.tol) then
         if (nid.eq.0)
     $      write(*,*)'box-Err elem distribusion is not uniform in ',s1
     $               ,dxmax,dxmin,dxmax-dxmin,tol
         ierr = 1
      endif

      dxref = (glmax(xx,nxyz*nelt) - glmin(xx,nxyz*nelt)) / nelx
      if (abs(dxmid-dxref).gt.tol) then
         if (nid.eq.0)
     $      write(*,*)'box-Err nel',s1,' and spacing mismatch'
     $               ,dxmid,dxref,nelx,abs(dxmid-dxref),tol
         ierr = 2
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine check_box2_length(nout, tol, ierr)
      implicit none
      include 'SIZE'
      include 'GEOM' ! xm1
      include 'PARALLEL' ! lglel
      include 'ZPER'
      real tol
      real dxmin, dymin, dzmin
     $   , dxmax, dymax, dzmax
     $   , dxlen, dylen, dzlen, dx, dy, dz
     $   , vlmin, vlmax, glmin, glmax
      integer nout, ierr, nt
     $      , e, eg, ex, ey, ez

      nt = lx1*ly1*lz1*nelt

      dxmin = -1e30
      dymin = -1e30
      dzmin = -1e30
      dxmax = 1e30
      dymax = 1e30
      dzmax = 1e30

      do e=1,nelt
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)

         dxlen = vlmax(xm1(1,1,1,e),nt) - vlmin(xm1(1,1,1,e),nt)
         dx = dxlen / nout
         if (ex.eq.1.OR.ex.eq.nelx) dx = dxlen / (nout-0.5)

         dylen = vlmax(ym1(1,1,1,e),nt) - vlmin(ym1(1,1,1,e),nt)
         dy = dylen / nout
         if (ey.eq.1.OR.ey.eq.nely) dy = dylen / (nout-0.5)

         dxmin = min(dxmin,dx)
         dxmax = max(dxmax,dx)
         dymin = min(dymin,dy)
         dymax = max(dymax,dy)
         
         if (ldim.eq.3) then
           if (ez.eq.1.OR.ez.eq.nelz) dz = dzlen / (nout-0.5)
           dzlen = vlmax(zm1(1,1,1,e),nt) - vlmin(zm1(1,1,1,e),nt)
           dz = dzlen / nout
           dzmin = min(dzmin,dz)
           dzmax = max(dzmax,dz)
         endif

      enddo
      dxmin = glmin(dxmin,1)
      dymin = glmin(dymin,1)
      dzmin = glmin(dzmin,1)
      dxmax = glmax(dxmax,1)
      dymax = glmax(dymax,1)
      dzmax = glmax(dzmax,1)

      if (dxmax-dxmin.lt.tol) ierr = 1
      if (dymax-dymin.lt.tol) ierr = 2
      if (dzmax-dzmin.lt.tol) ierr = 3

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_interp_matrix2(ni,no)
      implicit none
      include 'SIZE'
      include 'ZPER' ! nelx
      include 'PARALLEL' ! lglel

      integer ni, no, l, dd
      parameter (l=50)

      integer etype
      real jx, jxt, we, ze
      common /my_box_interp2/ jx(l*l,3),jxt(l*l,3),we(l*l),ze(l)
      common /my_box_interp2i/ etype(3,lelt)

      integer e,ex,ey,ez,eg,e3d(3),nelxyz(3),nt
      real dx, dx2, xlen, glmin, glmax

      if (no.gt.l) call exitti('box_interp2, memory$',no)

      nt = lx1*ly1*lz1*nelt

      ! set up matrices
      call zwgll(ze,we,ni)

      ! x-r direction
      call zuni_int_left(we,no)
      call gen_int_gz(jx(1,1),jxt(1,1),we,no,ze,ni)

      call zuni_int(we,no)
      call gen_int_gz(jx(1,2),jxt(1,2),we,no,ze,ni)

      call zuni_int_right(we,no)
      call gen_int_gz(jx(1,3),jxt(1,3),we,no,ze,ni)

      ! set up etype
      nelxyz(1) = nelx
      nelxyz(2) = nely
      nelxyz(3) = nelz
      do e=1,nelt
         eg = lglel(e)
         call get_exyz(e3d(1),e3d(2),e3d(3),eg,nelx,nely,nelz)

         do dd=1,ldim
            if (e3d(dd).eq.1) then
               etype(dd,e) = 1
            elseif (e3d(dd).eq.nelxyz(dd)) then
               etype(dd,e) = 3
            else
               etype(dd,e) = 2
            endif
         enddo

      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine box_outfld_rescale_mesh
c     Input: uniform box mesh sotred in xm1,ym1,zm1
c     Output: make boundary elements thinner for box_outfld (imode=2)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'RESTART'
      include 'ZPER'

      integer ierr, nout, e, eg, ex, ey, ez, nt
      real xlen, ylen, zlen, dx, dy, dz, glmin, glmax, x0, x1
     $   , exlen1, eylen1, ezlen1
     $   , exlen2, eylen2, ezlen2

      nt = lx1*ly1*lz1*nelt

      call check_unif_box_mesh(1, ierr)
      if (ierr.ne.0) then
         call exitti('mesh check failed$',ierr)
      else
         if (nio.eq.0) write(*,*)'box-chk uniform box mesh passed!'
      endif

      xlen = glmax(xm1,nt) - glmin(xm1,nt)
      ylen = glmax(ym1,nt) - glmin(ym1,nt)
      zlen = glmax(zm1,nt) - glmin(zm1,nt)

      nrg = lxo ! FIXME, nek5000 bug
      nout = nrg
      dx = xlen / (nelx * lxo - 1.0)
      dy = ylen / (nely * lxo - 1.0)
      dz = zlen / (nelz * lxo - 1.0)

      exlen1 = dx * nout
      exlen2 = dx * (nout-0.5)
      eylen1 = dy * nout
      eylen2 = dy * (nout-0.5)
      ezlen1 = dz * nout
      ezlen2 = dz * (nout-0.5)

      do e=1,nelt
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)

         x0 = 0.0
         if (ex.gt.1) x0 = exlen2
         if (ex.gt.2) x0 = x0 + (ex-2)*exlen1

         x1 = x0 + exlen1
         if (ex.eq.1.OR.ex.eq.nelx) x1 = x0 + exlen2
         call rescale_elem_x(xm1(1,1,1,e),x0,x1)

         x0 = 0.0
         if (ey.gt.1) x0 = eylen2
         if (ey.gt.2) x0 = x0 + (ey-2)*eylen1

         x1 = x0 + eylen1
         if (ey.eq.1.OR.ey.eq.nely) x1 = x0 + eylen2
         call rescale_elem_x(ym1(1,1,1,e),x0,x1)

         if (ldim.eq.3) then
            x0 = 0.0
            if (ez.gt.1) x0 = ezlen2
            if (ez.gt.2) x0 = x0 + (ez-2)*ezlen1

            x1 = x0 + ezlen1
            if (ex.eq.1.OR.ex.eq.nelx) x1 = x0 + ezlen2
            call rescale_elem_x(zm1(1,1,1,e),x0,x1)
         endif
      enddo

      call fix_geom
      call check_unif_box_mesh(2, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine rescale_elem_x(x,x0,x1)
      include 'SIZE'
      real x(1)

      n = lx1*ly1*lz1
      xmin = vlmin(x,n)
      xmax = vlmax(x,n)

      if (xmax.le.xmin) return

      scale = (x1-x0)/(xmax-xmin)
      do i=1,n
         x(i) = x0 + scale*(x(i)-xmin)
      enddo

      return
      end
