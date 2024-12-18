      program rdfcode 
c AGPL-v3.0 © 2024 [Elizabeth A. Ploetz and Paul E. Smith]

c This program is free software: you can redistribute it and/or modify
c it under the terms of the GNU Affero General Public License as published
c by the Free Software Foundation, either version 3 of the License, or
c (at your option) any later version.

c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU Affero General Public License for more details.

c You should have received a copy of the GNU Affero General Public License
c along with this program (see LICENSE file).

      implicit real*8 (a-h,o-z)
      parameter(nmax=100000,nbin=500,dist=5.00,ncnfgSkip=0)
      parameter(ngroup=4,maxatm=50,ntail=2)
c
      dimension x(nmax),y(nmax),z(nmax)
      dimension rrdf(ngroup,ngroup,nbin) !rrdf=counting as read traj. 
      dimension rdf(ngroup,ngroup)       !rdf=outputted rdf
      dimension w(ngroup,maxatm),ns(ngroup),ms(ngroup)
      dimension wtot(ngroup)
      dimension cx(maxatm),cy(maxatm),cz(maxatm)
      dimension ro(ngroup,ngroup),cno(ngroup,ngroup)   !ro=summing rdf over bins,
      dimension rnorm(ngroup,ngroup),rn(ngroup,ngroup) !rn counts #i-j pairs
      character ifile*50,funn*80,funn2*3
c 
      open(unit=10,file='rdf.dat',status='old') 
c
      read(10,*) nh
      ntot=0
      natm=0
      do i=1,ngroup
        read(10,*) funn2
        read(10,*) ns(i)
        ntot=ntot+ns(i)
        read(10,*) ms(i)
        if(ms(i).gt.maxatm) stop ' Too many atoms in mol'
        wtot(i)=0.0
        do j=1,ms(i)
          natm=natm+ns(i)
          read(10,*) w(i,j)
          wtot(i)=wtot(i)+w(i,j)
        end do
        read(10,*) funn 
        write(6,*)funn2,': ',i,ns(i),ms(i)
      end do
      if(ntot.gt.nmax) stop ' Too many mols'
      read(10,*) ntf
      read(10,*) nfiles
c
      dr    = dist/nbin
      dri   = 1.d00/dr
      pi    = acos(-1.d00)
      avol  = 0.0
      ncnfg = 0
c
      do irow=1,ngroup
        do icol=1,ngroup
          rdf(irow,icol)=0.0
          ro(irow,icol)=0.0
          cno(irow,icol)=0.0
          rnorm(irow,icol)=0.0
          do i=1,nbin
            rrdf(irow,icol,i)=0.0
          end do
        end do
      end do
c
      do ifiles=1,nfiles
        read(10,'(a)') ifile
        open(unit=20,file=ifile,status='old')
        do i=1,ncnfgSkip*(natm+nh+ntail+2)
          read(20,'(a80)') funn
        end do 
        do it=1,ntf
          do i=1,nh 
            read(20,'(a80)') funn
            write(6,'(a80)') funn
          end do
          read(20,'(6x,3f9.3)') xbox,ybox,zbox
          xbox=xbox*0.1
          ybox=ybox*0.1
          zbox=zbox*0.1
          write(6,'(6x,3f9.3)') xbox,ybox,zbox
          read(20,'(a80)') funn
          write(6,'(a80)') funn
          avol=avol+xbox*ybox*zbox
c
          imol=1
          do i=1,ngroup
c           if(it.eq.1)write(6,*)'read coord:',i,imol,ns(i)+imol-1
            do ii=imol,ns(i)+imol-1
              comx=0.0
              comy=0.0
              comz=0.0
              do j=1,ms(i)
                read(20,'(30x,3f8.3)') cx(j),cy(j),cz(j)
                comx=comx+w(i,j)*cx(j)
                comy=comy+w(i,j)*cy(j)
                comz=comz+w(i,j)*cz(j)
              end do
              comx=comx/wtot(i)
              comy=comy/wtot(i)
              comz=comz/wtot(i)
              x(ii)=comx*0.1
              y(ii)=comy*0.1
              z(ii)=comz*0.1
            end do
            imol=imol+ns(i)
          end do
          ncnfg = ncnfg + 1
c
          istart=1
          do irow=1,ngroup
            jstart=1
            do icol=1,ngroup
c             if(it.eq.1)write(6,*)'irow,icol',irow,icol
              rn(irow,icol)=0.0
              if(irow.eq.icol)then
                iend=istart+ns(irow)-2
              else
                iend=istart+ns(irow)-1
              end if
              do i=istart,iend
c               if(it.eq.1)write(6,*)'i loop:',istart,iend,i
                xi=x(i)
                yi=y(i)
                zi=z(i)
                if(irow.eq.icol)then
                  jstart2=i+1
                  jend=istart+ns(icol)-1
                else
                  jstart2=jstart
                  jend=jstart+ns(icol)-1
                end if
                do j=jstart2,jend
c                 if(it.eq.1)write(6,*)'j loop:',jstart2,jend,j
c                 if(it.eq.1)write(60,*)'index:',irow,icol,i,j
                  rn(irow,icol)=rn(irow,icol)+1.0
                  xij=x(j)-xi
                  yij=y(j)-yi
                  zij=z(j)-zi
                  xij=xij-nint(xij/xbox)*xbox
                  yij=yij-nint(yij/ybox)*ybox
                  zij=zij-nint(zij/zbox)*zbox 
                  rr=sqrt(xij*xij+yij*yij+zij*zij)
c                 if(rr.le.0.01)write(6,*)rr,irow,icol,i,j
                  if(rr.le.dist) then
                    ibn=rr*dri+1
                    rrdf(irow,icol,ibn)=rrdf(irow,icol,ibn)+1
                  endif
                end do
              end do
              jstart=jstart+ns(icol)
            end do
            istart=istart+ns(irow)
          end do
          do i=1,ntail
            read(20,'(a80)') funn 
          end do 
        end do !over a file
        close(unit=20)
      end do !over files
c 
      avol=avol/ncnfg
      write(6,*) 'NCONFIG = ',ncnfg
      do i=1,ngroup 
        write(30+i,'(a,f12.6)') '#',avol 
        write(60+i,'(a,f12.6)') '#',avol 
      end do
      do i=1,ngroup 
        do j=1,ngroup
          ro(i,j)=0.0
        end do
      end do
c
      do i=1,nbin
        r=(i-0.5)*dr
        rlow=(i-1)*dr
        rhigh=i*dr
        dv=4.0*pi*(rhigh**3-rlow**3)/3.0
        do irow=1,ngroup 
          do icol=1,ngroup
            ro(irow,icol)=ro(irow,icol)+rrdf(irow,icol,i)
            rnorm(irow,icol)=dv*rn(irow,icol)/avol
            rdf(irow,icol)=rrdf(irow,icol,i)/rnorm(irow,icol)/ncnfg
            cno(irow,icol)=real(ro(irow,icol))/ncnfg/ns(irow)
            if(irow.eq.icol)cno(irow,icol)=2*cno(irow,icol)
          end do
          if(i.eq.1)write(30,*)irow,(rn(irow,icol),icol=1,ngroup)
          write(30+irow,*)r,(rdf(irow,icol),icol=1,ngroup)
          write(60+irow,*)r,(cno(irow,icol),icol=1,ngroup)
        end do
      end do
c         
      stop
      end                
