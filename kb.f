      program kb
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

      implicit real*8(a-h,o-z)
      parameter(ngroup=4,ndat=500) 
      parameter(fudge=0.0051)
      dimension ns(ngroup),rblk(ngroup)    !rblk=rho_bulk=total system rho
      dimension cn0(ngroup,ngroup),G(ngroup,ngroup)
      dimension cn(ngroup,ngroup)
      dimension aG(ngroup,ngroup),Gr(ngroup,ngroup)
      dimension Gu(ngroup,ngroup),rbth(ngroup,ngroup) !rbth=rho_bath=rho of bath beyond sphere 
      dimension aN(ngroup,ngroup)
      character ofile*60,title*80,c1*1
      character funn*80,funn2*3
      character string*32
c
      open(unit=9,file='kb.dat',status='old')
c 
      read(9,'(a80)') title
      do i=1,ngroup
        do j=1,ngroup
          cn0(i,j)=0.0
          G(i,j)=0.0
          aG(i,j)=0.0
          Gr(i,j)=0.0
          Gu(i,j)=0.0
          rbth(i,j)=0.0
          aN(i,j)=0.0
        end do 
      end do 
      do i=1,ngroup
        read(9,*) funn2
        read(9,*) ns(i)
        read(9,*) funn 
        write(6,*)funn2,': ',i,ns(i)
      end do
      read(9,*) dr
      read(9,*) dist
      read(9,*) rl,ru 
      read(9,'(a)') ofile
c
      do i=1,ngroup
        read(60+i,*) c1,avol
        write(*,'(1x,f12.6)') avol 
        rblk(i)=ns(i)/avol !molecules/nm3
      end do
      open(unit=26,file=ofile,status='unknown')
c
      cons  = 6.022d+2 !(molec/nm3)/cons  = mole/cm3  
      cons2 = 0.6022   !(molec/nm3)/cons2 = mole/L = molarity
      pi    = dacos(-1.d00)
      cnst  = 4.0*pi/3.0
c
      count = 0.0
      do i = 1,ndat
        do ii=1,ngroup
          read(60+ii,*) r,(cn(ii,jj),jj=1,ngroup)
        end do
        dvol=avol-cnst*(r+0.5*dr)**3 !volume outside of sphere
        do ii=1,ngroup
          do jj=1,ngroup
            if(ii.eq.jj)rbth(ii,jj)=(ns(jj)-cn(ii,jj)-1)/dvol !bath concentration of j around i
            if(ii.ne.jj)rbth(ii,jj)=(ns(jj)-cn(ii,jj))/dvol   !bath concentration of j around i
          end do
        end do
        write(string,'("(a5,",I2,"(4x,I1,a1,I1,2x))")')ngroup*ngroup
        if(i.eq.1)write(*,string) 'r    ',((ii,'-',jj,jj=1,ngroup),
     *                                                 ii=1,ngroup)
        write(string,'("(1f7.3,",I2,"(f9.5))")')ngroup*ngroup
        write(*,string) r,((rbth(ii,jj)/cons2,jj=1,ngroup),ii=1,ngroup)
        if( r .le. dist )then
          do ii=1,ngroup
            do jj=1,ngroup
              Gr(ii,jj)=(cn(ii,jj)-cn0(ii,jj))/rbth(ii,jj)
     *                            -cnst*((r+0.5*dr)**3-(r-0.5*dr)**3)
              cn0(ii,jj)=cn(ii,jj)
              G(ii,jj)=G(ii,jj)+Gr(ii,jj)
              Gu(ii,jj)=G(ii,jj)*cons
              if(r.ge.(rl-fudge).and.r.le.(ru+fudge))then
                if(ii.eq.1.and.jj.eq.1)count  = count + 1.0
                aG(ii,jj)=aG(ii,jj)+Gu(ii,jj)
              end if
            end do
          end do
          do ii=1,ngroup
            write(40+ii,*) r,(Gu(ii,jj),jj=1,ngroup)
          end do
        end if
        do ii=1,ngroup
          do jj=1,ngroup
            Gu(ii,jj)=Gu(ii,jj)/1000.0 !convert cm3/mol to L/mol
          end do
        end do

        write(52,*)r,Gu(1,1),Gu(1,2),Gu(2,2),Gu(3,3),
     *               Gu(3,4),Gu(4,4),Gu(1,3),Gu(1,4),
     *               Gu(2,3),Gu(2,4)
      end do
      cnti = 1.0/count
      do i=1,ngroup
        rblk(i)=rblk(i)/cons2
        do j=1,ngroup
          rbth(i,j)=rbth(i,j)/cons2
          aG(i,j)=aG(i,j)*cnti
        end do
      end do
      do i=1,ngroup 
        do j=1,ngroup
c         aN(i,j)=rblk(j)*aG(i,j)*0.001   !excess coordination number, Nij, usual way
          aN(i,j)=rbth(i,j)*aG(i,j)*0.001 !using bath concentration
c         write(6,*)i,j,rblk(j),rbth(i,j)
        end do
      end do
      write(26,*) 'N:',count
      write(26,*) 'KBIs(cm3/mol): 1 row for each i, 1 column for each j'
      do i=1,ngroup
        write(26,*) i,(aG(i,j),j=1,ngroup)
      end do
      write(26,*) ''
      write(26,*) 'Nij: 1 row for each i, 1 column for each j'
      do i=1,ngroup
        write(26,*) i,(aN(i,j),j=1,ngroup)
      end do
      write(*,*) (i,rblk(i),i=1,ngroup)
c     stop
      end
