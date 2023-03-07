program hisread
!---------------------------------------------------------------------
!   read lammps trajectory file and wirte the trajectory in PDB format
!   written by Francesco Tavanti, March 2020
!   CNR-NANO S3, Modena, Italy
!
!   USAGE:
!   build the fortran executable ==> gfortran lammpstraj2pdb.f90
!   then ./a.out
!   insert the name of input lammps file
!   insert the name of output file in pdb format
!---------------------------------------------------------------------
implicit none
integer natmax,ncouplemax

parameter(natmax=50000)


integer it,hkey,pdbkey,t,nat,itemp,j,k,i,num(natmax)
integer nn,index(natmax),n,numero(natmax)
integer fluct,h,step,step_end,dump,eneyn,foryn
character*8 ctemp
character*2 atyp,chain(natmax),catena
character*2 atype(natmax)
character*3 resi(natmax)
character*70 ctemp2
character*80 title,tfilename

real*8 pos(3,natmax),vel(3,natmax),for(3,natmax)
real*8 temptot,raggio(natmax),speed,force
real*8 xa,xb,xc,ya,yb,yc,za,zb,zc,energy
real*8 fx(natmax),fy(natmax),fz(natmax),nothing
character*8 chskip
character*4 res(natmax),residuo
integer skip,nca,mass,charge,replica,index1,timestep,yn

character*1 option
character*4 atname
character*80 lammpsfile,pdbfile

!------------------------------------------------------------------
!     read the lammps file
!------------------------------------------------------------------
write(*,*)''
write(*,*)'*********************'
write(*,*)'  This script reads a lammps trajectory file and writes the trajectory in PDB format'
write(*,*)'  The trajectory could come both from serial and from parallel simulations'
write(*,*)'  written by Francesco Tavanti, March 2020'
write(*,*)'  CNR-NANO S3, Modena, Italy'
write(*,*)'*********************'
write(*,*)''

write(*,*)'Write name of lammps trajcetory file'
!read(*,*)lammpsfile

open(10,file='Ge50Se50_P.pos',status='old')
open(11,file='Ge50Se50_P.evp')

open(13,file='GeSe.label')
open(14,file='Ge50Se50_P.for')
open(15,file='conf_deepmd.xyz')

!write(*,*)'Do you want energy? 0 = NO , 1 = YES'
!read(*,*)eneyn

eneyn=1
!write(*,*)'Do you want forces? 0 = NO , 1 = YES'
!read(*,*)foryn
foryn=1
write(*,*)'Dump frames'
read(*,*)dump
write(*,*)'Do you want forces? 0=no 1=yes'
read(*,*)yn

do j=1,50
    read(13,'(a2)')atype(j)
enddo

!dump=10

step=1
step_end=100000
read(11,*)

xb= 21.4735
yb= 21.4735
zb= 21.4735
do it=1,step_end
    read(10,900)title
    read(11,*)nothing,nothing,nothing,nothing,nothing,energy
if((it.ge.step).and.(it.le.step_end).and.((modulo(it,dump)).eq.0))then

    write(15,'(i2)')50
    if(yn.eq.1)then
    write(15,'(a27,9f10.5,a1,a52,f11.4,a12)')'config_type=1000 Lattice=" ',xb/1.8897,0.0,0.0,0.0,yb/1.8897,0.0,0.0,0.0,zb/1.8897, &
    '"',' Properties=species:S:1:pos:R:3:forces:R:3 energy=',energy*27.2114,' pbc="T T T"'
    else
    write(15,'(a27,9f10.5,a1,a42,f11.4,a12)')'config_type=1000 Lattice=" ',xb/1.8897,0.0,0.0,0.0,yb/1.8897,0.0,0.0,0.0,zb/1.8897, &
        '"',' Properties=species:S:1:pos:R:3 energy=',energy*27.2114,' pbc="T T T"'
    endif
endif


read(14,*)title
do j=1,50
read(14,*)fx(j),fy(j),fz(j)
fx(j)=fx(j)*51.422086
fy(j)=fy(j)*51.422086
fz(j)=fz(j)*51.422086
if(foryn.eq.0)then
fx(j)=0
fy(j)=0
fz(j)=0
endif
read(10,*)(pos(k,j),k=1,3)
if((pos(1,j).gt.xb))then
    pos(1,j)=MOD(pos(1,j),xb)
elseif((pos(1,j).lt.0))then
    pos(1,j)=MOD(pos(1,j),xb)+xb
endif
if((pos(2,j).gt.yb))then
    pos(2,j)=MOD(pos(2,j),yb)
elseif((pos(2,j).lt.0))then
    pos(2,j)=MOD(pos(2,j),yb)+yb
endif
if((pos(3,j).gt.zb))then
    pos(3,j)=MOD(pos(3,j),zb)
elseif((pos(3,j).lt.0))then
    pos(3,j)=MOD(pos(3,j),zb)+zb
endif




!if(pos(1,j).gt.xb)then
!    pos(1,j)=pos(1,j)-xb
!elseif(pos(1,j).lt.0)then
!    pos(1,j)=pos(1,j)+xb
!elseif(pos(2,j).gt.yb)then
!    pos(2,j)=pos(2,j)-yb
!elseif(pos(2,j).lt.0)then
!    pos(2,j)=pos(2,j)+yb
!elseif(pos(3,j).gt.yb)then
!    pos(3,j)=pos(3,j)-zb
!elseif(pos(3,j).lt.0)then
!    pos(3,j)=pos(3,j)+zb
!else
!endif



if((it.ge.step).and.(it.le.step_end).and.((modulo(it,dump)).eq.0))then
    if(yn.eq.1)then
    write(15,'(a2,6f17.9)')atype(j),pos(1,j)/1.8897,pos(2,j)/1.8897,pos(3,j)/1.8897,fx(j),fy(j),fz(j)
    else
    write(15,'(a2,3f17.9)')atype(j),pos(1,j)/1.8897,pos(2,j)/1.8897,pos(3,j)/1.8897
    endif
endif
enddo


enddo

170  format(a6)

900  format(a80)

1070 format(a4,2x,i5,1x,a2,3x,a3,i6,4x,f8.3,f8.3,f8.3,f6.2,f6.2)

end
!------------------------------------------------------------------------

