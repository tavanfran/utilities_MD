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


integer it,hkey,key,t,nat,j,k,i,num(natmax)
integer nn,index(natmax),n,numero(natmax)
integer fluct,h,step,step_end,dump
character*8 ctemp
character*2 atyp,chain(natmax),catena
character*2 atype(natmax)
character*3 resi(natmax)
character*70 ctemp2
character*80 title,tfilename

real*8 pos(3,natmax),vel(3,natmax),for(3,natmax)
real*8 temptot,raggio(natmax),speed,force
real*8 xa,xb,xc,ya,yb,yc,za,zb,zc,energy,itemp
real*8 fx(natmax),fy(natmax),fz(natmax),e_val(natmax)
real*8 max_e_val,min_e_val,diff_pcnt,thrld
real*8 xx,xy,xz,yx,yy,yz,zx,zy,zz
character*44 nothing
character*4 res(natmax),residuo
integer skip,nca,mass,charge,replica,index1,timestep
integer n_config,yn,t2000,t1500,t1000,t800,t650,t500,t300

character*1 option
character*4 atname,config_type
character*32 cell_size
character*76 steps
character*80 lammpsfile,pdbfile

!------------------------------------------------------------------
!     read the lammps file
!------------------------------------------------------------------
write(*,*)''
write(*,*)'*********************'
write(*,*)'  This script reads a CP2K trajectory file and writes the trajectory in xyz format'
write(*,*)'  The output is consistent with the GAP_QUIP and DeepMD packages'
write(*,*)'  written by Francesco Tavanti, Feb 2023'
write(*,*)'  CNR-NANO S3, Modena, Italy'
write(*,*)'*********************'
write(*,*)''

!write(*,*)'Write name of lammps trajcetory file'
!read(*,*)lammpsfile

open(10,file='./DOP-pos-1.xyz',status='old')
!open(16,file='./pos.xyz',status='old')
open(11,file='./DOP-frc-1.xyz',status='old')
open(15,file='./DOP-1.stress',status='old')
open(12,file='conf.xyz')
open(14,file='conf.pdb')

write(*,*)'Insert number of total frames'
read(*,*)step_end
write(*,*)'Dump frames'
read(*,*)dump
write(*,*)'Do you want forces?'
write(*,*)'1 = energy, forces, stress       (ALL)' 
write(*,*)'2 = energy, forces, no stress    (NO STRESS)'
write(*,*)'3 = energy, no forces, stress    (NO FORCES)'
write(*,*)'4 = energy, no forces, no stress (ENERGY ONLY)'
write(*,*)'STRESS available from frame 496 on'
read(*,*)yn
!step_end=237
write(*,*)'Insert threshold value for energies cuotff on configurations in % (e.g. 10 for 10%)'
read(*,*)thrld
thrld=thrld/100
read(15,*)
n_config=0
step=1
t2000=0
t1500=0
t1000=0
t800=0
t650=0
t500=0
t300=0


do it=1,step_end
    
    if(it.le.516)then
        config_type='2000'
    elseif((it.gt.516).and.(it.le.611))then
        config_type='1500'
    elseif((it.gt.611).and.(it.le.696))then
        config_type='1000'
    elseif((it.gt.696).and.(it.le.790))then
        config_type='800 '
    elseif((it.gt.790).and.(it.le.878))then
        config_type='650 '
    elseif((it.gt.878).and.(it.le.978))then
        config_type='500 '
    elseif(it.gt.978)then
        config_type='300 '
    endif

    read(10,*)nat
    !read(16,*)nat
    read(11,*)nat
    !write(*,*)nat
    
    read(10,'(a44,f20.10)')nothing,energy
    !read(16,*)
    read(11,'(a44,f20.10)')nothing,energy

    if(it.ge.496)then
        read(15,'(i8,f12.3,9f20.10)')i,itemp,xx,xy,xz,yx,yy,yz,zx,zy,zz
        xx=xx*6.2417E-07
        xy=xy*6.2417E-07
        xz=xz*6.2417E-07
        yx=yx*6.2417E-07
        yy=yy*6.2417E-07
        yz=yz*6.2417E-07
        zx=zx*6.2417E-07
        zy=zy*6.2417E-07
        zz=zz*6.2417E-07
    else

  
    endif
    xb=20.31669
    yb=20.31669
    zb=20.31669

    e_val(it)=energy
    write(14,'(a6,3f9.3,3f7.2,a16)')'CRYST1',xb,yb,zb,90.0,90.0,90.0,' P 1           1'
    write(14,'(a5,i9)')'MODEL',it
    

            
    do j=1,nat
        read(11,*)atype(j),fx(j),fy(j),fz(j)
        fx(j)=fx(j)*51.42208619083232
        fy(j)=fy(j)*51.42208619083232
        fz(j)=fz(j)*51.42208619083232

        read(10,*)atype(j),(pos(k,j),k=1,3)
        !read(16,*)atype(j),(pos(k,j),k=1,3)
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
        

        write(14,1070)'ATOM',j,atype(j),' A ',j,(pos(k,j),k=1,3),&
&        0.0,0.0
        
 
        
        
    enddo
    if(it.gt.1)then
        max_e_val=max(abs(e_val(it)),abs(e_val(it-1)))
        min_e_val=min(abs(e_val(it)),abs(e_val(it-1)))
        diff_pcnt=(max_e_val-min_e_val)*100/max_e_val
        if(diff_pcnt.ge.thrld)then
            
            energy=energy*27.2114
            if((it.ge.step).and.(it.le.step_end).and.((modulo(it,dump)).eq.0))then
                n_config=n_config+1
            write(12,'(i6)')nat
            if(it.ge.496)then
                
                if(yn.eq.1)then
                write(12,'(a28,9f10.5,a1,a51,9f20.16,a9,f11.4,a12)')'config_type=liqu Lattice=" ', &
                xb,0.0,0.0,0.0,yb,0.0,0.0,0.0,zb,'"',' Properties=species:S:1:pos:R:3:forces:R:3 stress="', &
                xx,xy,xz,yx,yy,yz,zx,zy,zz,'" energy=',energy,' pbc="T T T"'
                elseif(yn.eq.2)then
                    write(12,'(a28,9f10.5,a1,a51,f11.4,a12)')'config_type=liqu Lattice=" ', &
                    xb,0.0,0.0,0.0,yb,0.0,0.0,0.0,zb,'"',' Properties=species:S:1:pos:R:3:forces:R:3 energy=', &
                    energy,' pbc="T T T"'
                elseif(yn.eq.3)then
                write(12,'(a28,9f10.5,a1,a40,9f20.16,a9,f11.4,a12)')'config_type=liqu Lattice=" ', &
                xb,0.0,0.0,0.0,yb,0.0,0.0,0.0,zb,'"',' Properties=species:S:1:pos:R:3 stress="', &
                xx,xy,xz,yx,yy,yz,zx,zy,zz,'" energy=',energy,' pbc="T T T"'
                elseif(yn.eq.4)then
                    write(12,'(a28,9f10.5,a1,a40,f11.4,a12)')'config_type=liqu Lattice=" ', &
                    xb,0.0,0.0,0.0,yb,0.0,0.0,0.0,zb,'"',' Properties=species:S:1:pos:R:3 energy=', &
                    energy,' pbc="T T T"'
                endif
            do j=1,nat
                if((yn.eq.1).or.(yn.eq.2))then
                    write(12,'(a2,1x,6f10.5)')atype(j),pos(1,j),pos(2,j),pos(3,j),fx(j),fy(j),fz(j)
                else
                    write(12,'(a2,1x,3f10.5)')atype(j),pos(1,j),pos(2,j),pos(3,j)
                endif
            enddo
            if(config_type.eq.'2000')then
                t2000=t2000+1
            elseif(config_type.eq.'1500')then
                t1500=t1500+1
            elseif(config_type.eq.'1000')then
                t1000=t1000+1
            elseif(config_type.eq.'800 ')then
                t800=t800+1
            elseif(config_type.eq.'650 ')then
                t650=t650+1
            elseif(config_type.eq.'500 ')then
                t500=t500+1
            else
                t300=t300+1
            endif
            else
            if((yn.eq.1).or.(yn.eq.2))then
                write(12,'(a28,9f10.5,a1,a51,f11.4,a12)')'config_type=liqu Lattice=" ', &
                xb,0.0,0.0,0.0,yb,0.0,0.0,0.0,zb,'"',' Properties=species:S:1:pos:R:3:forces:R:3 energy=', &
                energy,' pbc="T T T"'
            else
                write(12,'(a28,9f10.5,a1,a40,f11.4,a12)')'config_type=liqu Lattice=" ', &
                xb,0.0,0.0,0.0,yb,0.0,0.0,0.0,zb,'"',' Properties=species:S:1:pos:R:3 energy=', &
                energy,' pbc="T T T"'
            endif
            do j=1,nat
                if((yn.eq.1).or.(yn.eq.2))then
                    write(12,'(a2,1x,6f10.5)')atype(j),pos(1,j),pos(2,j),pos(3,j),fx(j),fy(j),fz(j)
                else
                    write(12,'(a2,1x,3f10.5)')atype(j),pos(1,j),pos(2,j),pos(3,j)
                endif
            enddo
            if(config_type.eq.'2000')then
                t2000=t2000+1
            elseif(config_type.eq.'1500')then
                t1500=t1500+1
            elseif(config_type.eq.'1000')then
                t1000=t1000+1
            elseif(config_type.eq.'800 ')then
                t800=t800+1
            elseif(config_type.eq.'650 ')then
                t650=t650+1
            elseif(config_type.eq.'500 ')then
                t500=t500+1
            else
                t300=t300+1
            endif
            endif
        endif
    endif
    else
        
    endif

    write(14,170)'ENDMDL'
enddo
write(*,'(a32,f10.4,a19,i5)')'Number of configurations within ',thrld*100,' % of energy are ',n_config

write(*,*)'Number of configs at each temperature:'
write(*,*)'2000K = ',t2000,float(t2000)/float(n_config)*100,'%'
write(*,*)'1500K = ',t1500,float(t1500)/float(n_config)*100,'%'
write(*,*)'1000K = ',t1000,float(t1000)/float(n_config)*100,'%'
write(*,*)'800K =  ',t800,float(t800)/float(n_config)*100,'%'
write(*,*)'650K =  ',t650,float(t650)/float(n_config)*100,'%'
write(*,*)'500K =  ',t500,float(t500)/float(n_config)*100,'%'
write(*,*)'300K =  ',t300,float(t300)/float(n_config)*100,'%'

170  format(a6)

900  format(a80)

1070 format(a4,2x,i5,1x,a2,3x,a3,i6,4x,f8.3,f8.3,f8.3,f6.2,f6.2)

end
!------------------------------------------------------------------------

