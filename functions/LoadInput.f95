module loadInput
	use general
	use configuration
	implicit none
	
	! GLOBAL variables
	real*8 switch,Eshift
	real*8,allocatable:: EmatU(:),EmatD(:)
	complex*16,allocatable:: VmatUU(:,:,:,:),VmatDD(:,:,:,:),VmatUD(:,:,:,:)
	complex*16,allocatable:: VmatUU2(:,:,:,:),VmatDD2(:,:,:,:),VmatUDU2(:,:,:),VmatDUD2(:,:,:)
	complex*16, allocatable:: SelfEneU(:,:),SelfEneD(:,:),VertU(:,:),VertD(:,:)
	logical ifLowMem
	
	
contains

subroutine loadInputs(filEU,filED,filVUU,filVDD,filVUD,pathfile,screenFolder,nVUU,nVDD,nVUD,nVU,nVD,screen)
	! loads basic inputs: SP energies and CME
	! also scales the inputs by required parameters
	implicit none
	! INPUTS: - names of files and folders (character)
	!		  - nVUU,nVDD,nVUD - numbers of CME to read
	!		  - nVU,nVD - numbers of SP states considered for CME - determines the size of CME tensor
	!		  - screen - dielectric constant
	character*8,intent(in):: filEU,filED,filVUU,filVUD,filVDD,screenFolder
	character*9,intent(in)::pathfile
	integer,intent(in):: nVUU,nVDD,nVUD,nVU,nVD
	real*8,intent(in):: screen
	
	integer i,x,y,z,w,dum
	real*8 Vr,Vi,E
	
	
	!!! READ SP ENERGIES
	allocate(EmatU(nEU),EmatD(nED))
	EmatU=0.d0; EmatD=0.d0;
	open(unit=21,file=pathfile//filEU,action='read')
	open(unit=22,file=pathfile//filED,action='read')
	
	do i=1,nEU
		read(21,*) E
		EmatU(i)=E		
	end do
	EmatU=EmatU-Eshift
	close(21)
	
	do i=1,nED
		read(22,*) E
		EmatD(i)=E		
	end do
	close(22)
	EmatD=EmatD-Eshift
	
	! if many SP energies are considered, low memory convention is assumed later in self energy calculation
	if (nEU.lt.130) then
		ifLowMem=.false.
		
		!!! READ CME for 2 electrons with spins UU, DD or UD
		allocate(VmatUU(nVU,nVU,nVU,nVU),VmatDD(nVD,nVD,nVD,nVD),VmatUD(nVU,nVD,nVD,nVU))
		VmatUU=(0.d0,0.d0); VmatUD=(0.d0,0.d0); VmatDD=(0.d0,0.d0)
		 
		open(unit=23,file=pathfile//screenFolder//filVUU,action='read')
		do i=1,nVUU
			read(23,*) x,y,z,w,Vr,Vi
			! scaling by dielectric constant and strength of Coulomb interaction
			VmatUU(x,y,z,w)=(Vr+im*Vi)/screen*switch 
			VmatUU(y,x,w,z)=(Vr+im*Vi)/screen*switch
		end do
		close(23)
		print *, VmatUU(1,2,2,1)
		
		open(unit=24,file=pathfile//screenFolder//filVDD,action='read')
		do i=1,nVDD
			read(24,*) x,y,z,w,Vr,Vi
			! scaling by dielectric constant and strength of Coulomb interaction
			VmatDD(x,y,z,w)=(Vr+im*Vi)/screen*switch
			VmatDD(y,x,w,z)=(Vr+im*Vi)/screen*switch
		end do
		close(24)
		
		open(unit=25,file=pathfile//screenFolder//filVUD,action='read')
		do i=1,nVUD
			read(25,*) x,y,z,w,Vr,Vi
			! scaling by dielectric constant and strength of Coulomb interaction
			VmatUD(x,y,z,w)=(Vr+im*Vi)/screen*switch
		end do
		close(25)
	else
	
		ifLowMem=.true.
		
	end if
	
end subroutine loadInputs

subroutine loadSelfEne1X(filEU,filED,pathfile,screenFolder,nVU,nVD,screen)
	! loads precalculated self energies for 1 exciton caculation
	! also scales the inputs by required parameters
	implicit none
	! INPUTS: - names of files and folders (character)
	!		  - number of elements to read
	!		  - screen - dielectric constant
	character*8,intent(in):: filEU,filED,screenFolder
	character*9,intent(in)::pathfile
	integer,intent(in):: nVU,nVD
	real*8,intent(in):: screen
	
	integer i,x,y,z,w,dum
	real*8 Vr,Vi,E
	
	print *, "nEU,nED:",nEU,nED
	! for many SP energies a low memory convention needs to be used
	if (nEU.lt.130) then
	
		allocate(SelfEneU(nEU,nEU),SelfEneD(nED,nED))
		SelfEneU=(0.d0,0.d0); SelfEneD=(0.d0,0.d0)
		open(unit=21,file=pathfile//screenFolder//filEU,action='read')
		open(unit=22,file=pathfile//screenFolder//filED,action='read')
		
		print *, "nVU",nVU
		do i=1,nVU
			read(21,*) x,y,Vr,Vi
			SelfEneU(x,y)=(Vr+im*Vi)/screen*switch
		end do
		close(21)
		
		do i=1,nVD
			read(22,*) x,y,Vr,Vi
			SelfEneD(x,y)=(Vr+im*Vi)/screen*switch
		end do
		close(22)
		
	else
	
		ifLowMem=.true.
	
	end if
	
end subroutine loadSelfEne1X

subroutine changeInput(fil,fil2,fil3,fil4,screenfolder,pathfile,nline,udborder)
	! transposes input CME from CME calculation format (1 file) 
	! to CI calculation format (3 files: UU, DD, UD)
	implicit none
	! INPUTS: - names of files and folders (character)
	!		  - nline - number of elements to read
	!		  - udborder - last number of the first spin, i.e. U:1-12, D:13-24, udborder=12
	character*8,intent(in):: fil,fil2,fil3,fil4,screenfolder
	character*9,intent(in):: pathfile
	integer,intent(in):: nline,udborder
	
	integer i,line(4),inds(4),ii,inds2(2),uu(2),dd(2),indu(4),indd(4)
	real*8 Vr,Vi,E,fac
	logical ifud(4)
	
	fac=1.d0
	open(unit=21,file=pathfile//screenfolder//fil,action='read')
	open(unit=22,file=pathfile//screenfolder//fil2,action='write',status='replace')
	open(unit=23,file=pathfile//screenfolder//fil3,action='write',status='replace')
	open(unit=24,file=pathfile//screenfolder//fil4,action='write',status='replace')
	
	
	print *, "udborder",udborder
	
	do i=1,nline
		read(21,*) line,Vr,Vi		
		ifud=line>udborder
		if (any(ifud)) then
			if (any(.not.ifud)) then
				if ((ifud(1).eqv.ifud(4)).and.(ifud(2).eqv.ifud(3))) then
					inds=0
					inds=pack([(ii,ii=1,4)],ifud)
					line(pack(inds,inds>0))=line(pack(inds,inds>0))-udborder
					if (ifud(1)) write(22,*) line,Vr/fac,Vi/fac	!UD					
				end if
				
			else
				write(23,*) line-udborder,Vr/fac,Vi/fac	!U
			end if
		else
			write(24,*) line,Vr/fac,Vi/fac !D	
			
		end if
	end do
	
	
	close(21)
	
	
end subroutine changeInput

subroutine loadQNumbers(Lfile,Vfile,pathfile)
	! loads angular momentum and valley quantum numbers
	implicit none
	! INPUTS: input files and folder
	character*8,intent(in):: Lfile,Vfile
	character*9,intent(in)::pathfile
	
	allocate(vecL(nEU),vecV(nEU))
	
	open(unit=121,file=pathfile//Lfile,action='read')
	open(unit=122,file=pathfile//Vfile,action='read')
	
	if (nEU.gt.100) then
		vecL=0
		vecV=0
	else
		read(121,*)vecL
		read(122,*)vecV
	end if
	
	close(121)
	close(122)


end subroutine loadQNumbers

subroutine loadSelfInputs(filVUU,filVDD,filVUD,pathfile,screenFolder,nVUU,nVDD,nVUD,nVU,nVD,screen,subfold)
	! loads basic inputs for a separate self energy caculation: SP energies and CME
	! also scales the inputs by required parameters
	implicit none
	! INPUTS: - names of files and folders (character)
	!		  - nVUU,nVDD,nVUD - numbers of CME to read
	!		  - nVU,nVD - numbers of SP states considered for CME - determines the size of CME tensor
	!		  - screen - dielectric constant
	character*8,intent(in):: filVUU,filVUD,filVDD,screenFolder
	character,intent(in)::pathfile*9,subfold*7
	integer,intent(in):: nVUU,nVDD,nVUD,nVU,nVD
	real*8,intent(in):: screen
	
	integer i,x,y,z,w,dum,ind3(3),exdir,uord
	real*8 Vr,Vi,E
	character folder*8
	integer(kind=8) iis(2)	
	
	folder=subfold//"/"
	open(unit=23,file=pathfile//screenFolder//folder//filVUU,action='read')
	open(unit=24,file=pathfile//screenFolder//folder//filVDD,action='read')
	open(unit=25,file=pathfile//screenFolder//folder//filVUD,action='read')
	
	! for many SP energies and size of CME tensor, it is better to index 
	! all CME needed for self energy calcualtion with 3 indices,
	! not 4, like for a general CI calculation
	if ((nVU.gt.130).or.(nVD.gt.130)) then
		
		ifLowMem=.true.
		allocate(VmatUU2(nVU,nVU,nVU,2),VmatDD2(nVD,nVD,nVD,2),VmatUDU2(nVU,nVD,nVU),VmatDUD2(nVD,nVU,nVD))
		VmatUU2=(0.d0,0.d0); VmatUDU2=(0.d0,0.d0); VmatDD2=(0.d0,0.d0); VmatDUD2=(0.d0,0.d0)
		
		do i=1,nVUU
			read(23,*) x,y,z,w,Vr,Vi
			call whichXYZWsameInd3([x,y,z,w],ind3,exdir,uord)
			VmatUU2(ind3(1),ind3(2),ind3(3),exdir+1)=(Vr+im*Vi)/screen*switch
		end do
		close(23)
		
		do i=1,nVDD
			read(24,*) x,y,z,w,Vr,Vi
			call whichXYZWsameInd3([x,y,z,w],ind3,exdir,uord)
			VmatDD2(ind3(1),ind3(2),ind3(3),exdir+1)=(Vr+im*Vi)/screen*switch
		end do
		close(24)
		
		do i=1,nVUD
			read(25,*) x,y,z,w,Vr,Vi
			call whichXYZWsameInd3([x,y,z,w],ind3,exdir,uord)
			if (uord.eq.1) then
				VmatUDU2(ind3(1),ind3(2),ind3(3))=(Vr+im*Vi)/screen*switch
			else
				VmatDUD2(ind3(1),ind3(2),ind3(3))=(Vr+im*Vi)/screen*switch
			end if
		end do
		close(25)
		
	
	else
		ifLowMem=.false.
		
		allocate(VmatUU(nVU,nVU,nVU,nVU),VmatDD(nVD,nVD,nVD,nVD),VmatUD(nVU,nVD,nVD,nVU))
		VmatUU=(0.d0,0.d0); VmatUD=(0.d0,0.d0); VmatDD=(0.d0,0.d0)
	
	
		do i=1,nVUU
			read(23,*) x,y,z,w,Vr,Vi
			VmatUU(x,y,z,w)=(Vr+im*Vi)/screen*switch
			VmatUU(y,x,w,z)=(Vr+im*Vi)/screen*switch
		end do
		close(23)
				
		do i=1,nVDD
			read(24,*) x,y,z,w,Vr,Vi
			VmatDD(x,y,z,w)=(Vr+im*Vi)/screen*switch
			VmatDD(y,x,w,z)=(Vr+im*Vi)/screen*switch
		end do
		close(24)
		
		do i=1,nVUD
			read(25,*) x,y,z,w,Vr,Vi
			VmatUD(x,y,z,w)=(Vr+im*Vi)/screen*switch
		end do
		close(25)
		
	end if
	
end subroutine loadSelfInputs

subroutine whichXYZWsameInd3(xyzw,ind3,exdir,uord)
	! translate xyzw CME indices into 3-component compound index for self energy use
	implicit none
	! INPUT: CME tensor xyzw indices
	integer, intent(in):: xyzw(4)
	! OUTPUTS: 3-index, exchange/direct info, spin of first state (for mixed spin elements)
	integer,intent(out):: ind3(3),exdir,uord

	
	integer x,y,z,w,difs(4),nums(4),i,ile0,exdirs(4),whs(4,3),uords(4)
	integer,allocatable:: inds0(:)
	x=xyzw(1); y=xyzw(2); z=xyzw(3); w=xyzw(4)
	
	! check which indices are the same to eliminate one
	difs=[x-z,x-w,y-z,y-w]
	whs(1,:)=[y,x,w]
	whs(2,:)=[y,x,z]
	whs(3,:)=[x,y,w]
	whs(4,:)=[x,y,z]
	exdirs=[1,0,0,1]
	uords=[0,0,1,1]
	
	nums=[(i,i=1,4)]
	ile0=count(difs.eq.0)
	allocate(inds0(ile0))	
	inds0=pack(nums,difs.eq.0)
	
	ind3=whs(inds0(1),:)
	exdir=exdirs(inds0(1))
	uord=uords(inds0(1))

end subroutine whichXYZWsameInd3



	
end module loadInput
	