module hamiltonian
	use configuration
	use loadInput
	implicit none
	
	! GLOBAL variables
	logical ifUseSelf,ifUseVertex,ifUseOffDiagVmat,ifUseOffDiagSelfLike
	logical ifSelfDirect,ifSelfPreCalc,ifScatterPreCalc,ifSubtractGSene
	
	
contains


subroutine preCalcSelfs1X(pathfile)
	! precalculates self energy components of a single excitation (eh pair) 
	! and outputs them to file for CI 
	! it will not be converged
	implicit none
	! INPUT: pathfile to write
	character*9,intent(in)::pathfile
	
	integer nu,nd,i,j,k
	complex*16 Vd,Vx
	
	! run over all filled VB states
	nu=NstVU
	nd=NstVD
	
	print *, "nEU,nED:",nEU,nED
	allocate(SelfEneU(nEU,nEU),SelfEneD(nED,nED))
	SelfEneU=(0.d0,0.d0); SelfEneD=(0.d0,0.d0)
	
	open(unit=333,file=pathfile//'selfPreCalcU.dat',action='write',status='replace')
	open(unit=334,file=pathfile//'selfPreCalcD.dat',action='write',status='replace')
	
	! loop over spin U
	do j=1,nEU
		do k=1,nEU
			do i=1,nu
				Vd=VmatUU(j,i,i,k)
				Vx=VmatUU(j,i,k,i)
				
				if (ifSelfDirect) then
					SelfEneU(j,k)=SelfEneU(j,k)+(Vd-Vx)
				else
					SelfEneU(j,k)=SelfEneU(j,k)+(-Vx)
				end if	
			end do
		
			do i=1,nd
				Vd=VmatUD(j,i,i,k)
				if (ifSelfDirect) then
					SelfEneU(j,k)=SelfEneU(j,k)+(Vd)
				end if
			end do
		
		write(333,*) j,k,real(SelfEneU(j,k)),imag(SelfEneU(j,k))
		
		end do		
	end do
	
	close(333)
	
	! loop over spin D
	do j=1,nED
		do k=1,nED
			do i=1,nd
				Vd=VmatDD(j,i,i,k)
				Vx=VmatDD(j,i,k,i)
				
				if (ifSelfDirect) then
					SelfEneD(j,k)=SelfEneD(j,k)+(Vd-Vx)
				else
					SelfEneD(j,k)=SelfEneD(j,k)+(-Vx)
				end if	
			end do
		
			do i=1,nu
				Vd=VmatUD(i,j,k,i)
				if (ifSelfDirect) then
					SelfEneD(j,k)=SelfEneD(j,k)+(Vd)
				end if
			end do
		
		
		write(334,*) j,k,real(SelfEneD(j,k)),imag(SelfEneD(j,k))
				
		end do		
	end do

	close(334)
	
end subroutine preCalcSelfs1X

subroutine calcSelfOnly1X(ind,pathfile,screenfolder,sufix) ! used externally for self energy convergence
	! precalculates self energy components of a single excitation (eh pair) 
	! and outputs them to file for CI 
	! it may be converged 
	implicit none
	! INPUTS: - folders
	!		  - sufix from main job file
	!		  - ind - starting index of SP states used in CI for which to calculate self energies
	character,intent(in)::pathfile*9,sufix*8,screenfolder*8
	integer,intent(in):: ind
	
	integer nu,nd,i,j,k,n1u,n2u,n1d,n2d,ui,di
	complex*16 Vd,Vx
	character jkstr*10
	
	! run over all filled states
	nu=NstVU
	nd=NstVD	
	print *, "nEU,nED:",nEU,nED
	allocate(SelfEneU(nEU,nEU),SelfEneD(nED,nED))
	SelfEneU=(0.d0,0.d0); SelfEneD=(0.d0,0.d0)
	
	open(unit=333,file=pathfile//screenfolder//'selfCalcU'//sufix//'.dat',action='write',status='replace')
	open(unit=334,file=pathfile//screenfolder//'selfCalcD'//sufix//'.dat',action='write',status='replace')
	
	! print progress: starting index, spin U loop, spin D loop
	print *, "starting ind",ind
	print *, "uloop"
	
	do j=ind,nEU
		do k=ind,nEU
			
			ui=654+j*10+k*100
			write(jkstr,'(A2,I0.3,A2,I0.3)') '_j',j,'_k',k
			open(unit=ui,file=pathfile//screenfolder//'selfUlog'//jkstr//sufix//'.dat',action='write',status='replace')
		
			! sum over spin U VB electrons
			do i=1,nu
				! use 4- or 3-index CME dependent on size
				if (ifLowMem) then
					Vd=VmatUU2(j,i,k,1)
					Vx=VmatUU2(j,i,k,2)
				else
					Vd=VmatUU(j,i,i,k)
					Vx=VmatUU(j,i,k,i)
				end if
				
				! use direct terms in self energy or not?
				if (ifSelfDirect) then
					SelfEneU(j,k)=SelfEneU(j,k)+(Vd-Vx)
				else
					SelfEneU(j,k)=SelfEneU(j,k)+(-Vx)
				end if			
				
				write(ui,*) i,real(Vx),real(SelfEneU(j,k))
				
			end do
		
			! sum over spin D VB electrons
			do i=1,nd
				! use 4- or 3-index CME dependent on size
				if (ifLowMem) then
					Vd=VmatUDU2(j,i,k)
				else
					Vd=VmatUD(j,i,i,k)
				end if
				
				! use direct terms in self energy or not?
				if (ifSelfDirect) then
					SelfEneU(j,k)=SelfEneU(j,k)+(Vd)
				end if
				
			end do
		
		write(333,*) j,k,real(SelfEneU(j,k)),imag(SelfEneU(j,k))
		
		
		close(ui)
		end do		
	end do
	
	close(333)
	
	print *, "dloop"	
	do j=ind,nED
		do k=ind,nED
						
			di=654+j*10+k*100
			write(jkstr,'(A2,I0.3,A2,I0.3)') '_j',j,'_k',k
			open(unit=di,file=pathfile//screenfolder//'selfDlog'//jkstr//sufix//'.dat',action='write',status='replace')
				
			! sum over spin U VB electrons
			do i=1,nd
				if (ifLowMem) then
					Vd=VmatDD2(j,i,k,1)
					Vx=VmatDD2(j,i,k,2)
				else
					Vd=VmatDD(j,i,i,k)
					Vx=VmatDD(j,i,k,i)
				end if
				
				if (ifSelfDirect) then
					SelfEneD(j,k)=SelfEneD(j,k)+(Vd-Vx)
				else
					SelfEneD(j,k)=SelfEneD(j,k)+(-Vx)
				end if			
				
				write(di,*) i,real(Vx),real(SelfEneD(j,k))
				
			end do
		
			! sum over spin D VB electrons
			do i=1,nu
				if (ifLowMem) then
					Vd=VmatDUD2(j,i,k)
				else
					Vd=VmatUD(i,j,k,i)
				end if
				if (ifSelfDirect) then
					SelfEneD(j,k)=SelfEneD(j,k)+(Vd)
				end if
				
			end do
		
		
		write(334,*) j,k,real(SelfEneD(j,k)),imag(SelfEneD(j,k))	
				
		close(di)
				
		end do		
	end do

	close(334)	


end subroutine calcSelfOnly1X

function offDiagTwoCo(coU,coD,diffU,diffD,UorD)result(V)
	! calculates the Hamiltonian matrix element for configurations differing by 2 electron states
	implicit none
	! INPUTS: starting configurations and indices that differ them from final, spin UorD
	integer,intent(in):: diffU(:,:),diffD(:,:),coU(:),coD(:) ! coU,coD given for left configuration
	integer,intent(in):: UorD !-1,0,+1 for D,both,U
	! OUTPUT: Hamiltonian matrix element
	complex*16 V
	
	integer i,j,nu,nd,il,ir,diffl(2),diffr(2),ilel,iler,fac
	complex*16 Vd,Vx
	
	if (UorD.eq.-1) then
		! sorting indices to work out the sign
		diffl=sort(diffD(1,:))
		diffr=sort(diffD(2,:))
		Vd=VmatDD(diffl(1),diffl(2),diffr(2),diffr(1))
		Vx=VmatDD(diffl(1),diffl(2),diffr(1),diffr(2))
		
		! working out the occupied states in the bra and ket
		ilel=count((coD>diffl(1)).and.(coD<diffl(2)).and.(coD.ne.diffl(1)).and.(coD.ne.diffl(2)))
		iler=count((coD>diffr(1)).and.(coD<diffr(2)).and.(coD.ne.diffl(1)).and.(coD.ne.diffl(2)))
		! sign factor based on the odd or even number of operator swaps in the bra and ket 
		if (mod(ilel,2)==mod(iler,2)) then
			fac=1
		else
			fac=-1
		end if
		
		V=fac*(Vd-Vx)
		
	else if (UorD.eq.1) then
		! sorting indices to work out the sign
		diffl=sort(diffU(1,:))
		diffr=sort(diffU(2,:))
		Vd=VmatUU(diffl(1),diffl(2),diffr(2),diffr(1))
		Vx=VmatUU(diffl(1),diffl(2),diffr(1),diffr(2))
		
		! working out the occupied states in the bra and ket
		ilel=count((coU>diffl(1)).and.(coU<diffl(2)).and.(coU.ne.diffl(1)).and.(coU.ne.diffl(2)))
		iler=count((coU>diffr(1)).and.(coU<diffr(2)).and.(coU.ne.diffl(1)).and.(coU.ne.diffl(2)))
		! sign factor based on the odd or even number of operator swaps in the bra and ket 
		if (mod(ilel,2)==mod(iler,2)) then
			fac=1
		else
			fac=-1
		end if
		
		V=fac*(Vd-Vx)
		
	else if (UorD.eq.0) then
		diffl=[diffU(1,1),diffD(1,1)]
		diffr=[diffU(2,1),diffD(2,1)]
		Vd=VmatUD(diffl(1),diffl(2),diffr(2),diffr(1))
		
		! working out the occupied states in the bra and ket
		ilel=count((coU>diffl(1)).and.(coU.ne.diffl(1)))+count((coD<diffl(2)).and.(coD.ne.diffl(2)))
		iler=count((coU>diffr(1)).and.(coU.ne.diffl(1)))+count((coD<diffr(2)).and.(coD.ne.diffl(2)))
		! sign factor based on the odd or even number of operator swaps in the bra and ket 
		if (mod(ilel,2)==mod(iler,2)) then
			fac=1
		else
			fac=-1
		end if		
		
		V=fac*Vd
		
	else
		V=(0.d0,0.d0)
	end if

end function offDiagTwoCo

function offDiagOneCo(coU,coD,diff,UorD)result(V)
	! calculates the Hamiltonian matrix element for configurations differing by 1 electron state
	implicit none
	! INPUTS: starting configurations and indices that differ them from final, spin UorD
	integer,intent(in):: diff(:,:),coU(:),coD(:) ! coU,coD given for left configuration
	logical,intent(in):: UorD
	! OUTPUT: Hamiltonian matrix element
	complex*16 V
	
	integer i,j,nu,nd,il,ir,facu,facd,ilel,iler
	complex*16 Vd,Vx
	
	nu=NelU
	nd=NelD
	il=diff(1,1)
	ir=diff(2,1)
	V=(0.d0,0.d0)
	
	! working out the occupied states in the bra and ket for spin U
	ilel=count(coU>il);	iler=count((coU>ir) .and. (coU.ne.il))
	! sign factor based on the odd or even number of operator swaps in the bra and ket 
		if (mod(ilel,2)==mod(iler,2)) then
		facu=1
	else
		facu=-1
	end if
	
	! working out the occupied states in the bra and ket for spin D
	ilel=count(coD>il);	iler=count((coD>ir) .and. (coD.ne.il))
	! sign factor based on the odd or even number of operator swaps in the bra and ket 
		if (mod(ilel,2)==mod(iler,2)) then
		facd=1
	else
		facd=-1
	end if
	
	! for 1 electron difference we have to sum over all unchanged occupied states
	do i=1,nu		
		if (UorD.and.(coU(i).ne.il)) then	
			Vd=VmatUU(coU(i),il,ir,coU(i))
			Vx=VmatUU(coU(i),il,coU(i),ir)
			V=V+facu*(Vd-Vx)
			
		else if (.not.UorD) then
			Vd=VmatUD(coU(i),il,ir,coU(i))
			V=V+facd*Vd
			
		end if		
	end do
	
	do i=1,nd		
		if ((.not.UorD).and.(coD(i).ne.il)) then			
			Vd=VmatDD(coD(i),il,ir,coD(i))
			Vx=VmatDD(coD(i),il,coD(i),ir)
			V=V+facd*(Vd-Vx)
			
		else if (UorD) then
			Vd=VmatUD(il,coD(i),coD(i),ir)
			V=V+facu*Vd
			
		end if		
	end do	
	
end function OffDiagOneCo

function diagCo(coU,coD)result(V)
	! calculates diagonal Hamiltonian matrix element due to Coulomb interaction 
	! it's self energy and vertex correction, but non converged
	implicit none
	! INPUT: configurations U D
	integer,intent(in):: coU(:),coD(:)
	! OUTPUT: Hamiltonian matrix element
	complex*16 V
	
	integer nu,nd,i,j
	complex*16 Vd,Vx
	
	nu=NelU
	nd=NelD
	V=(0.d0,0.d0)
	
	! sum over all occupied states - all electrons within the configuration
	do i=1,nu
		do j=1,nu
			if (i.ne.j) then
				Vd=VmatUU(coU(i),coU(j),coU(j),coU(i))
				Vx=VmatUU(coU(i),coU(j),coU(i),coU(j))
				
				! if taking direct terms for self energy?
				if (ifSelfDirect) then
					V=V+(Vd-Vx)/2
				else
					V=V+(-Vx)/2 
				end if
				
			end if
		end do
		do j=1,nd
			Vd=VmatUD(coU(i),coD(j),coD(j),coU(i))
			if (ifSelfDirect) then
				V=V+(Vd)
			end if
		end do
	end do
	
	do i=1,nd
		do j=1,nd
			if (i.ne.j) then
				Vd=VmatDD(coD(i),coD(j),coD(j),coD(i))
				Vx=VmatDD(coD(i),coD(j),coD(i),coD(j))
				if (ifSelfDirect) then
					V=V+(Vd-Vx)/2
				else
					V=V+(-Vx)/2
				end if
				
			end if
		end do		
	end do
end function diagCo

function diagSPCo(coU,coD)result(E)
	! calculates diagonal SP Hamiltonian matrix element
	implicit none
	! INPUT: configurations U D
	integer,intent(in):: coU(:),coD(:)
	! OUTPUT: Hamiltonian matrix element
	real*8 E
	
	real*8 Eu,Ed
	
	E=0.d0
	
	Eu=0.d0
	Ed=0.d0
	if (noConfU>0) Eu=sum(EmatU(coU))
	if (noConfD>0) Ed=sum(EmatD(coD))
	E=Eu+Ed
end function diagSPCo

subroutine makeSelfs(selfs,vertcorr)
	! calculate the vector of self energies and vertex correction 
	! to fill the diagonal of the Hamilotnian matrix
	!!! WARNING: only works for 1 excitation currently
	implicit none
	! OUTPUTS: self energies and vertex correction vectors
	complex*16,intent(out),allocatable:: selfs(:),vertcorr(:)
	
	integer siz,i,j,i2,j2,ileu,iled,ile,ind1,UD,siz2,Xile
	integer co1u(NelU),co2u(NelU),co1d(NelD),co2d(NelD),XileU,XileD
	integer co1u01(NstVU+NstCU),co2u01(NstVU+NstCU)
	integer co1d01(NstVD+NstCD),co2d01(NstVD+NstCD)
	logical UorD
	
	! because we select desired number of excitations,
	! some configurations get rejected, so tha Hamiltonian matrix size is smaller
	siz=noConfU*noConfD
	siz2=count(ifXconfUD)	
	print *, "siz,siz2",siz,siz2
	if (siz.gt.siz2) siz=siz2
	
	allocate(selfs(siz),vertcorr(siz))
	selfs=(0.d0,0.d0); vertcorr=(0.d0,0.d0)
	
	! loop over all configurations
	do i=1,noConfU
		co1u=configsU(i,:)
		co1u01=configs01U(i,:)
		XileU=count(ehconfU(i,1,:)>0)			
		
		do j=1,noConfD
			co1d=configsD(j,:)
			co1d01=configs01D(j,:)	
			XileD=count(ehconfD(j,1,:)>0)
			
			ind1=infoConfUD(i,j,1)
			Xile=XileU+XileD
			
			! check if current configuration fulfills our requirements for excitations
			if (ifXconfUD(i,j)) then			
			
				if (Xile==0) then
					selfs(ind1)=(0.d0,0.d0)				
					
				else
					if (XileU==0) then
						UorD=.false.
						selfs(ind1)=selfdiag_1X(ehconfD(j,:,:),UorD)
						vertcorr(ind1)=vertexCorrect_1X(ehconfD(j,:,:),UorD)
					else if (XileD==0) then
						UorD=.true.
						selfs(ind1)=selfdiag_1X(ehconfU(i,:,:),UorD)
						vertcorr(ind1)=vertexCorrect_1X(ehconfU(i,:,:),UorD)
					else
						
						!!fill this out for more than 1 excitations
						
					end if	
				end if
						
			end if		
		end do			
	end do
			
end subroutine makeSelfs

subroutine makeSelfsU(selfs,vertcorr)
	! calculate the vector of self energies and vertex correction 
	! to fill the diagonal of the Hamilotnian matrix, spin polarised system
	!!! WARNING: only works for 1 excitation currently
	implicit none
	! OUTPUTS: self energies and vertex correction vectors
	complex*16,intent(out),allocatable:: selfs(:),vertcorr(:)
	
	integer siz,i,j,i2,j2,ileu,iled,ile,ind1,UD,siz2,ii,Xile
	integer co1u(NelU),co2u(NelU),co1d(NelD),co2d(NelD),XileU,XileD
	integer co1u01(NstVU+NstCU),co2u01(NstVU+NstCU)
	integer co1d01(NstVD+NstCD),co2d01(NstVD+NstCD)
	logical UorD
	
	! because we select desired number of excitations,
	! some configurations get rejected, so tha Hamiltonian matrix size is smaller
	siz=noConfU
	siz2=count(ifXconfUD)	
	print *, "siz,siz2",siz,siz2
	if (siz.gt.siz2) siz=siz2
	
	allocate(selfs(siz),vertcorr(siz))
	selfs=(0.d0,0.d0); vertcorr=(0.d0,0.d0)
	
	! loop over all configurations
	do i=1,noConfU
		co1u=configsU(i,:)
		co1u01=configs01U(i,:)
		XileU=count(ehconfU(i,1,:)>0)	
			
		! check if current configuration fulfills our requirements for excitations
		if (ifXconfUD(i,1)) then			
			ind1=infoConfUD(i,1,1)
		
			if (XileU==0) then
				selfs(ind1)=(0.d0,0.d0)	
			else				
				UorD=.true.
				! currently only works for 1 excitation
				selfs(ind1)=selfdiag_1X(ehconfU(i,:,:),UorD)	
				vertcorr(ind1)=vertexCorrect_1X(ehconfU(j,:,:),UorD)				
			end if		
			
		end if	
	end do
			

end subroutine makeSelfsU

subroutine makeSelfsD(selfs,vertcorr)
	! calculate the vector of self energies and vertex correction 
	! to fill the diagonal of the Hamilotnian matrix, spin polarised system
	!!! WARNING: only works for 1 excitation currently
	implicit none
	! OUTPUTS: self energies and vertex correction vectors
	complex*16,intent(out),allocatable:: selfs(:),vertcorr(:)
	
	integer siz,i,j,i2,j2,ileu,iled,ile,ind1,UD,siz2,ii,Xile
	integer co1u(NelU),co2u(NelU),co1d(NelD),co2d(NelD),XileU,XileD
	integer co1u01(NstVU+NstCU),co2u01(NstVU+NstCU)
	integer co1d01(NstVD+NstCD),co2d01(NstVD+NstCD)
	logical UorD
	
	
	! because we select desired number of excitations,
	! some configurations get rejected, so tha Hamiltonian matrix size is smaller
	siz=noConfD
	siz2=count(ifXconfUD)	
	print *, "siz,siz2",siz,siz2!,siz.gt.siz2
	if (siz.gt.siz2) siz=siz2
	! print *, "siz,siz2",siz,siz2,siz.gt.siz2
	
	allocate(selfs(siz),vertcorr(siz))
	selfs=(0.d0,0.d0); vertcorr=(0.d0,0.d0)
	
	! loop over all configurations
	do i=1,noConfD
		co1d=configsD(i,:)
		co1d01=configs01D(i,:)
		XileD=count(ehconfD(i,1,:)>0)	
			
		! check if current configuration fulfills our requirements for excitations
		if (ifXconfUD(i,1)) then			
			ind1=infoConfUD(i,1,1)
			
			if (XileD==0) then
				selfs(ind1)=(0.d0,0.d0)	
			else				
				UorD=.false.
				! currently only works for 1 excitation
				selfs(ind1)=selfdiag_1X(ehconfD(i,:,:),UorD)	
				vertcorr(ind1)=vertexCorrect_1X(ehconfD(j,:,:),UorD)				
			end if		
			
		end if	
	end do

end subroutine makeSelfsD

function selfdiag_1X(ehco,UorD)result(V)
	! calculates self energy of a single configuration for 1 eh excitation case
	implicit none
	! INPUTS: configuration, spin U or D
	integer, intent(in):: ehco(:,:)
	logical, intent(in):: UorD
	! OUTPUT: self energy
	complex*16 V,Vi
	
	integer i,j,n	
	
	n=count(ehco(1,:)>0) !still, this code is for n=1
	V=(0.d0,0.d0)
	
	if (UorD) then
		do i=1,n
			Vi=SelfEneU(ehco(2,i),ehco(2,i))-SelfEneU(ehco(1,i),ehco(1,i))	
			V=V+Vi
		end do
	else
		do i=1,n
			Vi=SelfEneD(ehco(2,i),ehco(2,i))-SelfEneD(ehco(1,i),ehco(1,i))	
			V=V+Vi
		end do	
	end if


end function selfdiag_1X

function vertexCorrect_1X(ehco,UorD)result(V)
	! calculates vertex correction for a single configuration for 1 eh excitation
	implicit none
	! INPUTS: configuration, spin U or D
	integer, intent(in):: ehco(:,:)
	logical, intent(in):: UorD
	! OUTPUT: vertex correction
	complex*16 V,Vd,Vx
	
	integer i,j,n
	logical ifprint
		
	n=count(ehco(1,:)>0)
	V=(0.d0,0.d0)
	
	if (UorD) then
		do i=1,n
			Vd=VmatUU(ehco(1,i),ehco(2,i),ehco(2,i),ehco(1,i))
			Vx=VmatUU(ehco(1,i),ehco(2,i),ehco(1,i),ehco(2,i))
			V=V-(Vd-Vx)	
		end do
	else
		do i=1,n
			Vd=VmatDD(ehco(1,i),ehco(2,i),ehco(2,i),ehco(1,i))
			Vx=VmatDD(ehco(1,i),ehco(2,i),ehco(1,i),ehco(2,i))
			V=V-(Vd-Vx)
		end do	
	end if

end function vertexCorrect_1X

function scatterOneOffEH_1X(ehco1,ehco2,UorD)result(V)
	! calculates scattering term for a apir of configurations 
	! for 1 eh excitation, includes self-energy-like terms and simple scattering terms
	implicit none
	! INPUTS: two configurations, spin U or D
	integer,intent(in):: ehco1(:,:),ehco2(:,:) ! coU,coD given for left configuration
	logical,intent(in):: UorD
	! OUTPUT: self-energy-like scattering term
	complex*16 V
	
	integer i,j,e1,e2,h1,h2,fac,ilel,iler
	complex*16 Vi,Vj
	logical ifprint,EorH	
	
	V=(0.d0,0.d0)
	
	! configurations can differ only by one particle - E or H
	e1=ehco1(2,1)
	h1=ehco1(1,1)
	e2=ehco2(2,1)
	h2=ehco2(1,1)
	if (e1==e2) then
		EorH=.false.
	else
		EorH=.true.
	end if
	
	! bra and ket contents to figure out the sign of the element
	if (UorD) then
		ilel=NstVU-h1
		iler=NstVU-h2
	else
		ilel=NstVD-h1
		iler=NstVD-h2
	end if
	
	! sign factor based on the odd or even number of operator swaps in the bra and ket 
	fac=1
	if (mod(ilel,2).ne.mod(iler,2)) fac=-1		
	
	! value of the elements
	if (UorD) then	
		Vi=(0.d0,0.d0)
		if (EorH) then
			Vi=SelfEneU(e1,e2)
		else
			Vi=-SelfEneU(h2,h1)
		end if
		
		Vj=VmatUU(e1,h2,h1,e2)-VmatUU(e1,h2,e2,h1)
		if (ifUseOffDiagVmat) V=V-Vj ! simple scattering
		if (ifUseOffDiagSelfLike) V=V+Vi ! self-energy-like scattering
		
	else if (.not.UorD) then
		Vi=(0.d0,0.d0)
		if (EorH) then
			Vi=SelfEneD(e1,e2)
		else
			Vi=-SelfEneD(h2,h1)
		end if
		
		Vj=VmatDD(e1,h2,h1,e2)-VmatDD(e1,h2,e2,h1)
		if (ifUseOffDiagVmat) V=V-Vj ! simple scattering
		if (ifUseOffDiagSelfLike) V=V+Vi ! self-energy-like scattering
		
	end if		
	
	V=V*fac

end function scatterOneOffEH_1X



end module hamiltonian
	