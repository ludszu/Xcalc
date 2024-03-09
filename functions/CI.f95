module CI
	use hamiltonian
	use configuration
	
	implicit none
	
	
contains


function makeCIham()result(ham)
	! makes a Hamiltonian matrix
	implicit none
	! OUTPUT: Hamiltonian matrix
	complex*16,allocatable:: ham(:,:)
	
	integer siz,i,j,i2,j2,ileu,iled,ile,ind1,ind2,UD,siz2,ii
	integer co1u(NelU),co2u(NelU),co1d(NelD),co2d(NelD),co1u01(NstVU+NstCU)
	integer co2u01(NstVU+NstCU),co1d01(NstVD+NstCD),co2d01(NstVD+NstCD)
	integer,allocatable:: ehco1u(:,:),ehco2u(:,:),ehco1d(:,:),ehco2d(:,:)
	integer,allocatable:: diffu(:,:),diffd(:,:)
	complex*16 Vdiag,Voff1,Voff2,GSeneV
	complex*16, allocatable:: selfs(:),verts(:)
	real*8 Ediag,Efreez,GSeneSP
	logical UorD,ifOffCI
	
	
	! because we select desired number of excitations,
	! some configurations get rejected, so tha Hamiltonian matrix size is smaller
	siz=noConfU*noConfD
	siz2=count(ifXconfUD)	
	print *, "siz,siz2",siz,siz2
	if (siz.gt.siz2) siz=siz2	
	allocate(ham(siz,siz))
	ham=(0.d0,0.d0)
	
	! if self energies are precalculated, a diagonal vector can be made
	! otherwise - each time we calculate on the fly
	if (ifSelfPreCalc) then
		call makeSelfs(selfs,verts)
	end if
	
	! precalculate GS energy to subtract
	GSeneSP=diagSPCo(GSconfU,GSconfD)
	if (.not.ifLowMem) GSeneV=diagCo(GSconfU,GSconfD)	
	
	print *,"GSeneSP:",GSeneSP
	if (.not.ifLowMem) print *,"GSeneV:",GSeneV
		
	ifOffCI=.false.
	if (ifUseOffDiagVmat.and.ifUseOffDiagSelfLike) ifOffCI=.true.
	
	! loop over all configurations, each spin separately
	print *, "LOOP"
	do i=1,noConfU
		co1u=configsU(i,:)
		co1u01=configs01U(i,:)
		ehco1u=ehconfU(i,:,:)
						
		do i2=1,noConfU
			co2u=configsU(i2,:)
			co2u01=configs01U(i2,:)
			ehco2u=ehconfU(i2,:,:)
			! compare configurations U
			call compareConfigs(co1u01,co2u01,ileu,diffu)
					
			if (i2==1) write(103,*) co1u,-co1d
			
			! if difference between configurations is ok, proceed
			if (ileu.ne.-1) then
				do j=1,noConfD
					co1d=configsD(j,:)
					co1d01=configs01D(j,:)	
					ehco1d=ehconfD(j,:,:)
					ind1=infoConfUD(i,j,1)
					
					! if # of excitations for configuration ij is ok, proceed
					if (ifXconfUD(i,j)) then
					
						do j2=1,noConfD
							co2d=configsD(j2,:)
							co2d01=configs01D(j2,:)	
							ehco2d=ehconfD(j2,:,:)
							ind2=infoConfUD(i2,j2,1)
							! compare configurations D
							call compareConfigs(co1d01,co2d01,iled,diffd)									
							
							! if # of excitations for configuration i2j2 is ok, proceed
							if (ifXconfUD(i2,j2)) then
							
								! if difference between configurations is ok, proceed
								if (iled.ne.-1) then
									ile=ileu+iled
									
									!!! DIAGONAL TERM
									if (ile==0) then
										Ediag=diagSPCo(co1u,co1d) !SP energy
										if (ifSubtractGSene) Ediag=Ediag-GSeneSP ! subtract GS SP energy													
										
										if (ifSelfPreCalc) then
											Vdiag=(0.d0,0.d0)
											if (ifUseSelf) Vdiag=Vdiag+selfs(ind1) ! precalculated self energy
											if (ifUseVertex) Vdiag=Vdiag+verts(ind1) ! precalculated vertex corr.													
										else
											Vdiag=(0.d0,0.d0)
											if (ifUseSelf.and.ifUseVertex) Vdiag=diagCo(co1u,co1d) ! on the fly self energy+ vertex corr.
											if (ifSubtractGSene) Vdiag=Vdiag-GSeneV	! subtract GS Coulomb energy												
										end if
																							
										ham(ind1,ind2)=(Ediag+Vdiag)
										
									!!! OFF DIAGONAL TERM - 1 ELECTRON DIFFERENCE
									else if (ile==1) then
										Voff1=(0.d0,0.d0)
										! spin UP
										if (ileu==1) then
											UorD=.true. 
											if (ifScatterPreCalc) then 
												Voff1=scatterOneOffEH_1X(ehco1u,ehco2u,UorD) ! precalculated offdiagonal scattering
											else
												if (ifOffCI) Voff1=offDiagOneCo(co1u,co1d,diffu,UorD) ! on the fly off diagonal scattering
											end if	

										! spin DN
										else if (iled==1) then
											UorD=.false.	
											if (ifScatterPreCalc) then 
												Voff1=scatterOneOffEH_1X(ehco1d,ehco2d,UorD) ! precalculated offdiagonal scattering
											else
												if (ifOffCI) Voff1=offDiagOneCo(co1u,co1d,diffd,UorD) ! on the fly offdiagonal scattering
											end if				
										else
											Voff1=(0.d0,0.d0)
										end if		
										
										ham(ind1,ind2)=Voff1
										
									!!! OFF DIAGONAL TERM - 2 ELECTRON DIFFERENCE												
									else if (ile==2) then
										! spin U & U
										if (ileu==2) then
											UD=1
											Voff2=offDiagTwoCo(co1u,co1d,diffu,diffd,UD)
										! spin D & D
										else if (iled==2) then
											UD=-1
											Voff2=offDiagTwoCo(co1u,co1d,diffu,diffd,UD)
										! spin 1U & 1D
										else if (ileu==1 .and. iled==1) then
											UD=0
											Voff2=offDiagTwoCo(co1u,co1d,diffu,diffd,UD)
										else
											Voff2=(0.d0,0.d0)
										end if	
										if (ifUseOffDiagVmat) ham(ind1,ind2)=Voff2 ! on the fly 2-electron scattering
										
									else
										ham(ind1,ind2)=(0.d0,0.d0)
									end if
								
								end if
																									
												
							end if
							
							
							deallocate(diffd)
						end do							
								
								
					end if	
								
								
								
				end do
			end if
			
			deallocate(diffu)
		end do
	end do
			
	
	print *, "filled ham"
	write(104,*) real(ham)
	write(104,*) imag(ham)


end function makeCIham 
 
function makeCIhamU()result(ham)
	! makes a Hamiltonian matrix for spin polarised case U
	implicit none
	! OUTPUT: Hamiltonian matrix for spin polarised case U
	complex*16,allocatable:: ham(:,:)
	
	integer siz,i,j,i2,j2,ileu,iled,ile,ind1,ind2,UD
	integer co1u(NelU),co2u(NelU),siz2
	integer co1u01(NstVU+NstCU),co2u01(NstVU+NstCU)
	integer,allocatable:: ehco1u(:,:),ehco2u(:,:)
	integer,allocatable:: diffu(:,:)
	complex*16 Vdiag,Voff1,Voff2,Vfreez,GSeneV,GSeneVCI
	complex*16, allocatable:: selfs(:),verts(:)
	real*8 Ediag,Efreez,GSeneSP
	logical UorD,ifOffCI
	
	! because we select desired number of excitations,
	! some configurations get rejected, so tha Hamiltonian matrix size is smaller
	siz=noConfU
	siz2=count(ifXconfUD)	
	print *, "siz,siz2",siz,siz2!,siz.gt.siz2
	if (siz.gt.siz2) siz=siz2
	allocate(ham(siz,siz))
	ham=(0.d0,0.d0)
	
	! if self energies are precalculated, a diagonal vector can be made
	! otherwise - each time we calculate on the fly
	if (ifSelfPreCalc) then
		call makeSelfsU(selfs,verts)
	end if
	
	! precalculate GS energy to subtract
	GSeneSP=diagSPCo(GSconfU,GSconfU)
	GSeneV=diagCo(GSconfU,GSconfU)
	
	print *,"GSeneSP,GSeneV:",GSeneSP,GSeneV
	
	ifOffCI=.false.
	if (ifUseOffDiagVmat.and.ifUseOffDiagSelfLike) ifOffCI=.true.
	
	! loop over all configurations, each spin separately
	print *, "LOOP"
	do i=1,noConfU
		co1u=configsU(i,:)
		co1u01=configs01U(i,:)
		ehco1u=ehconfU(i,:,:)
		ind1=infoConfUD(i,1,1)
	
		! if # of excitations for configuration i is ok, proceed
		if (ifXconfUD(i,1)) then
		
			do i2=1,noConfU
				co2u=configsU(i2,:)
				co2u01=configs01U(i2,:)
				ehco2u=ehconfU(i2,:,:)
				ind2=infoConfUD(i2,1,1)
				
				! if # of excitations for configuration i2 is ok, proceed
				if (ifXconfUD(i2,1)) then
				
					! compare configurations
					call compareConfigs(co1u01,co2u01,ileu,diffu)
					
					! if difference between configurations is ok, proceed
					if (ileu.ne.-1) then
						
						ile=ileu
						!!! DIAGONAL TERM
						if (ile==0) then
							Ediag=diagSPCo(co1u,co1u) ! SP energy
							if (ifSubtractGSene) Ediag=Ediag-GSeneSP ! subtract GS SP energy
							
							if (ifSelfPreCalc) then
								Vdiag=(0.d0,0.d0)
								if (ifUseSelf) Vdiag=Vdiag+selfs(ind1) ! precalculated self energy
								if (ifUseVertex) Vdiag=Vdiag+verts(ind1) ! precalculated vertex corr.
								
							else
								Vdiag=(0.d0,0.d0)
								if (ifUseSelf.and.ifUseVertex) Vdiag=diagCo(co1u,co1u) ! self energy and vertex corr. on the fly
								if (ifSubtractGSene) Vdiag=Vdiag-GSeneV ! subtract GS Coulomb energy
														
							end if
							
							ham(ind1,ind2)=(Ediag+Vdiag)
														
						!!! OFF DIAGONAL TERM - 1 ELECTRON DIFFERENCE												
						else if (ile==1) then
							Voff1=(0.d0,0.d0)
							if (ileu==1) then
								UorD=.true. 
								if (ifScatterPreCalc) then 
									Voff1=scatterOneOffEH_1X(ehco1u,ehco2u,UorD) ! precalculated offdiagonal scattering
								else
									if (ifOffCI) Voff1=offDiagOneCo(co1u,co1u,diffu,UorD) ! on the fly offdiagonal scattering
								end if	
								
							else
								Voff1=(0.d0,0.d0)
							end if		
							ham(ind1,ind2)=Voff1
														
						!!! OFF DIAGONAL TERM - 2 ELECTRON DIFFERENCE												
						else if (ile==2) then
							if (ileu==2) then
								UD=1
								Voff2=offDiagTwoCo(co1u,co1u,diffu,diffu,UD) ! on the fly 2-electron scattering
							else
								Voff2=(0.d0,0.d0)
							end if	
							if (ifUseOffDiagVmat) ham(ind1,ind2)=Voff2
						else
							ham(ind1,ind2)=(0.d0,0.d0)
						end if
								
								
					end if
					
					deallocate(diffu)
					
				end if
					
					
					
			end do
			
		end if
	end do
	
	write(104,*) real(ham)
	write(104,*) imag(ham)


end function makeCIhamU 
 
function makeCIhamD()result(ham)
	! makes a Hamiltonian matrix for spin polarised case D
	implicit none
	! OUTPUT: Hamiltonian matrix for spin polarised case D
	complex*16,allocatable:: ham(:,:)
	
	integer siz,i,j,i2,j2,ileu,iled,ile,ind1,ind2,UD
	integer co1d(NelD),co2d(NelD),siz2
	integer co1d01(NstVD+NstCD),co2d01(NstVD+NstCD)
	integer,allocatable:: ehco1d(:,:),ehco2d(:,:)
	integer,allocatable:: diffd(:,:)
	complex*16 Vdiag,Voff1,Voff2,GSeneV,GSeneVCI
	real*8 Ediag,GSeneSP
	logical UorD,ifOffCI
	complex*16, allocatable:: selfs(:),verts(:)
	
	! if self energies are precalculated, a diagonal vector can be made
	! otherwise - each time we calculate on the fly
	siz=noConfD
	siz2=count(ifXconfUD)	
	print *, "siz,siz2",siz,siz2!,siz.gt.siz2
	if (siz.gt.siz2) siz=siz2
	allocate(ham(siz,siz))
	ham=(0.d0,0.d0)
	
	! if self energies are precalculated, a diagonal vector can be made
	! otherwise - each time we calculate on the fly
	if (ifSelfPreCalc) then
		call makeSelfsD(selfs,verts)
	end if
	
	! precalculate GS energy to subtract
	GSeneSP=diagSPCo(GSconfD,GSconfD)
	GSeneV=diagCo(GSconfD,GSconfD)
	print *,"GSeneSP,GSeneV:",GSeneSP,GSeneV
	
	ifOffCI=.false.
	if (ifUseOffDiagVmat.and.ifUseOffDiagSelfLike) ifOffCI=.true.
	
	! loop over all configurations, each spin separately
	print *, "LOOP"
	do j=1,noConfD
		co1d=configsD(j,:)
		co1d01=configs01D(j,:)	
		ehco1d=ehconfD(j,:,:)
		ind1=j
		write(103,*) -co1d
			
		! if # of excitations for configuration i is ok, proceed
		if (ifXconfUD(i,1)) then
		
			do j2=1,noConfD
				co2d=configsD(j2,:)
				co2d01=configs01D(j2,:)	
				ehco2d=ehconfD(j2,:,:)
				ind2=j2
				
				! if # of excitations for configuration i is ok, proceed
				if (ifXconfUD(j,1)) then
				
					! compare configurations
					call compareConfigs(co1d01,co2d01,iled,diffd)
										
					
					! if difference between configurations is ok, proceed
					if (iled.ne.-1) then
						ile=iled
						
						!!! DIAGONAL TERM
						if (ile==0) then
							Ediag=diagSPCo(co1d,co1d) ! SP energy
							if (ifSubtractGSene) Ediag=Ediag-GSeneSP ! subtract GS SP energy
							
							if (ifSelfPreCalc) then
								Vdiag=(0.d0,0.d0)
								if (ifUseSelf) Vdiag=Vdiag+selfs(ind1) ! precalculated self energy
								if (ifUseVertex) Vdiag=Vdiag+verts(ind1) ! ! precalculated vertex corr.
								
							else
								Vdiag=(0.d0,0.d0)
								if (ifUseSelf.and.ifUseVertex) Vdiag=diagCo(co1d,co1d) ! self energy and vertex corr. on the fly
								if (ifSubtractGSene) Vdiag=Vdiag-GSeneV ! subtract GS Coulomb energy
							end if
							
							ham(ind1,ind2)=(Ediag+Vdiag)
							
						!!! OFF DIAGONAL TERM - 1 ELECTRON DIFFERENCE												
						else if (ile==1) then
							Voff1=(0.d0,0.d0)
							if (iled==1) then
								UorD=.false.	
								if (ifScatterPreCalc) then 
									Voff1=scatterOneOffEH_1X(ehco1d,ehco2d,UorD) ! precalculated 1-electron scattering
								else
									if (ifOffCI) Voff1=offDiagOneCo(co1d,co1d,diffd,UorD) ! on the fly 1-electron scattering
								end if
								
								
							else
								Voff1=(0.d0,0.d0)
							end if		
							ham(ind1,ind2)=Voff1
						
						!!! OFF DIAGONAL TERM - 2 ELECTRON DIFFERENCE												
						else if (ile==2) then
							if (iled==2) then
								UD=-1
								Voff2=offDiagTwoCo(co1d,co1d,diffd,diffd,UD)! on the fly 2-electron scattering
							else
								Voff2=(0.d0,0.d0)
							end if	
							if (ifUseOffDiagVmat) ham(ind1,ind2)=Voff2
						else
							ham(ind1,ind2)=(0.d0,0.d0)
						end if
					
					end if
					
								
					deallocate(diffd)
					
				end if
					
			end do
			
		end if
			
	end do

	
	write(104,*) real(ham)
	write(104,*) imag(ham)

end function makeCIhamD 
 
subroutine diagCIham(UorD,eigvals,wfs,ifwf)
	! diagonalises Hamiltonian matrix
	implicit none
	! INPUTS: spin U or D, flag to print wavefunctions
	integer,intent(in):: UorD
	logical,intent(in):: ifwf
	! OUTPUTS: eigenvalues, wavefunctions
	real*8,allocatable,intent(out):: eigvals(:)
	complex*16,allocatable,intent(out):: wfs(:,:)
	
	complex*16,allocatable:: ham(:,:)
	
	! LAPACK variables
	CHARACTER*1 :: JOBZ, UPLO
	INTEGER :: N, LDA, LWORK, INFO
	COMPLEX*16,ALLOCATABLE:: A(:,:),WORK(:)
	REAL*8,ALLOCATABLE:: W(:),RWORK(:)
	
	! make Hamiltonian matrix
	if (UorD==-1) then
		ham=makeCIhamD()	
	else if (UorD==1) then
		ham=makeCIhamU()		
	else if (UorD==0) then	
		ham=makeCIham()
	end if
	
	! LAPACK prep
	if (ifwf) then
		JOBZ='V'
	else
		JOBZ='N'
	end if
	UPLO='U'
	N=size(ham,1)
	LDA=N
	LWORK=2*N-1
	ALLOCATE(WORK(LWORK),A(LDA,N),W(N),RWORK(3*N-2))
	
	! output Hamiltonian
	write(11,*) real(ham)
	write(11,*) imag(ham)
	A=ham
	
	deallocate(ham)
	
	! LAPACK DIAGONALISATION
	CALL ZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
		
		
	allocate(eigvals(N))
	eigvals=W
	
	if (ifwf) then
		allocate(wfs(N,N))
		wfs=A	
	end if

end subroutine  
 

end module CI