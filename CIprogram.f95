program CI_prog
	use hamiltonian
	use loadInput
	use configuration
	use general
	use CI
	
	implicit none
	
	integer i,ile,nVUU,nVUU2,nV,nVUU3
	integer NEtotal,UorD,hamSize,nVDD,nVUD,udborder,j,ix,siz,nVs
	character pathfile*9,foldername*8,swsufix*6,screenType*1,sufix*10
	character*8 fileE,fileV,fileVv,fileV2,fileV3,fileE2,fileE3
	character*8 fileS,fileS2,filPhase,Lfil,Vfil,screenFolder
	real :: start, finish
	real*8 Sz,screen
	complex*16,allocatable:: ham(:,:),wfs(:,:)
	real*8,allocatable:: eigvals(:)
	logical ifwf,preCalc
	
	!!! SETUP	
	
	! work folder and job type
	foldername='wse2_ehX'
	pathfile=foldername//'/'
	preCalc=.false. ! prepare CME input or do CI calculation
	screenType="c" ! c for Coulomb, k for Keldysh
	
	! electron system setup options
	switch=1.d0 ! strength of Coulomb interaction for gradual switching
	screen=1.0 ! dielectric constant
	NEtotal=12 ! total electron number
	Sz=0.d0/2.d0 ! total Sz of electron system	
	nEU=22	! number of SP energy levels spin U
	nED=22	! number of SP energy levels spin D
	Eshift=0.0824 ! shift all SP energies wrt to TB calculation
	
	! electron system setup
	NtotelU=int(NEtotal/2.d0+Sz) ! total spin U electron number 
	NtotelD=int(NEtotal/2.d0-Sz) ! total spin D electron number
	NstVU=NtotelU	! number of filled VB SP states spin U
	NstVD=NtotelD	! number of filled VB SP states spin D
	NstCU=nEU-NstVU	! number of empty CB SP states spin U
	NstCD=nED-NstVD ! number of empty CB SP states spin D
	maxX=1	! max allowed number of electrons excited up, if =NEtotal it's full CI
	minX=1	! min allowed number of electrons excited up
	print *, 'NE,Sz:',NEtotal,Sz,'maxX:',maxX
	print *, 'NelU,NelD:',NtotelU,NtotelD,', NstVU,NstCU:',NstVU,NstCU,', NstVD,NstCD:',NstVD,NstCD
	
	! which Hamiltonian terms to use
	ifUseOffDiagSelfLike=.false. ! offdiagonal terms with summation over filled states
	ifUseOffDiagVmat=.false. ! other offdiagonal terms
	ifUseVertex=.false. ! vertex correction
	ifUseSelf=.false. ! self energy
	ifSelfDirect=.false. ! direct part of self energy, if it's screened
	ifSubtractGSene=.true. ! calculate all wrt to filled GS energy
	
	! has the self energy and scattering terms been precalculated?
	ifSelfPreCalc=.true.
	ifScatterPreCalc=.true.
	
	!!! PREP	
		
	! file names prep
	fileE='Emat.dat'
	fileV='VSOU.dat'
	fileV2='VSOD.dat'
	fileV3='VSUD.dat'
	fileVv='V_SO.dat'
	fileE2='EmaU.dat'
	fileE3='EmaD.dat'
	fileS='selU.dat'
	fileS2='selD.dat'
	Lfil='Lfil.dat'
	Vfil='Vfil.dat'
	
	! output files prep
	write(sufix,'(A2,I2,A3,F3.1)') '_N',NEtotal,'_Sz',abs(Sz)
	write(swsufix,'(A3,F3.1)') '_sw',switch
	
	open(unit=101,file=pathfile//'wfsCI'//sufix//swsufix//'.dat',action='write',status='replace')
	open(unit=102,file=pathfile//'eigCI'//sufix//swsufix//'.dat',action='write',status='replace')
	open(unit=103,file=pathfile//'confs'//sufix//swsufix//'.dat',action='write',status='replace')
	open(unit=1033,file=pathfile//'cosEH'//sufix//swsufix//'.dat',action='write',status='replace')
	open(unit=104,file=pathfile//'hamCI'//sufix//swsufix//'.dat',action='write',status='replace')
	open(unit=11,file=pathfile//'testCI.dat',action='write',status='replace')
	
	! count all time passed
	call cpu_time(start)
    
	! pick screening folder for inputs
	if (screenType=='k') then
		screenfolder='keldysh/'
	else if (screenType=='c') then
		screenfolder='coulomb/'
	end if
	
	if (preCalc) then
	
		!!!changing CME file from CME calculation to CI calculation
		
		udborder=nEU ! when spin U/D SP states are counted together which number ends enumerating U?
		nVs=937024 ! number of all CME to read from input
		! outputs new files to the same folder
		call changeInput(fileVv,fileV3,fileV,fileV2,screenfolder,pathfile,nVs,udborder)
	
	else
		
		!!! BUILDING CONFIGURATIONS
	
		! load angular momentum and valley indices from file
		call loadQNumbers(Lfil,Vfil,pathfile)
		! make and print all possible configurations
		call makeConfigs()
		call printConfigs()
		
		
		!!! LOADING INPUTS
		
		! read number of CME from formatted input files to read
		open(unit=110011,file=pathfile//screenfolder//'/nVs.dat',action='read')
		read(110011,*) nVUU,nVDD,nVUD
		close(110011)
		
		! use formatted CME files to input CME
		nV=nEU !number of SP energy levels used for CME - determines the CME tensor size
		call loadInputs(fileE2,fileE3,fileV,fileV2,fileV3,pathfile,screenFolder,nVUU,nVDD,nVUD,nV,nV,screen)
		print *, "CME load done"
		
		
		! load or precalculate self energy
		if (ifSelfPreCalc) then			
			! opportunity to input converged self energies
		
			nVs=484 ! number of self energy entries to input
			call loadSelfEne1X(fileS,fileS2,pathfile,screenFolder,nVs,nVs,screen)
			print *, "self load done"
		else			
			! these self energies will not be converged
			
			call preCalcSelfs1X(pathfile)
			print *, "self precalc done"
		end if
			
				
		!!! MAIN CI DONE HERE
		
		! only continue to CI once the self energy is precalculated 
		! from step above (nonconverged) or from a separate procedure (possibly converged)
		if (ifSelfPreCalc) then
		
			ifwf=.true. ! if printng wavefunctions
			
			! predetermine Hamiltonian size
			if ((NtotelU>0).and.(NtotelD>0)) then
				UorD=0		
				hamSize=noConfU*noConfD
			else if (NtotelU>0) then
				UorD=1
				hamSize=noConfU
			else if (NtotelD>0) then
				UorD=-1
				hamSize=noConfD
			end if
			print *, 'UorD:',UorD,'ham size:',hamSize
			
			! fill and diagonalise Hamiltonian
			call diagCIham(UorD,eigvals,wfs,ifwf)
			
			! print outputs
			write(102,*) eigvals
			print *, eigvals
			if (ifwf) then				
				write(101,*) real(wfs)
				write(101,*)
				write(101,*) imag(wfs)
			end if
			
			
		end if
		
		
		
		! close open output files
		close(11)
		close(101)
		close(102)
		close(103)
		close(1033)
		close(104)
		
		! finish counting time
		call cpu_time(finish)
		print '("Time = ",f10.3," seconds.")',finish-start
		
	end if
	  
end program
