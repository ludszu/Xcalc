module configuration
	use general
	implicit none
	
	! GLOBAL variables
	integer:: NtotelU,NtotelD,NstVU,NstVD,NstCD,NstCU,maxX,nEU,nED,minX
	integer NelU,NelD,NfreezU,NfreezD,noConfU,noConfD,sizHam
	integer,allocatable:: vecL(:),vecV(:),configLVU(:,:),configLVD(:,:),ehconfU(:,:,:),ehconfD(:,:,:)
	integer,allocatable:: GSconfU(:),GSconfD(:),GSconf01U(:),GSconf01D(:)
	integer,allocatable:: infoConfUD(:,:,:)
	integer,allocatable:: configsU(:,:),configsD(:,:),configs01U(:,:),configs01D(:,:)
	logical,allocatable:: ifXconfUD(:,:)
	

contains


subroutine makeConfigs()
	implicit none
	integer,allocatable:: confU(:,:),confD(:,:),XsU(:),XsD(:)
	integer,allocatable:: confUtru(:,:),confDtru(:,:),tempU(:,:),tempD(:,:)
	integer,allocatable:: indU(:),indD(:),tempU2(:,:),tempD2(:,:),XileU(:),XileD(:)
	integer,allocatable:: diff(:,:)
	logical,allocatable:: tagU(:),tagD(:)
	integer i,j,NU,ND,sU,sD,ile,xi,xj,xx,ii,iji,Lij,Vij,ileUx,ileDx,NNU,NND
	
	! number of filled states has to be lower or equal to number of total electrons
	if ((NstVU<=NtotelU).and.(NstVD<=NtotelD)) then
		! prep
		NelU=NstVU
		NelD=NstVD
		NU=NstVU+NstCU
		ND=NstVD+NstCD
		NND=NstVD*NstCD+1
		NNU=NstVU*NstCU+1
		
		! create all combinations for spins U and D separately
		if (NelU>0) call combs(NU,NelU,maxX,confU,tagU,XsU)
		if (NelD>0) call combs(ND,NelD,maxX,confD,tagD,XsD)		
		sU=size(confU,1)
		sD=size(confD,1)
		noConfU=0
		noConfD=0
		
		! check which combinations have the required number of excitations
		if (NelU>0) noConfU=count(tagU)
		if (NelD>0) noConfD=count(tagD)
		print *, 'noConfU,noConfD:' ,noConfU,noConfD		
		
		if (NelU>0)  then
			allocate(tempU(sU,(NelU+1)),indU(noConfU))
			tempU(:,1:NelU)=confU
			tempU(:,NelU+1)=XsU
			indU=pack([(i,i=1,sU)],tagU)
			allocate(tempU2(noConfU,NelU+1))
			tempU2=tempU(indU,:)
			deallocate(tempU,indU)		
			allocate(confUtru(noConfU,NelU+1))
			confUtru=sortrows(tempU2,NelU+1)
			deallocate(tempU2)
		end if
		
		if (NelD>0) then
			allocate(tempD(sD,(NelD+1)),indD(noConfD))
			tempD(:,1:NelD)=confD
			tempD(:,NelD+1)=XsD
			indD=pack([(i,i=1,sD)],tagD)
			allocate(tempD2(noConfD,NelD+1))
			tempD2=tempD(indD,:)
			deallocate(tempD,indD)		
			allocate(confDtru(noConfD,NelD+1))
			confDtru=sortrows(tempD2,NelD+1)
			deallocate(tempD2)
		end if
	
		if (NelU>0) allocate(configsU(noConfU,NelU),configLVU(noConfU,2))
		if (NelD>0) allocate(configsD(noConfD,NelD),configLVD(noConfD,2))
		if (NelU>0) configsU=confUtru(:,1:NelU)
		if (NelD>0) configsD=confDtru(:,1:NelD)
		
		if (NelU>0) XileU=confUtru(:,NelU+1)
		if (NelD>0) XileD=confDtru(:,NelD+1)
		
		! GS configuration
		if (NelU>0) GSconfU=configsU(1,1:NelU)
		if (NelD>0) GSconfD=configsD(1,1:NelD)

		! store configurations (1,3,5) and ocupations (1,0,1,0,1)
		if (NelU>0) allocate(configs01U(noConfU,NU))
		if (NelD>0) allocate(configs01D(noConfD,ND))
		
		if (NelU>0) then
			configs01U=0
			do i=1,noConfU
				configs01U(i,confUtru(i,1:NelU))=1		
			end do
		endif
		if (NelD>0) then
			configs01D=0
			do i=1,noConfD
				configs01D(i,confDtru(i,1:NelD))=1		
			end do
		end if
		
		if (NelD>0) deallocate(confDtru)
		if (NelU>0) deallocate(confUtru)
		
		! GS occupation
		if (NelU>0) GSconf01U=configs01U(1,:)
		if (NelD>0) GSconf01D=configs01D(1,:)

		! store also indices of excitations: init->fin, H->E
		! (1,3,5)=(1,0,1,0,1)=(1,1,1,0,0)->(1,0,1,0,1)=H:2->E:5
		if (NelU>0) then
			allocate(ehconfU(noConfU,2,maxX))
			ehconfU=0.d0
			do i=2,noConfU
				call giveEHconfInd(GSconf01U,configs01U(i,:),ile,diff)
				ehconfU(i,1,:)=diff(1,:)
				ehconfU(i,2,:)=diff(2,:)				
			end do
		endif

		if (NelD>0) then
			allocate(ehconfD(noConfD,2,maxX))
			ehconfD=0.d0
			do i=2,noConfD
				call giveEHconfInd(GSconf01D,configs01D(i,:),ile,diff)
				ehconfD(i,1,:)=diff(1,:)
				ehconfD(i,2,:)=diff(2,:)				
			end do
		endif

		! store quantum numbers: angular momentum and valley index
		if (NelU>0) configLVU=0 
		if (NelD>0) configLVD=0
		do i=1,noConfU
			do j=1,NelU
				configLVU(i,1)=configLVU(i,1)+vecL(configsU(i,j))
				configLVU(i,2)=configLVU(i,2)+vecV(configsU(i,j))
			end do
		end do
		do i=1,noConfD
			do j=1,NelD
				configLVD(i,1)=configLVD(i,1)-vecL(configsD(i,j))
				configLVD(i,2)=configLVD(i,2)-vecV(configsD(i,j))
			end do
		end do

		! make an info viarbale with all details on configurations
		if ((NelD>0).and.(NelU>0)) then
		
			allocate(ifXconfUD(noConfU,noConfD),infoConfUD(noConfU,noConfD,4))
			ifXconfUD=.false.
			
			ii=1
			do i=1,noConfU
				xi=XileU(i)
				do j=1,noConfD
					xj=XileD(j)
					
					xx=xi+xj
					infoConfUD(i,j,2)=xx
					Lij=configLVU(i,1)+configLVD(j,1)
					Vij=configLVU(i,2)+configLVD(j,2)
					infoConfUD(i,j,3:4)=[Lij,Vij]
					if ((xx.le.maxX).and.(xx.ge.minX)) then
						ifXconfUD(i,j)=.true.
						infoConfUD(i,j,1)=ii
						ii=ii+1
					else
						infoConfUD(i,j,1)=0					
					end if
					
				end do				
			end do
			
				
		else if (NelU>0) then
			allocate(ifXconfUD(noConfU,1),infoConfUD(noConfU,1,4))
			infoConfUD=0
			ifXconfUD=.true.
			where ((XileU.gt.maxX).or.(XileU.lt.minX))	ifXconfUD(:,1)=.false.
			ileUx=count(ifXconfUD)
			
			ii=1
			do i=1,noConfU
				if ((XileU(i).le.maxX).and.(XileU(i).ge.minX)) then
					infoConfUD(i,1,1)=ii
					ii=ii+1
				end if
			end do
			
			infoConfUD(:,1,2)=XileU
			infoConfUD(:,1,3)=configLVU(:,1)
			infoConfUD(:,1,4)=configLVU(i,2)
		else if (NelD>0) then
			allocate(ifXconfUD(noConfD,1),infoConfUD(noConfD,1,4))
			infoConfUD=0
			ifXconfUD=.true.
			where ((XileD.gt.maxX).or.(XileD.lt.minX))	ifXconfUD(:,1)=.false.
			ileDx=count(ifXconfUD)
			
			ii=1
			do i=1,noConfD
				if ((XileD(i).le.maxX).and.(XileD(i).ge.minX)) then
					infoConfUD(i,1,1)=ii
					ii=ii+1
				end if
			end do
			
			infoConfUD(:,1,2)=XileD
			infoConfUD(:,1,3)=configLVD(:,1)
			infoConfUD(:,1,4)=configLVD(i,2)
		end if
		
		sizHam=count(ifXconfUD)		
		print *, "sizHam",sizHam
		
	else
		print *, 'NstVU or NstVD too big'
	end if

end subroutine makeConfigs

subroutine giveEHconfInd(co101,co201,ile,diff)
	! finds number and indices of excitations for a given occupation co201 wrt co101
	implicit none
	! INPUT: two occupations, usually one is the GS
	integer,intent(in):: co101(:),co201(:)
	! OUTPUT: number of excitations ile, indices of excitations diff
	integer,intent(out):: ile
	integer,allocatable,intent(out)::diff(:,:)
	
	integer,allocatable:: codiff(:)
	integer n1,n2,i
		
	n1=size(co101)
	n2=size(co201)
	
	if (n1==n2) then
		allocate(codiff(n1))
		codiff=co101-co201
		ile=count(codiff.ne.0)/2
		
		allocate(diff(2,ile))
		diff(1,:)=pack([(i,i=1,n1)],codiff>0)
		diff(2,:)=pack([(i,i=1,n1)],codiff<0)
		
	else
		ile=-1
	end if

end subroutine giveEHconfInd

subroutine compareConfigs(co101,co201,ile,diff)
	! conpares two electron configurations: 
	! finds number of electrons on different SP states and indices of these states
	implicit none
	! INPUT: two occupations
	integer,intent(in):: co101(:),co201(:)
	! OUTPUT: number of different electron states ile and indices diff
	integer,intent(out):: ile
	integer,allocatable,intent(out)::diff(:,:)
	
	integer,allocatable:: codiff(:)
	integer n1,n2,i
	
	n1=size(co101)
	n2=size(co201)
	
	if (n1==n2) then
		allocate(codiff(n1))
		codiff=co101-co201
		ile=count(codiff.ne.0)/2
		if (ile==0) then
			allocate(diff(1,1))
			diff=0
		else if (ile==1) then
			allocate(diff(2,1))
			diff(1,:)=pack([(i,i=1,n1)],codiff>0)
			diff(2,:)=pack([(i,i=1,n1)],codiff<0)
		else if (ile==2) then
			allocate(diff(2,2))
			diff(1,:)=pack([(i,i=1,n1)],codiff>0)
			diff(2,:)=pack([(i,i=1,n1)],codiff<0)
		else
			allocate(diff(1,1))
			diff=0
		end if
	else
		ile=-1
	end if

end subroutine  compareConfigs

subroutine printConfigs()
	implicit none
	character lstr*5,vstr*5,xstr*5,ehstr1*13,formt*29,ixstr*1
	character confstr*31,ehstr2*13,ehstr*28
	integer ix,i,j
	
	
	print *, " "
		print *, 'configUD:'
		do i=1,noConfU
			do j=1,noConfD
				if (ifXconfUD(i,j)) then
					ix=maxX
					write(ixstr,'(I1)') maxX
					write(confstr,'(A2,6I2,A5,6I2)') "U:",configsU(i,:),"   D:",configsD(j,:)
					formt="(A4,"//ixstr//"I2,A5,"//ixstr//"I2)"
					write(ehstr1,formt) "h U:",ehconfU(i,1,1:ix)," e U:",ehconfU(i,2,1:ix)
					write(ehstr2,formt)"h D:",ehconfD(j,1,1:ix)," e D:",ehconfD(j,2,1:ix)
					ehstr=ehstr1//"  "//ehstr2 
	
					write(1033,*) ehconfD(j,1,1:ix),ehconfU(i,1,1:ix),ehconfD(j,2,1:ix),ehconfU(i,2,1:ix)
	
					write(xstr,'(A2,I3)') "X=",infoConfUD(i,j,2)
					write(lstr,'(A2,I3)') "L=",infoConfUD(i,j,3)
					write(vstr,'(A2,I3)') "V=",infoConfUD(i,j,4)
					print *, infoConfUD(i,j,1),":   ",confstr,"   ",xstr,"   ",ehstr
				end if
			end do
		end do
	

end subroutine printConfigs





end module configuration