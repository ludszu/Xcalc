module general
	implicit none
	! parameters
	complex*16,parameter:: im=(0.0,1.0)
	real*8,parameter:: pi=3.1415926
	
contains
	
function choose(n, k)
    implicit none
	integer :: choose
    integer, intent(in) :: n, k
    
 
    integer :: imax, i, imin
 
    if ( (n < 0 ) .or. (k < 0 ) ) then
       print *, "negative in choose"
       choose = 0    
    else
       if ( n < k ) then
          choose = 0
       else if ( n == k ) then
          choose = 1
       else
          imax = max(k, n-k)
          imin = min(k, n-k)
          choose = 1
          do i = imax+1, n
             choose = choose * i
          end do
          do i = 2, imin
             choose = choose / i
          end do
       end if
    end if
end function choose
      
subroutine combs(n, k, nX, co,tags,Xs)
    implicit none
	integer, intent(in) :: n, k,nX
    integer,allocatable, intent(out) :: co(:,:),Xs(:)
	logical,allocatable,intent(out):: tags(:)
	
    integer :: i, j, s, ix, kx, hm, t
    
	
    hm = choose(n, k)
    
 
    allocate(co(hm,k),tags(hm),Xs(hm))
    
	tags=.true.
	Xs=0
    do i = 1, hm
	   ix = i-1; kx = k
	   do s = 1, n
		  if ( kx == 0 ) then
			  exit
		  else
			  t = choose(n-s, kx-1)
			  if ( ix < t ) then
				 co(i,k-kx+1) = s
				 kx = kx - 1
			  else
				 ix = ix - t			 
			  end if
		  end if
	   end do
	   if (k>nX) then
			if (co(i,k-nX)>k) tags(i)=.false.
	   end if
	   Xs(i)=count(co(i,:)>k)
    end do
 
end subroutine combs
 
function unique(vec)result(uniq)
    implicit none
    integer,intent(in):: vec(:)
    integer, dimension(:),allocatable:: uniq,uni
	integer :: i, min_val, max_val
    
	
	i=0	
	allocate(uni(size(vec)))
    min_val = minval(vec)-1
    max_val = maxval(vec)
	uni=0
    do while (min_val<max_val)
        i = i+1
        min_val = minval(vec, mask=vec>min_val)
        uni(i) = min_val
	enddo
	allocate(uniq(i))
	uniq=uni(1:i)
	
end function unique
	
subroutine uniqueInd(vec,uniq,indA,indB)
    implicit none
    real*8,intent(in):: vec(:)
    integer, dimension(:),allocatable,intent(out):: indA,indB
	real*8,allocatable:: uniq(:),uni(:)
	integer, dimension(:),allocatable::ind1
	integer :: i, min_ind,m
	real*8 min_val, max_val
	logical,allocatable:: temp(:)
    
	
	i=0	
	m=size(vec)
	allocate(uni(m),ind1(m),indB(m),temp(m))
    min_val = minval(vec)-1
    max_val = maxval(vec)
	uni=0
	ind1=0
	indB=0
    do while (min_val<max_val)
        i = i+1
        min_ind = minloc(vec,dim=1, mask=vec>min_val)
		min_val = minval(vec, mask=vec>min_val)
		temp=vec==min_val
		where (temp) indB=i
        uni(i) = min_val
		ind1(i) = min_ind
		 
    enddo
	allocate(uniq(i),indA(i))
	uniq=uni(1:i)
	indA=ind1(1:i)
	
end subroutine uniqueInd	
	
function binomial(n,k)result(binom)
	implicit none
	integer, intent(in):: n,k
	integer binom
	
	integer i,prod
	
	if ((k.ne.0).and.(k.lt.n)) then
		prod=1
		do i=1,(n-k)
			prod=prod*(k+i)/i	
		end do
		binom=prod
	else if ((k.eq.0).or.(k.eq.n)) then 
		binom=1
	else
		binom=0
	end if
	
end function binomial
	
function sub2ind(siz,row,col)result(ind)

	implicit none
	integer,intent(in)::siz(2)
	integer,intent(in),dimension(:):: row,col
	
	integer,allocatable,dimension(:):: ind
	
	allocate(ind(size(row)))
	ind=(col-1)*siz(1)+row

end function sub2ind

function sub2ind_mD(siz,rowcols)result(ind)

	implicit none
	integer,intent(in)::siz(:)
	integer,intent(in),dimension(:,:):: rowcols
	
	integer(kind=8),allocatable,dimension(:):: ind
	integer mdim,ile,i,x,y,z,w
	integer(kind=8) n1,n2,n3,indi
	
	mdim=size(siz,1)
	ile=size(rowcols,1)
	
	n1=siz(1)
	n2=siz(1)*siz(2)
	n3=siz(1)*siz(2)*siz(3)
	
	allocate(ind(ile))
	
	do i=1,ile
		x=rowcols(i,1)
		y=rowcols(i,2)
		z=rowcols(i,3)
		w=rowcols(i,4)
		
		indi=(x-1)*n3+(y-1)*n2+(z-1)*n1+w
		ind(i)=indi
	
	end do
	
end function sub2ind_mD

subroutine ind2sub(siz,inde,row,col)

	implicit none
	integer,intent(in)::siz(2)
	integer,intent(in),dimension(:):: inde
	
	integer,allocatable,dimension(:),intent(out):: row,col
	
	allocate(row(size(inde)),col(size(inde)))
	col=(inde-1)/siz(1)+1
	row=inde-(col-1)*siz(1)
	
end subroutine ind2sub	

subroutine ind2sub_mD(siz,inde,subs)
	implicit none
	integer,intent(in)::siz(:)
	integer(kind=8),intent(in),dimension(:):: inde
	
	integer,allocatable,dimension(:,:),intent(out):: subs
	
	integer mdim,ile,i,x,y,z,w
	integer(kind=8) n1,n2,n3,ii
	
	mdim=size(siz,1)
	ile=size(inde,1)
	
	allocate(subs(ile,mdim))
	subs=0
		
	do i=1,ile
	
		ii=inde(i)
		
		n1=siz(1)
		n2=siz(1)*siz(2)
		n3=siz(1)*siz(2)*siz(3)
		
		
		x=(ii-1)/n3+1
		y=(ii-(x-1)*n3-1)/n2+1
		z=(ii-(x-1)*n3-(y-1)*n2-1)/n1+1
		w=ii-(x-1)*n3-(y-1)*n2-(z-1)*n1
		
		subs(i,:)=[x,y,z,w]
				
	end do

end subroutine ind2sub_mD

function sort(A) result(B)
	implicit none
    integer,intent(in):: A(:)
	integer,allocatable:: B(:)
	
    integer:: nsize,irow,krow,siz,buf

    siz=size(A)
	
	allocate(B(siz))
	B=A	
		
    do irow = 1, siz
		krow = minloc(B(irow:siz),dim=1) + irow - 1
		buf = B(irow)
		B(irow) = B(krow)
        B(krow) = buf
    enddo
		
end function sort

function sortrows(A,Ncol) result(B)
	implicit none
    integer,intent(in):: A(:,:),Ncol
	integer,allocatable:: B(:,:)
	
    integer:: nsize,irow,krow,shap(2)
	integer,allocatable:: buf(:)

    shap=shape(A)
	
	allocate(buf(shap(2)))
	allocate(B(shap(1),shap(2)))
	B=A
	
		
    do irow = 1, shap(1)
		krow = minloc( B( irow:shap(1), Ncol ), dim=1 ) + irow - 1
		buf( : )     = B( irow, : )
		B( irow, : ) = B( krow, : )
        B( krow, : ) = buf( : )
    enddo
		
end

subroutine cart2polar(x,y,r,fi)
implicit none
	real*8,intent(in):: x,y	
	real*8,intent(out):: r,fi
	
	
	if (x<0) then
		fi=atan(y/x)+pi
	else if (x>0) then
		fi=atan(y/x)
	else if (x.eq.0) then
		if (y>0) then
			fi=pi/2
		else if (y<0) then
			fi=-pi/2
		else if (y.eq.0) then
			fi=0.0
		end if
	end if
	
	r=sqrt(x**2+y**2)
	fi=fi*180.0/pi

end subroutine cart2polar

function kron(a,b)	result(ab)
	implicit none
	real*8,intent(in):: a(:,:)
	complex*16,intent(in):: b(:,:)		
	complex*16, allocatable:: ab(:,:)	
	integer r,ra,rb,c,ca,cb,i,j
	
	ra = ubound(a,dim = 1)	
	ca = ubound(a,dim = 2)	
	rb = ubound(b,dim = 1)	
	cb = ubound(b,dim = 2)	
	 
	allocate (ab(ra*rb,ca*cb))	
	ab=(0.d0,0.d0)
	r = 0		
	do i = 1,ra	
		c = 0		
		do j = 1,ca		
			ab(r + 1:r + rb,c + 1:c + cb) = a(i,j)*b
			c = c + cb		
		end do		
		r = r + rb		
	end do	
		  
end function kron	

function diagMat(vec,d) result(mat)
	implicit none
	real*8,intent(in):: vec(:)
	integer,intent(in):: d
	real*8,allocatable:: mat(:,:)
	
	integer m,n,i,ad
	
	
	m=size(vec)
	ad=abs(d)
	n=ad+m
	allocate(mat(n,n))
	mat=0.0
	
	if (d.ge.0) then
		do i=1,m
			mat(i,ad+i)=vec(i)		
		end do
	else
		do i=1,m
			mat(ad+i,i)=vec(i)		
		end do	
	end if
	
end function diagMat

	
 
 
end module general
