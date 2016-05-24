subroutine calcvar(nx,n,px,qr,qrank,qpivot,cov,tmpcov,work)

implicit double precision (a-h,o-z)
integer n,px,qrank,qpivot(px)
double precision qr(nx,px),cov(px,px), tmpcov(px,px),work(1)
    
double precision dsum
integer i,j,km
#  compute the unscaled coviance matrix for the linear coefficients
do i=1,qrank{
	do j=1,qrank{
		tmpcov(i,j)=0d0
		cov(i,j)=qr(i,j)
		}
	tmpcov(i,i)=1e0
	}
info=0

#
# notice I have put px in as the row dim of cov
#
call dbksl(cov,px,qrank,tmpcov,px,info)

do i=1,qrank{
	do j=i,qrank{
		dsum=0e0
		km=max(i,j)
		for(k=km;k<=qrank;k=k+1){
			dsum=dsum+tmpcov(i,k)*tmpcov(j,k)
			}
		tmpcov(i,j)=dsum
		tmpcov(j,i)=dsum
		}
	}
# no
do i=1,qrank{
	do j=1,qrank{		
                 cov(i,j)=tmpcov(i,j)
		}
	}
return
end

