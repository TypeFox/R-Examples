#parzen, wei and ying's bootstrap
subroutine pwy(m,n,k,m5,n2,a,c,b,t,toler,ift,x,e,s, wa,wb)
double precision b(m),a(k,n),x(n,k)
double precision wa(m5,n2),wb(m),e(m),c(m,n)
double precision t,toler
integer m,n,k,m5,n2,ift
integer s(m)
do i=1,k{
	call dcopy(n,a(i,1),k,c(m,1),m)
	call rq0(m,n,m5,n2,c,b,t,toler,ift,x(1,i),e,s,wa,wb)
	}
return
end
#ratfor outer loop for xy-pairs rq bootstrap
#notation is horrendous 
#   ratfor   R-function
#______________________
#	m -> n  number of original obs
#	n -> p  number of parameters
#	k -> R  number of BS replications
#	mofn -> m  number of BS observations
#
subroutine xys(mofn,m,n,k,mofn5,n2,a,b,tau,toler,ift,x,e,s, wa,wb,aa,bb,ss)
double precision b(m),a(m,n),x(n,k)
double precision wa(mofn5,n2),wb(mofn)
double precision aa(mofn,n),bb(mofn),e(mofn)
double precision tau,toler
integer ss(mofn,k),s(mofn),mofn,m,n,k,mofn5,n2,ift(k)
do i=1,k {
	do ii=1,mofn{
		bb(ii)=b(ss(ii,i))
		do jj=1,n{
			aa(ii,jj)=a(ss(ii,i),jj)
			}
		}
	call rq0(mofn,n,mofn5,n2,aa,bb,tau,toler,ift(i),x(1,i),e,s,wa,wb)
	}
return
end
# Weighted (Bose) Bootstrap version
subroutine wxy(m,n,k,m5,n2,a,b,tau,toler,ift,x,e,s,wa,wb,aa,bb,w)
double precision b(m),a(m,n),x(n,k)
double precision w(m,k),wa(m5,n2),wb(m)
double precision aa(m,n),bb(m),e(m)
double precision tau,toler
integer s(m),m,n,k,m5,n2,ift(k)
do i=1,k {
	do ii=1,m{
		bb(ii)=b(ii)*w(ii,i)
		do jj=1,n{
			aa(ii,jj)=a(ii,jj)*w(ii,i)
			}
		}
	call rq0(m,n,m5,n2,aa,bb,tau,toler,ift(i),x(1,i),e,s,wa,wb)
	}
return
end

#does a matrix multiply to make Y matrix for heqf bootstrap
subroutine heqfy(n,p,r,x,b,y)
integer n,p,r
double precision x(n,p),b(p,n,r),y(n,r)
do i=1,r{
	do j=1,n{
		y(j,i)=ddot(p,x(j,1),n,b(1,j,i),1)
		}
	}
return
end
