subroutine pack(n,p,match,x,xbar)
integer n,p,match(n)
double precision x(n),xbar(n)
	do i=1,p
		xbar(i)=0e0
	do i=1,n
		xbar(match(i))=xbar(match(i))+x(i)
return
end

subroutine suff(n,p,match,x,y,w,xbar,ybar,wbar,work)
integer n,p,match(n)
double precision x(n),xbar(n),y(n),ybar(n),w(n),wbar(n),work(n)
call pack(n,p,match,w,wbar)
do i=1,n
	xbar(match(i))=x(i)
do i=1,n
	work(i)=y(i)*w(i)
call pack(n,p,match,work,ybar)
	do i=1,p{
		if(wbar(i)>0e0) 
			ybar(i)=ybar(i)/wbar(i) 
		else ybar(i)=0e0

	}
return
end
subroutine unpack(n,p,match,xbar,x)
integer n,p, match(n)
double precision x(n),xbar(p+1)
if(p<n)xbar(p+1)=0e0
do i = 1,n
	x(i)=xbar(match(i))
return
end

