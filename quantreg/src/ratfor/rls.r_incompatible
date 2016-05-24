#This is a simple recursive least squares routine Reference:  Harvey TSM p. 100
#
subroutine rls(n,p,x,y,b,A,Ax)
integer n,p
double precision x(p,n),y(n),b(p,n),A(p,p),Ax(p)
double precision zero,one,mone,f,r,ddot
data one/1.d0/
data mone/-1.d0/
data zero/0.d0/
#
#On input:
#
#	A = crossprod(x[1:p,1:p))^{-1} 
#	b(,p) = A crossprod(x[1:p,1:p],y[1:p])  
#
do i = (p+1),n {
	call dgemv('N',p,p,one,A,p,x(1,i),1,zero,Ax,1)
	f = one + ddot(p,x(1,i),1,Ax,1)
	r = (y(i)-ddot(p,x(1,i),1,b(1,i-1),1))/f
	call daxpy(p,one,b(1,i-1),1,b(1,i),1)
	call daxpy(p,r,Ax,1,b(1,i),1)
	call dger(p,p,mone/f,Ax,1,Ax,1,A,p)
	}
return
end
		

