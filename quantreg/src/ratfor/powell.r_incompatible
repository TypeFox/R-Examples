# Revised Version of Fitzenberger's Steepest Descent Algorithm for the Powell Estimator
#
# Roger Koenker:  last revision: 6 February 2008 
#
# Problem:  || min{Ax,c} - b ||_tau  = min!
#
subroutine powell(n,p,p2,A,b,c,x,tau,h,f,U,s,g,d,xh,maxit,nflag)

integer n,p,p2
double precision x(p),A(n,p),b(n),c(n)
double precision f(n),U(p,p),s(n),g(p2),d(p),xh(p)
double precision zero, one, mone, step, tau,pow
integer h(p),hin,hout,k,it,inset,maxit,nflag

PARAMETER(zero = 0.0d0, one = 1.d0, mone = -1.d0)

it = 0

repeat {
it = it + 1
if(it > 1) # Otherwise, assume U is A(h)^{-1}
	call pivot(n,p,h,hin,hout,A,U,d,xh,nflag)
if(nflag > 0)
	{nflag = nflag + 2; return}
do i = 1,p{
	xh(i) = b(h(i))
	}

call dgemv('N',p,p,one,U,p,xh,1,zero,x,1)
call dgemv('N',n,p,one,A,n,x,1,zero,f,1)
do i = 1,n{
	if(inset(p,i,h) > 0 | f(i) > c(i))
		s(i) = zero
	else if(b(i) < f(i))
		s(i) = one - tau
	else
		s(i) = - tau
	}
call dgemv('T',n,p,one,A,n,s,1,zero,xh,1)
call dgemv('T',p,p,one,U,p,xh,1,zero,g,1)
do i = 1,p {
	if(f(h(i)) < c(h(i))) 
		if(b(h(i)) < c(h(i)))
			g(i + p) = - g(i) + one - tau
		else
			g(i + p) = - g(i) - tau
	else
		g(i + p) = - g(i) + tau
	g(i) = g(i) + one - tau
	}
k = idmin(p2,g,1)
if(g(k) >= 0 | it > maxit)
	break
call dscal(p,zero,d,1)
if(k <= p)
	call daxpy(p,one,U(1,k),1,d,1)
else{
	k = k - p
	call daxpy(p,mone,U(1,k),1,d,1)
	}
call dgemv('N',n,p,one,A,n,d,1,zero,s,1)
do i = 1,n {
	call dcopy(p,x,1,xh,1)
	step = (b(i) - f(i))/s(i)
	call daxpy(p,step,d,1,xh,1)
	s(i) = pow(n,p,xh,A,b,c,tau)
	}
hin = idmin(n,s,1)
if(inset(p,hin,h) > 0) 
	{nflag = 2; break}
hout = h(k)
}
if(it > maxit) nflag = 1
return
end

#####################################################
subroutine pivot(n,p,h,hin,hout,A,B,u,v,eflag)

integer n,p,h(p),hin,hout,inset,k,eflag
double precision A(n,p),B(p,p),u(p),v(p)
double precision zero,one

PARAMETER(zero = 0.d0, one =  1.d0) 

eflag = 0
k = inset(p,hout,h)
if(k == 0)
	{eflag = 1; return}
if(inset(p,hin,h) > 0)
	{eflag = 2; return}
if(hin < 1 | hin > n)
	{eflag = 3; return}

call dcopy(p,A(hin,1),n,v,1)
call dgemv('T',p,p,one,B,p,v,1,zero,u,1)
call dcopy(p,B(1,k),1,v,1)
do j = 1,p{
	do i = 1,p{
		if(j == k)
			B(i,j) = B(i,j)/u(k)
		else
			B(i,j) = B(i,j) - (u(j)/u(k)) * v(i)
		}
	}
h(k) = hin
return
end


#####################################################
integer function inset(p,k,h)
integer p,k,h(p)

do inset = 1,p{
	if(h(inset) == k) return
	}
inset = 0
return
end

#####################################################
double precision function pow(n,p,x,A,b,c,tau)
integer n,p
double precision x(p),A(n,p),b(n),c(n)
double precision tau,u,zero,rho,fit

PARAMETER(zero= 0.d0) 

pow = zero
do i = 1,n{
	fit = ddot(p,A(i,1),n,x,1)
	u = b(i)  - min(fit,c(i))
	pow = pow + rho(u, tau)
	}
return
end
#####################################################
double precision function rho(u,tau)
double precision u,tau,one

PARAMETER(one = 1.d0) 

if(u < 0)
	rho = u * (tau - one)
else
	rho = u * tau
return
end
