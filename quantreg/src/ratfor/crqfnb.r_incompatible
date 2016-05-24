# Peng-Huang Censored Quantile Regression
subroutine crqfnb(n,p,a1,c1,n1,X,y,c,B,g,m,r,s,d,u,wn,wp,info)

# Input
#	a1 = p by n1 design matrix (X[c,] transposed)
#	c1 = n1 response vector -y[c] 
#	X  = n by p design matrix
#	y  = n response vector
#	c  = n censoring indicator
#	g  = m grid vector of taus
#	d  = residual n-vector, initialized == 1 
#	s  = initialized cumhaz vector (rep(0,n))
#	u  = upper bound n-vector, initalized == 1
# Workspace
#	r  = rhs p-vector
#	wn = work array n1 by 9
#	wp = work array p by p+3
# Output
#	B  = p by m matrix of estimated coefficients
#	m  = final column dimension of of B
# Note
#	d  is used both for residuals and as a work vector for rqfnb.
#	u  needs to be reinitialized to ones at each iteration
	
integer n,p,n1,m,info,nit(3)
double precision a1(p,n1),c1(n),X(n,p),y(n),c(n),B(p,m),g(m)
double precision wn(n1,9),wp(p,p+3),r(p),s(n),d(n),u(n)
double precision zero,half,one,beta,eps,dH

parameter( zero  =  0.0d0)
parameter( half  =  0.5d0)
parameter( one   =  1.0d0)
parameter( beta  =  0.99995d0)
parameter( eps   =  1.0d-8)

do k = 2,m {
	dH = -log(one - g(k)) + log(one - g(k-1))
	do i = 1,n {
		u(i) = one
		wn(i,1) = half # initialize dual vector
		if(d(i) >= zero) s(i) = s(i) + dH
		d(i) = c(i) - s(i)
		}
	call dgemv('T',n,p,one,X,n,d,1,zero,r,1)
	call rqfnb(n1,p,a1,c1,r,d,u,beta,eps,wn,wp,nit,info)
	if(info != 0) break
	call dcopy(p,wp,1,B(1,k-1),1)
	call dcopy(n,y,1,d,1)
	call dgemv('N',n,p,one,X,n,B(1,k-1),1,one,d,1)
	}
m = k-1
return
end
