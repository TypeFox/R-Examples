# This is a primal-dual log barrier form of the
# interior point LP solver of Lustig, Marsten and Shanno ORSA J Opt 1992.
# It is a projected Newton primal-dual logarithmic barrier method which uses
# the predictor-corrector approach of Mehrotra for the mu steps.
# For the sake of brevity we will call it a Frisch-Newton algorithm.
# The primary difference between this code and the previous version fna.r is
# that we assume a feasible starting point so p,d,b feasibility gaps = 0.
# Problem:
# 	min c'x s.t. Ax=b, 0<=x<=u
# 
# The linear system we are trying to solve at each interation is:
# 		 Adx = 0
# 	     dx + ds = 0
# 	A'dy -dw +dz = 0
# 	   Xdz + Zdx = me - XZe - DXDZe
# 	   Sdw + Wds = me - SWe - DSDWe 
# But the algorithm proceeds in two steps the first of which is to solve:
# 		 Adx = 0
# 	     dx + ds = 0
# 	A'dy -dw +dz = 0
# 	   Xdz + Zdx =  - XZe 
# 	   Sdw + Wds =  - SWe
# and then to make some refinement of mu and modify the implied Newton step.
# Denote dx,dy,dw,ds,dz as the steps for the respective variables x,y,w,s,z, and
# the corresponding upper case letters are the diagonal matrices respectively.
# 
# To illustrate the use of the function we include a calling routine to
# compute solutions to the linear quantile regression estimation problem.
# See the associated S function rqfn for details on the calling sequence.
# On input:
# 	a is the p by n matrix X'
# 	y is the n-vector of responses
# 	u is the n-vector of upper bounds 
#       d is an n-vector of ones
#       wn is an n-vector of ones in the first n elements
# 	beta is a scaling constant, conventionally .99995
# 	eps is a convergence tolerance, conventionally 1d-07
# On output:
# 	a,y are unaltered
# 	wp contains the solution  coefficient vector in the first p elements
# 	wn contains the residual vector in the first n elements
# 
# 
subroutine rqfn(n,p,a,y,rhs,d,u,beta,eps,wn,wp,aa,nit,info)
integer n,p,info,nit(3)
double precision a(p,n),y(n),rhs(p),d(n),u(n),wn(n,10),wp(p,p+3),aa(p,p)
double precision one,beta,eps
parameter( one = 1.0d0)
call fna(n,p,a,y,rhs,d,u,beta,eps,wn(1,1),wn(1,2),
	wp(1,1),wn(1,3),wn(1,4),wn(1,5), wn(1,6),
	wp(1,2),wn(1,7),wn(1,8),wn(1,9),wn(1,10),wp(1,3), wp(1,4),aa,nit,info)
return
end
subroutine fna(n,p,a,c,b,d,u,beta,eps,x,s,y,z,w,
		dx,ds,dy,dz,dw,dsdw,dxdz,rhs,ada,aa,nit,info)

integer n,p,pp,i,info,nit(3)
double precision a(p,n),c(n),b(p)
double precision zero,one,mone,big,ddot,dmax1,dmin1,dasum
double precision deltap,deltad,beta,eps,cx,by,uw,uz,mu,mua,acomp,rdg,g
double precision x(n),u(n),s(n),y(p),z(n),w(n),d(n),rhs(p),ada(p,p)
double precision aa(p,p),dx(n),ds(n),dy(p),dz(n),dw(n),dxdz(n),dsdw(n)

parameter( zero  = 0.0d0)
parameter( half  = 0.5d0)
parameter( one   = 1.0d0)
parameter( mone   = -1.0d0)
parameter( big  = 1.0d+20)

# Initialization:  We try to follow the notation of LMS
# On input we require:
# 
# 	c = n-vector of marginal costs (-y in the rq problem)
# 	a = p by n matrix of linear constraints (x' in rq)
# 	b = p-vector of rhs ((1-tau)x'e in rq)
#       u = upper bound vector ( e in rq)
# 	beta = barrier parameter, LMS recommend .99995
# 	eps = convergence tolerance, LMS recommend 10d-8
# 	
# 	the integer vector nit returns iteration counts
# 	the integer info contains an error code from the Cholesky in stepy
# 		info = 0 is fine
# 		info < 0 invalid argument to dposv 
# 		info > 0 singular matrix
nit(1)=0
nit(2)=0
nit(3)=n
pp=p*p
# Start at the OLS estimate for the parameters
call dgemv('N',p,n,one,a,p,c,1,zero,y,1)
call stepy(n,p,a,d,y,aa,info)
if(info != 0) return
# Save sqrt of aa' for future use for confidence band
do i=1,p{
	do j=1,p
		ada(i,j)=zero
        ada(i,i)=one
	}
call dtrtrs('U','T','N',p,p,aa,p,ada,p,info)
call dcopy(pp,ada,1,aa,1)
# put current residual vector in s (temporarily)
call dcopy(n,c,1,s,1)
call dgemv('T',p,n,mone,a,p,y,1,one,s,1)
# Initialize remaining variables
# N.B. x must be initialized on input: for rq as (one-tau) in call coordinates
do i=1,n{
	d(i)=one
	if(dabs(s(i)) < eps){
		z(i) = dmax1( s(i),zero) + eps
		w(i) = dmax1(-s(i),zero) + eps
		}
	else {
		z(i) = dmax1( s(i),zero)
		w(i) = dmax1(-s(i),zero) 
		}
	s(i)=u(i)-x(i)
	}
cx = ddot(n,c,1,x,1)
by = ddot(p,b,1,y,1)
uw = dasum(n,w,1)
uz = dasum(n,z,1)
# rdg =  (cx - by + uw)/(one + dabs( by - uw))
# rdg =  (cx - by + uw)/(one + uz + uw)
rdg =  (cx - by + uw)
while(rdg > eps) {
	nit(1)=nit(1)+1
	do i =1,n{
		d(i) = one/(z(i)/x(i) + w(i)/s(i))
		ds(i)=z(i)-w(i)
		dx(i)=d(i)*ds(i)
		}
	call dgemv('N',p,n,one,a,p,dx,1,zero,dy,1)#rhs
	call dcopy(p,dy,1,rhs,1)#save rhs
	call stepy(n,p,a,d,dy,ada,info)
	if(info != 0) return
	call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1)
	deltap=big
	deltad=big
	do i=1,n{
		dx(i)=d(i)*ds(i)
		ds(i)=-dx(i)
		dz(i)=-z(i)*(dx(i)/x(i) + one)
		dw(i)=w(i)*(dx(i)/s(i) - one)
		dxdz(i)=dx(i)*dz(i)
		dsdw(i)=ds(i)*dw(i)
		if(dx(i)<0)deltap=dmin1(deltap,-x(i)/dx(i))
		if(ds(i)<0)deltap=dmin1(deltap,-s(i)/ds(i))
		if(dz(i)<0)deltad=dmin1(deltad,-z(i)/dz(i))
		if(dw(i)<0)deltad=dmin1(deltad,-w(i)/dw(i))
		}
	deltap=dmin1(beta*deltap,one)
	deltad=dmin1(beta*deltad,one)
	if(deltap*deltad<one){
		nit(2)=nit(2)+1
		acomp=ddot(n,x,1,z,1)+ddot(n,s,1,w,1)
		g=acomp+deltap*ddot(n,dx,1,z,1)+
			deltad*ddot(n,dz,1,x,1)+ 
			deltap*deltad*ddot(n,dz,1,dx,1)+
			deltap*ddot(n,ds,1,w,1)+
			deltad*ddot(n,dw,1,s,1)+ 
			deltap*deltad*ddot(n,ds,1,dw,1)
		mu=acomp/dfloat(2*n)
		mua=g/dfloat(2*n)
		mu=mu*(mua/mu)**3
		#if(acomp>1) mu=(g/dfloat(n))*(g/acomp)**2
		#else mu=acomp/(dfloat(n)**2)
		do i=1,n{
			dz(i)=d(i)*(mu*(1/s(i)-1/x(i))+
				dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i))
			}
		call dswap(p,rhs,1,dy,1)
		call dgemv('N',p,n,one,a,p,dz,1,one,dy,1)#new rhs
		call dpotrs('U',p,1,ada,p,dy,p,info)
		call daxpy(p,mone,dy,1,rhs,1)#rhs=ddy
		call dgemv('T',p,n,one,a,p,rhs,1,zero,dw,1)#dw=A'ddy
		deltap=big
		deltad=big
		do i=1,n{
			dx(i)=dx(i)-dz(i)-d(i)*dw(i)
			ds(i)=-dx(i)
			dz(i)=mu/x(i) - z(i)*dx(i)/x(i) - z(i) - dxdz(i)/x(i)
			dw(i)=mu/s(i) - w(i)*ds(i)/s(i) - w(i) - dsdw(i)/s(i)
			if(dx(i)<0)deltap=dmin1(deltap,-x(i)/dx(i))
			else deltap=dmin1(deltap,-s(i)/ds(i))
			if(dz(i)<0)deltad=dmin1(deltad,-z(i)/dz(i))
			if(dw(i)<0)deltad=dmin1(deltad,-w(i)/dw(i))
			}
		deltap=dmin1(beta*deltap,one)
		deltad=dmin1(beta*deltad,one)
		}
	call daxpy(n,deltap,dx,1,x,1)
	call daxpy(n,deltap,ds,1,s,1)
	call daxpy(p,deltad,dy,1,y,1)
	call daxpy(n,deltad,dz,1,z,1)
	call daxpy(n,deltad,dw,1,w,1)
	cx=ddot(n,c,1,x,1)
	by=ddot(p,b,1,y,1)
	uw = dasum(n,w,1)
	uz = dasum(n,z,1)
	#rdg=(cx-by+uw)/(one+dabs(by-uw))
	#rdg=(cx-by+uw)/(one+uz+uw)
	rdg=(cx-by+uw)
	}
# return residuals in the vector x
call daxpy(n,mone,w,1,z,1)
call dswap(n,z,1,x,1)
return
end
