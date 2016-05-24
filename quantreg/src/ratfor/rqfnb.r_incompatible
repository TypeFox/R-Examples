subroutine rqfnb(n,p,a,y,rhs,d,u,beta,eps,wn,wp,nit,info)
integer n,p,info,nit(3)
double precision a(p,n),y(n),rhs(p),d(n),u(n),wn(n,9),wp(p,p+3)
double precision one,beta,eps
parameter( one = 1.0d0)
call lpfnb(n,p,a,y,rhs,d,u,beta,eps,wn(1,1),wn(1,2),
        wp(1,1),wn(1,3),wn(1,4),wn(1,5),wn(1,6),
        wp(1,2),wn(1,7),wn(1,8),wn(1,9),wp(1,3),wp(1,4),nit,info)
return
end
# This is a revised form of my primal-dual log barrier form of the
# interior point LP solver based on Lustig, Marsten and Shanno ORSA J Opt 1992.
# It is a projected Newton primal-dual logarithmic barrier method which uses
# the predictor-corrector approach of Mehrotra for the mu steps.
# For the sake of brevity we will call it a Frisch-Newton algorithm.
# Problem:
# 	min c'x s.t. Ax=b, 0<=x<=u
# 
# Denote dx,dy,dw,ds,dz as the steps for the respective variables x,y,w,s,z
# 
subroutine lpfnb(n,p,a,c,b,d,u,beta,eps,x,s,y,z,w,
	dx,ds,dy,dz,dw,dr,rhs,ada,nit,info)

integer n,p,pp,i,info,nit(3),maxit
double precision a(p,n),c(n),b(p)
double precision zero,one,mone,big,ddot,dmax1,dmin1,dxdz,dsdw
double precision deltap,deltad,beta,eps,mu,gap,g
double precision x(n),u(n),s(n),y(p),z(n),w(n),d(n),rhs(p),ada(p,p)
double precision dx(n),ds(n),dy(p),dz(n),dw(n),dr(n)

parameter( zero  = 0.0d0)
parameter( one   = 1.0d0)
parameter( mone  = -1.0d0)
parameter( big  = 1.0d+20)
parameter( maxit  = 50)

# Initialization:  We follow the notation of LMS
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
# Start at the OLS estimate for the dual vector y
call dgemv('N',p,n,one,a,p,c,1,zero,y,1)
do i=1,n
	d(i)=one
call stepy(n,p,a,d,y,ada,info)
if(info != 0) return
# put current residual vector in s (temporarily)
call dcopy(n,c,1,s,1)
call dgemv('T',p,n,mone,a,p,y,1,one,s,1)
# Initialize remaining variables
do i=1,n{
	if(dabs(s(i))<eps){
		z(i)=dmax1(s(i), zero) + eps
		w(i)=dmax1(-s(i),zero) + eps
		}
	else {
		z(i)=dmax1(s(i), zero) 
		w(i)=dmax1(-s(i),zero)
		}
	s(i)=u(i)-x(i)
	}
gap = ddot(n,z,1,x,1)+ddot(n,w,1,s,1)
while(gap > eps && nit(1)<maxit) {
	nit(1)=nit(1)+1
	do i = 1,n{
		d(i) = one/(z(i)/x(i) + w(i)/s(i))
		ds(i)=z(i)-w(i)
		dz(i)=d(i)*ds(i)
		}
	call dcopy(p,b,1,dy,1)#save rhs
	call dgemv('N',p,n,mone,a,p,x,1,one,dy,1) 
	call dgemv('N',p,n,one,a,p,dz,1,one,dy,1) 
	call dcopy(p,dy,1,rhs,1)#save rhs
	call stepy(n,p,a,d,dy,ada,info)
	if(info != 0) return
	call dgemv('T',p,n,one,a,p,dy,1,mone,ds,1) #ds -> A'dy - ds
	deltap=big
	deltad=big
	do i=1,n{
		dx(i)=d(i)*ds(i)
		ds(i)=-dx(i)
		dz(i)=-z(i)*(dx(i)/x(i) + one)
		dw(i)=-w(i)*(ds(i)/s(i) + one)
		if(dx(i)<0)deltap=dmin1(deltap,-x(i)/dx(i))
		if(ds(i)<0)deltap=dmin1(deltap,-s(i)/ds(i))
		if(dz(i)<0)deltad=dmin1(deltad,-z(i)/dz(i))
		if(dw(i)<0)deltad=dmin1(deltad,-w(i)/dw(i))
		}
	deltap=dmin1(beta*deltap,one)
	deltad=dmin1(beta*deltad,one)
	if(min(deltap,deltad) < one){
		nit(2)=nit(2)+1
		# Update mu
		mu = ddot(n,x,1,z,1)+ddot(n,s,1,w,1)
		g = mu + deltap*ddot(n,dx,1,z,1)+
			deltad*ddot(n,dz,1,x,1) +
			deltap*deltad*ddot(n,dz,1,dx,1)+
			deltap*ddot(n,ds,1,w,1)+
			deltad*ddot(n,dw,1,s,1) +
			deltap*deltad*ddot(n,ds,1,dw,1)
		mu = mu * ((g/mu)**3) /dfloat(2*n)
		# Compute modified step
		do i=1,n{
			dr(i)=d(i)*(mu*(1/s(i)-1/x(i))+
				dx(i)*dz(i)/x(i)-ds(i)*dw(i)/s(i))
			}
		call dswap(p,rhs,1,dy,1)
		call dgemv('N',p,n,one,a,p,dr,1,one,dy,1)# new rhs
		call dpotrs('U',p,1,ada,p,dy,p,info)# backsolve for dy
		call dgemv('T',p,n,one,a,p,dy,1,zero,u,1)#ds=A'ddy
		deltap=big
		deltad=big
		do i=1,n{
			dxdz =  dx(i)*dz(i)
			dsdw =  ds(i)*dw(i)
			dx(i)=  d(i)*(u(i)-z(i)+w(i))-dr(i)
			ds(i)= -dx(i)
			dz(i)= -z(i)+(mu - z(i)*dx(i) - dxdz)/x(i)
			dw(i)= -w(i)+(mu - w(i)*ds(i) - dsdw)/s(i)
			if(dx(i)<0)deltap=dmin1(deltap,-x(i)/dx(i))
			if(ds(i)<0)deltap=dmin1(deltap,-s(i)/ds(i))
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
	gap = ddot(n,z,1,x,1)+ddot(n,w,1,s,1)
	}
# return residuals in the vector x
call daxpy(n,mone,w,1,z,1)
call dswap(n,z,1,x,1)
return
end
subroutine stepy(n,p,a,d,b,ada,info)
integer n,p,pp,i,info
double precision a(p,n),b(p),d(n),ada(p,p),zero
parameter( zero = 0.0d0)
# Solve the linear system ada'x=b by Choleski -- d is diagonal
# Note that a isn't altered, and on output ada returns the upper
# triangle Choleski factor, which can be reused, eg with blas dtrtrs
pp=p*p
do j=1,p
	do k=1,p
		ada(j,k)=zero
do i=1,n
	call dsyr('U',p,d(i),a(1,i),1,ada,p)
call dposv('U',p,1,ada,p,b,p,info)
return
end
