subroutine rqfnc(n1,n2,p,a1,y,a2,r,rhs,d1,d2,u,beta,eps,wn1,wn2,wp,nit,info)
integer n1,n2,p,info,nit(3)
double precision a1(p,n1),a2(p,n2),y(n1),r(n2),rhs(p),d1(n1),d2(n2),u(n1)
double precision wn1(n1,9),wn2(n2,6),wp(p,p+3)
double precision one,beta,eps
parameter(one = 1.0d0)
call lpfnc(n1,n2,p,a1,y,a2,r,rhs,d1,d2,u,beta,eps,wn1(1,1),wn2(1,1),wn1(1,2),
	wp(1,1),wn1(1,3),wn2(1,2),wn1(1,4),wn1(1,5),wn2(1,3),wn1(1,6),
        wp(1,2),wn1(1,7),wn2(1,4),wn1(1,8),wn1(1,9),wn2(1,5),wn2(1,6),
	wp(1,3),wp(1,4),nit,info)
return
end
# This is a revised form of the primal-dual log barrier form of the
# interior point LP solver based on Lustig, Marsten and Shanno ORSA J Opt 1992.
# It generalizes rqfnb and lpfnb to accomodate inequality constraints 
# It is a projected Newton primal-dual logarithmic barrier method which uses
# the predictor-corrector approach of Mehrotra for the mu steps.
# For the sake of brevity we will call it a Frisch-Newton algorithm.
# Problem:
# 	min c1'x1 + c2'x2  s.t. A1x1 + A2x2 = b, 0<=x1<=u, 0<=x2
# 
# Denote dx1,dx2,dy,dw,ds,dz1,dz2 as the steps for x1,x2,y,w,s,z1,z2
# 
subroutine lpfnc(n1,n2,p,a1,c1,a2,c2,b,d1,d2,u,beta,eps,x1,x2,s,
	y,z1,z2,w,dx1,dx2,ds,dy,dz1,dz2,dw,dr1,dr2,r2, rhs,ada,nit,info)

integer n1,p,i,info,nit(3),maxit
double precision a1(p,n1),a2(p,n2),c1(n1),c2(n2),b(p)
double precision zero,one,mone,big,ddot,dmax1,dmin1,dxdz1,dxdz2,dsdw
double precision deltap,deltad,beta,eps,mu,gap,g
double precision x1(n1),x2(n2),u(n1),s(n1),y(p),z1(n1),z2(n2),w(n1)
double precision d1(n1),d2(n2),rhs(p),ada(p,p)
double precision dx1(n1),dx2(n2),ds(n1),dy(p),dz1(n1),dz2(n2),dw(n1)
double precision dr1(n1),dr2(n2),r2(n2)

parameter(zero  = 0.0d0)
parameter(one   = 1.0d0)
parameter(mone   = -1.0d0)
parameter(big  = 1.0d+20)
parameter(maxit  = 500)

# Initialization:  We try to follow the notation of LMS
# On input we require:
# 
# 	c1 = n1-vector of marginal costs (-y in the rq problem)
# 	a1 = p by n1 matrix of equality constraints (x' in rq)
# 	c2 = n2-vector of marginal costs (-r in the rq problem)
# 	a2 = p by n2 matrix of inequality constraints (R' in rq)
# 	b = p-vector of rhs ((1-tau)x'e in rq)
#       u = upper bound vector ( 1_p in rq)
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
nit(3)=n1
# Start at the OLS estimate for the dual vector y
call dgemv('N',p,n1,one,a1,p,c1,1,zero,y,1) 
do i=1,n1
	d1(i)=one  
do i=1,n2{
	d2(i)=zero
	z2(i)=one
	}
call stepy2(n1,n2,p,a1,d1,a2,d2,y,ada,info)
if(info != 0) return
# put current residual vector r1 in s (temporarily)
call dcopy(n1,c1,1,s,1)
call dgemv('T',p,n1,mone,a1,p,y,1,one,s,1)
# Initialize remaining variables
do i=1,n1{
	if(dabs(s(i)) < eps){
		z1(i)=dmax1(s(i),zero)+eps
		w(i)=dmax1(-s(i),zero)+eps
		}
	else {
		z1(i)=dmax1(s(i),zero)
		w(i)=dmax1(-s(i),zero)
		}
	s(i)=u(i)-x1(i)
	}
gap = ddot(n1,z1,1,x1,1)+ddot(n2,z2,1,x2,1)+ddot(n1,w,1,s,1)
# 
# Main Loop
#
while(gap > eps && nit(1)<maxit) {
	nit(1)=nit(1)+1
	call dcopy(n2,c2,1,r2,1)
	call dgemv('T',p,n2,mone,a2,p,y,1,one,r2,1) 
	call dcopy(p,b,1,dy,1)
	call dgemv('N',p,n1,mone,a1,p,x1,1,one,dy,1) 
	call dgemv('N',p,n2,mone,a2,p,x2,1,one,dy,1) 
	do i = 1,n1{
		d1(i)=one/(z1(i)/x1(i) + w(i)/s(i)) 
		ds(i)=z1(i)-w(i) 
		dz1(i)=d1(i)*ds(i) 
		}
	do i = 1,n2{
		d2(i)=x2(i)/z2(i) 
		dz2(i)=d2(i)*r2(i) 
		}
	call dgemv('N',p,n1,one,a1,p,dz1,1,one,dy,1) 
	call dgemv('N',p,n2,one,a2,p,dz2,1,one,dy,1) 
	call dcopy(p,dy,1,rhs,1) #save rhs
	call stepy2(n1,n2,p,a1,d1,a2,d2,dy,ada,info)
	if(info != 0) return
	call dgemv('T',p,n1,one,a1,p,dy,1,mone,ds,1) 
	deltap=big
	deltad=big
	do i=1,n1{
		dx1(i)=d1(i)*ds(i)
		ds(i)=-dx1(i)
		dz1(i)=-z1(i)*(dx1(i)/x1(i) + one)
		dw(i)=-w(i)*(ds(i)/s(i) + one)
		if(dx1(i)<0)deltap=dmin1(deltap,-x1(i)/dx1(i))
		if(ds(i)<0)deltap=dmin1(deltap,-s(i)/ds(i))
		if(dz1(i)<0)deltad=dmin1(deltad,-z1(i)/dz1(i))
		if(dw(i)<0)deltad=dmin1(deltad,-w(i)/dw(i))
		}
	call dcopy(n2,r2,1,dx2,1)
	call dgemv('T',p,n2,one,a2,p,dy,1,mone,dx2,1) #dx2 = A2'dy - r2
	do i=1,n2{
		dx2(i)=d2(i)*dx2(i)
		dz2(i)=-z2(i)*(dx2(i)/x2(i) + one)
		if(dx2(i)<0)deltap=dmin1(deltap,-x2(i)/dx2(i))
		if(dz2(i)<0)deltad=dmin1(deltad,-z2(i)/dz2(i))
		}
	deltap=dmin1(beta*deltap,one)
	deltad=dmin1(beta*deltad,one)
	if(min(deltap,deltad) < one){
		nit(2)=nit(2)+1
		# Update mu
		mu = ddot(n1,x1,1,z1,1)+ddot(n2,x2,1,z2,1)+ddot(n1,s,1,w,1)
		g = mu + deltap*ddot(n1,dx1,1,z1,1)+
			deltad*ddot(n1,dz1,1,x1,1)+
			deltap*deltad*ddot(n1,dz1,1,dx1,1)+
			deltap*ddot(n2,dx2,1,z2,1)+
			deltad*ddot(n2,dz2,1,x2,1)+
			deltap*deltad*ddot(n2,dz2,1,dx2,1)+
			deltap*ddot(n1,ds,1,w,1)+
			deltad*ddot(n1,dw,1,s,1) +
			deltap*deltad*ddot(n1,ds,1,dw,1)
		mu = mu * ((g/mu)**3) /(dfloat(2*n1)+dfloat(n2))
		# Compute modified step
		do i=1,n1{
			dsdw =  ds(i)*dw(i)
			dr1(i)=d1(i)*(mu*(one/s(i)-one/x1(i))+
				dx1(i)*dz1(i)/x1(i)-dsdw/s(i))
			} 
		do i=1,n2{
			dr2(i)=d2(i)*(dx2(i)*dz2(i)/x2(i)-mu/x2(i))
			} 
		call dswap(p,rhs,1,dy,1)
		call dgemv('N',p,n1,one,a1,p,dr1,1,one,dy,1)
		call dgemv('N',p,n2,one,a2,p,dr2,1,one,dy,1)
		call dpotrs('U',p,1,ada,p,dy,p,info)
		call dgemv('T',p,n1,one,a1,p,dy,1,zero,u,1)
		deltap=big
		deltad=big
		do i=1,n1{
			dsdw =  ds(i)*dw(i)
			dxdz1 =  dx1(i)*dz1(i)
			dx1(i) =  d1(i)*(u(i)-z1(i)+w(i))-dr1(i)
			ds(i) = -dx1(i)
			dz1(i) = -z1(i)+(mu - z1(i)*dx1(i) - dxdz1)/x1(i)
			dw(i) = -w(i)+(mu - w(i)*ds(i) - dsdw)/s(i)
			if(dx1(i)<0)deltap=dmin1(deltap,-x1(i)/dx1(i))
			if(ds(i)<0)deltap=dmin1(deltap,-s(i)/ds(i))
			if(dz1(i)<0)deltad=dmin1(deltad,-z1(i)/dz1(i))
			if(dw(i)<0)deltad=dmin1(deltad,-w(i)/dw(i))
			}
		call dgemv('T',p,n2,one,a2,p,dy,1,zero,u,1)#u=A2'dy
		do i=1,n2{
			dxdz2 =  dx2(i)*dz2(i)
			dx2(i) =  d2(i)*(u(i)-r2(i))-dr2(i)
			dz2(i) = -z2(i)+(mu - z2(i)*dx2(i) - dxdz2)/x2(i)
			if(dx2(i)<0)deltap=dmin1(deltap,-x2(i)/dx2(i))
			if(dz2(i)<0)deltad=dmin1(deltad,-z2(i)/dz2(i))
			}
		deltap=dmin1(beta*deltap,one)
		deltad=dmin1(beta*deltad,one)
		}
	call daxpy(n1,deltap,dx1,1,x1,1)
	call daxpy(n2,deltap,dx2,1,x2,1)
	call daxpy(n1,deltap,ds,1,s,1)
	call daxpy(p,deltad,dy,1,y,1)
	call daxpy(n1,deltad,dz1,1,z1,1)
	call daxpy(n2,deltad,dz2,1,z2,1)
	call daxpy(n1,deltad,dw,1,w,1)
	gap = ddot(n1,z1,1,x1,1)+ddot(n2,z2,1,x2,1)+ddot(n1,w,1,s,1)
	}
# return residuals in the vector x1
call daxpy(n1,mone,w,1,z1,1)
call dswap(n1,z1,1,x1,1)
return
end
subroutine stepy2(n1,n2,p,a1,d1,a2,d2,b,ada,info)
integer n1,n2,p,i,j,k,info
double precision a1(p,n1),a2(p,n2),b(p),d1(n1),d2(n2),ada(p,p),zero
parameter(zero = 0.0d0)
# Solve the linear system ada'x=b by Choleski -- d is diagonal
# Note that a isn't altered, and on output ada returns the upper
# triangle Choleski factor, which can be reused, eg with blas dtrtrs
do j=1,p
	do k=1,p
		ada(j,k)=zero
do i=1,n1
	call dsyr('U',p,d1(i),a1(1,i),1,ada,p)
do i=1,n2
	call dsyr('U',p,d2(i),a2(1,i),1,ada,p)
call dposv('U',p,1,ada,p,b,p,info)
return
end
