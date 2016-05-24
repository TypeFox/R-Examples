subroutine sspl(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,
		dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,
		xwy,hs,sg,abd,p1ip,ier)
implicit double precision(a-h,o-z) 

# A Cubic B-spline Smoothing routine.
#
#          The algorithm minimises:
#
#      (1/n) * sum w(i)* (y(i)-s(i))**2 + lambda* int ( s"(x) )**2 dx
#
#	for each of p response variables in y
#Input

#  x(n)		vector containing the ordinates of the observations
#  y(ldy,p)	matrix (n x p) of responses (ldy can be greater than n)
#  w(n)		vector containing the weights given to each data point
#  n    	number of data points
#  ldy		leading dimension of y
#  p		number of columns in y	
#  knot(nk+4)	vector of knot points defining the cubic b-spline basis.
#  nk		number of b-spline coefficients to be estimated
#			nk <= n+2
#  method 	method for selecting amount of smoothing, lambda
#		1 = fixed lambda 
#		2 = fixed df
#		3 = gcv
#		4 = cv
#  tol		used in Golden Search routine
#  wp(p)	weights, length p,  used to combine cv or gcv in 3 or 4 above
#  ssy(p)	offsets for weighted sum of squares for y; can be all zero,
#			else should be the variability lost due to collapsing
#			onto unique values
#  dfoff	offset df used in gcv calculations (0 is good default)
#  dfmax	maximum value for df allowed when gcv or cv are used
#  		routine simply returns the value at dfmax if it was exceeded
#  cost		cost per df (1 is good default)
#Input/Output
#  lambda	penalised likelihood smoothing parameter
#  df		trace(S)
#Output
#  cv		omnibus cv criterion
#  gcv		omnibus gcv criterion (including penalty and offset)
#  coef(nk,p)	vector of spline coefficients
#  s(ldy,p)	matrix of smoothed y-values
#  lev(n)	vector of leverages
# Working arrays/matrix
#  xwy(nk,p)	X'Wy
#  hs(nk,4)	the diagonals of the X'WX matrix
#  sg(nk,4)   	the diagonals of the Gram matrix
#  abd(4,nk)	[ X'WX+lambda*SIGMA] in diagonal form
#  p1ip(4,nk)	inner products between columns of L inverse
#  ier          error indicator
#                  ier = 0 ___  everything fine
#                  ier = 1 ___  spar too small or too big
#                               problem in cholesky decomposition
double precision x(n),y(ldy,p),w(n),knot(nk+4),tol,wp(p),ssy(p),
	 	dfoff,dfmax,cost,lambda,df,cv,gcv,coef(nk,p),s(ldy,p),lev(n),
		xwy(nk,p),hs(nk,4),sg(nk,4),abd(4,nk),p1ip(4,nk)
integer         n,p,ldy,nk,method,ier

#  Compute SIGMA, X'WX, X'WY, trace, ratio, s0, s1.

# SIGMA-> sg[]
# X'WX -> hs[]
# X'WY -> xwy[]
call sgram(sg(1,1),sg(1,2),sg(1,3),sg(1,4),knot,nk)
call stxwx2(x,y,w,n,ldy,p,knot,nk,xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,4))

# Compute estimate

if(method==1) {# Value of lambda supplied

	call sslvr2(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,
		dfoff,cost,lambda,df,cv,gcv,coef,s,lev,
		xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,4),
		sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
	}

else {# Use Forsythe, Malcom and Moler routine to minimise criterion
	call fmm(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,
		dfoff,cost,lambda,df,cv,gcv,coef,s,lev,
		xwy,hs,sg,abd,p1ip,ier)
	if(method>2&df>dfmax){
		df=dfmax
		call fmm(x,y,w,n,ldy,p,knot,nk,2,tol,wp,ssy,
		dfoff,cost,lambda,df,cv,gcv,coef,s,lev,
		xwy,hs,sg,abd,p1ip,ier)
		}
	}
return
end

subroutine fmm(xs,ys,ws,n,ldy,nvar,knot,nk,method,tol,wp,ssy,
		dfoff,cost,lambda,df,cv,gcv,coef,s,lev,
		xwy,hs,sg,abd,p1ip,ier)
double precision xs(n),ys(ldy,nvar),ws(n),knot(nk+4),tol,wp(nvar),ssy(nvar),
	 	dfoff,cost,lambda,df,cv,gcv,coef(nk,nvar),s(ldy,nvar),lev(n),
		xwy(nk,nvar),hs(nk,4),sg(nk,4),abd(4,nk),p1ip(4,nk)
integer         n,ldy,nvar,nk,method,ier

# Local variables

double precision  t1,t2,ratio,
                  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w,
                  fu,fv,fw,fx,x,targdf,
		  ax,bx

integer i
ax=1d-10  	#used to be lspar
bx=1.5		#used to be uspar
t1=0. ; t2=0.
targdf=df
do i=3,nk-3 { t1 = t1 + hs(i,1) }
do i=3,nk-3 { t2 = t2 + sg(i,1) }
ratio = t1/t2
#
#  an approximation  x  to the point where f  attains a minimum  on
#  the interval  (ax,bx)  is determined.
#
#
#  input..
#
#  ax	 left endpoint of initial interval
#  bx	 right endpoint of initial interval
#  f	 function subprogram which evaluates  f(x)  for any  x
#	 in the interval  (ax,bx)
#  tol	 desired length of the interval of uncertainty of the final
#	 result ( .ge. 0.0)
#
#
#  output..
#
#  fmin  abcissa approximating the point where	f  attains a minimum
#
#
#      the method used is a combination of  golden  section  search  and
#  successive parabolic interpolation.	convergence is never much slower
#  than  that  for  a  fibonacci search.  if  f  has a continuous second
#  derivative which is positive at the minimum (which is not  at  ax  or
#  bx),  then  convergence  is	superlinear, and usually of the order of
#  about  1.324....
#      the function  f	is never evaluated at two points closer together
#  than  eps*dabs(fmin) + (tol/3), where eps is	approximately the square
#  root  of  the  relative  machine  precision.   if   f   is a unimodal
#  function and the computed values of	 f   are  always  unimodal  when
#  separated by at least  eps*dabs(x) + (tol/3), then  fmin  approximates
#  the abcissa of the global minimum of  f  on the interval  ax,bx  with
#  an error less than  3*eps*dabs(fmin) + tol.  if   f	is not unimodal,
#  then fmin may approximate a local, but perhaps non-global, minimum to
#  the same accuracy.
#      this function subprogram is a slightly modified	version  of  the
#  algol  60 procedure	localmin  given in richard brent, algorithms for
#  minimization without derivatives, prentice - hall, inc. (1973).
#
#
#      double precision  a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
#      double precision  fu,fv,fw,fx,x
#
#  c is the squared inverse of the golden ratio
#
      c = 0.5*(3. - dsqrt(5d0))
#
#  eps is approximately the square root of the relative machine
#  precision.
#
      eps = 1d0
   10 eps = eps/2d0
      tol1 = 1d0 + eps
      if (tol1 .gt. 1d0) go to 10
      eps = dsqrt(eps)
#
#  initialization
#
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0
	lambda = ratio*16.**(-2. + x*(6.))
	call sslvr2(xs,ys,ws,n,ldy,nvar,knot,nk,method,tol,wp,ssy,
		dfoff,cost,lambda,df,cv,gcv,coef,s,lev,
		xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,4),
		sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
switch(method){
	case 2:
		fx=3d0+(targdf-df)**2
	case 3:
		fx=gcv
	case 4:
		fx=cv
	}
      fv = fx
      fw = fx
#
#  main loop starts here
#
   20 xm = 0.5*(a + b)
      tol1 = eps*dabs(x) + tol/3d0
      tol2 = 2d0*tol1
#
#  check stopping criterion
#
      if (dabs(x - xm) .le. (tol2 - 0.5*(b - a))) go to 90
#
# is golden-section necessary
#
      if (dabs(e) .le. tol1) go to 40
#
#  fit parabola
#
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.00*(q - r)
      if (q .gt. 0.0) p = -p
      q =  dabs(q)
      r = e
      e = d
#
#  is parabola acceptable
#
   30 if (dabs(p) .ge. dabs(0.5*q*r)) go to 40
      if (p .le. q*(a - x)) go to 40
      if (p .ge. q*(b - x)) go to 40
#
#  a parabolic interpolation step
#
      d = p/q
      u = x + d
#
#  f must not be evaluated too close to ax or bx
#
      if ((u - a) .lt. tol2) d = dsign(tol1, xm - x)
      if ((b - u) .lt. tol2) d = dsign(tol1, xm - x)
      go to 50
#
#  a golden-section step
#
   40 if (x .ge. xm) e = a - x
      if (x .lt. xm) e = b - x
      d = c*e
#
#  f must not be evaluated too close to x
#
   50 if (dabs(d) .ge. tol1) u = x + d
      if (dabs(d) .lt. tol1) u = x + dsign(tol1, d)

		lambda = ratio*16.**(-2. + u*(6.))
	call sslvr2(xs,ys,ws,n,ldy,nvar,knot,nk,method,tol,wp,ssy,
		dfoff,cost,lambda,df,cv,gcv,coef,s,lev,
		xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,4),
		sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)
switch(method){
	case 2:
		fu=3d0+(targdf-df)**2
	case 3:
		fu=gcv
	case 4:
		fu=cv
	}
#
#  update  a, b, v, w, and x
#
      if (fu .gt. fx) go to 60
      if (u .ge. x) a = x
      if (u .lt. x) b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 if (u .lt. x) a = u
      if (u .ge. x) b = u
      if (fu .le. fw) go to 70
      if (w .eq. x) go to 70
      if (fu .le. fv) go to 80
      if (v .eq. x) go to 80
      if (v .eq. w) go to 80
      go to 20
   70 v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 v = u
      fv = fu
      go to 20
#
#  end of main loop
#
   90 continue 
if(method==2){call sslvr2(xs,ys,ws,n,ldy,nvar,knot,nk,1,tol,wp,ssy,
	dfoff,cost,lambda,df,cv,gcv,coef,s,lev,
	xwy,hs(1,1),hs(1,2),hs(1,3),hs(1,4),
	sg(1,1),sg(1,2),sg(1,3),sg(1,4),abd,p1ip,ier)}
	
return 
end
subroutine stxwx2(x,z,w,k,ldy,pz,xknot,n,y,hs0,hs1,hs2,hs3)
implicit double precision(a-h,o-z) 
double precision  z(ldy,pz),w(k),x(k),xknot(n+4),
		 y(n,pz),hs0(n),hs1(n),hs2(n),hs3(n),
		 eps,vnikx(4,1),work(16) # local
integer  k,n,pz,ldy,
j,i,pp,ileft,mflag  # local

# Initialise the output vectors

do i=1,n { 
	hs0(i)=0d0  
	hs1(i)=0d0
	hs2(i)=0d0
	hs3(i)=0d0
	do j=1,pz { y(i,j)=0d0 }
	}

# Compute X'WX -> hs0,hs1,hs2,hs3  and X'WZ -> y

eps = .1d-9
do i=1,k {
	call interv(xknot(1),(n+1),x(i),ileft,mflag)
	if(mflag== 1) {
		if(x(i)<=(xknot(ileft)+eps)){ileft=ileft-1}
		else{return}
	}
	call bsplvd (xknot,4,x(i),ileft,work,vnikx,1)

	j= ileft-4+1
	do pp=1,pz {y(j,pp) = y(j,pp)+w(i)*z(i,pp)*vnikx(1,1)}
	hs0(j)=hs0(j)+w(i)*vnikx(1,1)**2
	hs1(j)=hs1(j)+w(i)*vnikx(1,1)*vnikx(2,1)
	hs2(j)=hs2(j)+w(i)*vnikx(1,1)*vnikx(3,1)
	hs3(j)=hs3(j)+w(i)*vnikx(1,1)*vnikx(4,1)

	j= ileft-4+2
	do pp=1,pz {y(j,pp) =  y(j,pp)+w(i)*z(i,pp)*vnikx(2,1)}
	hs0(j)=hs0(j)+w(i)*vnikx(2,1)**2
	hs1(j)=hs1(j)+w(i)*vnikx(2,1)*vnikx(3,1)
	hs2(j)=hs2(j)+w(i)*vnikx(2,1)*vnikx(4,1)

	j= ileft-4+3
	do pp=1,pz {y(j,pp) =  y(j,pp)+w(i)*z(i,pp)*vnikx(3,1)}
	hs0(j)=hs0(j)+w(i)*vnikx(3,1)**2
	hs1(j)=hs1(j)+w(i)*vnikx(3,1)*vnikx(4,1)

	j= ileft-4+4
	do pp=1,pz {y(j,pp) =  y(j,pp)+w(i)*z(i,pp)*vnikx(4,1)}
	hs0(j)=hs0(j)+w(i)*vnikx(4,1)**2  
	}

return
end
subroutine sslvr2(x,y,w,n,ldy,p,knot,nk,method,tol,wp,ssy,
		dfoff,cost,lambda,df,cv,gcv,coef,sz,lev,
		xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3,
		abd,p1ip,info)

implicit double precision(a-h,o-z) 
double precision x(n),y(ldy,p),w(n),knot(nk+4),tol,wp(p),ssy(p),
	 	dfoff,cost,lambda,df,cv,gcv,coef(nk,p),sz(ldy,p),lev(n),
		xwy(nk,p),
		hs0(nk),hs1(nk),hs2(nk),hs3(nk),
		sg0(nk),sg1(nk),sg2(nk),sg3(nk),
		abd(4,nk),p1ip(4,nk)
integer         n,p,ldy,nk,method,info
#local storage
double precision b0,b1,b2,b3,eps,vnikx(4,1),work(16),
	 xv,bvalue,rss,tssy
integer  ld4,i,icoef,ileft,ilo,j,mflag
logical fittoo
fittoo= (method!=2)

         ilo = 1 ; eps = .1d-10 ; ld4=4

# Purpose : Solves the smoothing problem and computes the
#           criterion functions (CV and GCV).
# The coeficients of estimated smooth

if(fittoo){
	do i=1,nk { do j=1,p {coef(i,j) = xwy(i,j) } }
	}
do i=1,nk { abd(4,i)   = hs0(i)+lambda*sg0(i) }
do i=1,(nk-1) { abd(3,i+1) = hs1(i)+lambda*sg1(i) }
do i=1,(nk-2) { abd(2,i+2) = hs2(i)+lambda*sg2(i) }
do i=1,(nk-3) { abd(1,i+3) = hs3(i)+lambda*sg3(i) }

call dpbfa(abd,ld4,nk,3,info)
if(info.ne.0) { return } 
if(fittoo){
	do j=1,p{ call dpbsl(abd,ld4,nk,3,coef(1,j)) }

# Value of smooths at the data points

	icoef = 1
	do i=1,n { 
		xv = x(i)
		do j=1,p{ sz(i,j) = bvalue(knot,coef(1,j),nk,4,xv,0) }
		}
	}
# Compute the criterion functions

# Get Leverages First

call sinrp2(abd,ld4,nk,p1ip)

do i=1,n { 
	xv = x(i)
	call interv(knot(1),(nk+1),xv,ileft,mflag)
	if(mflag==-1) { ileft = 4   ; xv = knot(4)+eps }
	if(mflag==1)  { ileft = nk  ; xv = knot(nk+1)-eps }
	j=ileft-3
	call bsplvd(knot,4,xv,ileft,work,vnikx,1)
	b0=vnikx(1,1);b1=vnikx(2,1);b2=vnikx(3,1);b3=vnikx(4,1)
	lev(i) = (p1ip(4,j)*b0**2   + 2.*p1ip(3,j)*b0*b1 +
				   2.*p1ip(2,j)*b0*b2 + 2.*p1ip(1,j)*b0*b3 +
		  p1ip(4,j+1)*b1**2 + 2.*p1ip(3,j+1)*b1*b2 +
				   2.*p1ip(2,j+1)*b1*b3 +
		  p1ip(4,j+2)*b2**2 + 2.*p1ip(3,j+2)*b2*b3 +
		  p1ip(4,j+3)*b3**2 )*w(i)    
	}
# Evaluate Criteria
rss = 0d0 ; df = 0d0 ; sumw=0d0;gcv=0d0;cv=0d0;
do i=1,n { 
	df  = df  + lev(i)
	}
if(fittoo){
	do i=1,n { 
		sumw  = sumw  + w(i)
		do j=1,p{
			rss = rss + w(i)*wp(j)*(y(i,j)-sz(i,j))**2
			cv = cv +w(i)*wp(j)*((y(i,j)-sz(i,j))/(1-lev(i)))**2
			} 
		}
tssy=0d0
do j=1,p{tssy=tssy+wp(j)*ssy(j)}

gcv=((rss+tssy)/sumw)/((1d0-((dfoff+df-1)*cost+1)/sumw)**2)

#note: the weights should sum to n (the number of original observations)

cv=(cv+tssy)/sumw
	}
#lev includes the weights
#Note that this version of cv omits ALL observations at
#tied x values, since the data are already collapsed here
return 
end
subroutine sinrp2(abd,ld4,nk,p1ip)
implicit double precision(a-h,o-z) 
double precision abd(ld4,nk),p1ip(ld4,nk),
	wjm3(3),wjm2(2),wjm1(1),c0,c1,c2,c3

integer	ld4,nk,i,j

	# Purpose :  Computes Inner Products between columns of L(-1)
	#	     where L = abd is a Banded Matrix with 3 subdiagonals

	#		A refinement of Elden's trick is used.
	#	Coded by Finbarr O'Sullivan
wjm3(1)=0d0; wjm3(2)=0d0; wjm3(1)=0d0
wjm2(1)=0d0; wjm2(2)=0d0
wjm1(1)=0d0

do i=1,nk { 
	j=nk-i+1
	c0 = 1d0/abd(4,j)
	if(j<=nk-3) {
		c1 = abd(1,j+3)*c0
		c2 = abd(2,j+2)*c0
		c3 = abd(3,j+1)*c0 
		}
	else if(j==nk-2) {
		c1 = 0d0
		c2 = abd(2,j+2)*c0
		c3 = abd(3,j+1)*c0 
		}
	else if(j==nk-1) {
		c1 = 0d0
		c2 = 0d0
		c3 = abd(3,j+1)*c0 
		}
	else if(j==nk) {
		c1 = 0d0
		c2 = 0d0
		c3 = 0d0
	}
	p1ip(1,j) = 0d0- (c1*wjm3(1)+c2*wjm3(2)+c3*wjm3(3))
	p1ip(2,j) = 0d0- (c1*wjm3(2)+c2*wjm2(1)+c3*wjm2(2))
	p1ip(3,j) = 0d0- (c1*wjm3(3)+c2*wjm2(2)+c3*wjm1(1))

	p1ip(4,j) = c0**2 + 
		c1**2*wjm3(1)+2.*c1*c2*wjm3(2)+2.*c1*c3*wjm3(3) +
		c2**2*wjm2(1)+2.*c2*c3*wjm2(2) +
		c3**2*wjm1(1)

	wjm3(1)=wjm2(1) ; wjm3(2)=wjm2(2) ; wjm3(3)=p1ip(2,j)
	wjm2(1)=wjm1(1) ; wjm2(2)=p1ip(3,j); wjm1(1)=p1ip(4,j)
	}
return
end
subroutine suff2(n,p,ny,match,y,w,ybar,wbar,ssy,work)
integer match(n),n,ny,p,i,j
double precision y(n,ny),ybar(p+1,ny),w(n),wbar(p+1),ssy(ny),work(n)
double precision tsum
#ssy is the within response variability that is lost by collapsing
call pack(n,p,match,w,wbar)
do j=1,ny{
	do i=1,n
		work(i)=y(i,j)*w(i)
	call pack(n,p,match,work,ybar(1,j))
	do i=1,p{
		if(wbar(i)>0d0) 
			ybar(i,j)=ybar(i,j)/wbar(i) 
		else ybar(i,j)=0d0
		}
	call unpack(n,p,match,ybar(1,j),work)
	tsum=0d0
	do i=1,n
		tsum=tsum+ w(i)*(y(i,j)-work(i))**2
	ssy(j)=tsum
	}
return
end
subroutine sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,
		dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,
		xrange,work,ier)
double precision x(n),y(n,p),w(n),knot(nk+4),tol,wp(p),
	 	dfoff,dfmax,cost,lambda,df,cv,gcv,coef(1),s(n,p),lev(nef),
		xrange(2),work(1)
#workspace must be (2*p+2)*nefp1 + (p+16)*nk + n +p
integer         n,p,nk,method,ier,nef, nefp1, n2,match(1)
logical center
double precision  xmiss,sigtol,xdiff,temp
if(nef==0){# match has not been initialized
	xmiss=1d20
	sigtol=1d-5
	call namat(x,match,n,nef,work,work(n+1),xmiss,sigtol)
	xrange(1)=work(1) #work is actually the sorted unique xs
	xrange(2)=work(nef)
	}
else{
	do i=1,n {work(match(i))=x(i)}
	}
xdiff=xrange(2)-xrange(1)
do i=1,nef {work(i)=(work(i)-xrange(1))/xdiff}
if(nk==0){
	call sknotl(work,nef,knot,nk)
	nk=nk-4
	}
if(dfmax > dble(nk))dfmax=dble(nk)
if(cost>0){
	temp=dble(n-dble(center))/cost - dfoff
	if(dfmax>temp)dfmax=temp
	}
nefp1=nef+1
n2=nefp1*(2*p+2)+1
call sspl1(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,nefp1,center,
	dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,
	work(1),work(nefp1+1), 		#xin,yin
	work(nefp1*(p+1)+1),work(nefp1*(p+2)+1), 	#win, sout
	work(n2),work(n2+p*nk),work(n2+(p+4)*nk), 	#xwy, hs,sg
	work(n2+(p+8)*nk),work(n2+(p+12)*nk),work(n2+(p+16)*nk),
	work(n2+(p+16)*nk+p),ier)

return
end

#Memory management subroutine
subroutine sspl1(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,nefp1,center,
		dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,lev,
		xin,yin,win,sout,
		xwy,hs,sg,abd,p1ip,ssy,work,ier)
double precision x(n),y(n,p),w(n),knot(nk+4),tol,wp(p),
	 	dfoff,dfmax,cost,lambda,df,cv,gcv,coef(nk,p),s(n,p),lev(nef),
		xin(nefp1),yin(nefp1,p),win(nefp1),sout(nefp1,p),
		xwy(nk,p),hs(nk,4),sg(nk,4),abd(4,nk),p1ip(4,nk),
		ssy(p),work(n)
integer         n,p,nefp1,nk,method,ier,match(n),nef
logical center
double precision  sbar, wmean,tdfoff
call suff2(n,nef,p,match,y,w,yin,win,ssy,work)
if(center){
	if(cost>0){tdfoff=dfoff-1/cost}
}
call sspl(xin,yin,win,nef,nefp1,p,knot,nk,method,tol,wp,ssy,
		tdfoff,dfmax,cost,lambda,df,cv,gcv,coef,sout,lev,
		xwy,hs,sg,abd,p1ip,ier)
	
#now unpack the results
do j=1,p{
	call unpack(n,nef,match,sout(1,j),s(1,j))
	if(center){
		sbar=wmean(nef,sout(1,j),win)
		do i=1,n{ s(i,j)=s(i,j)-sbar}
	}
}
if(center)df=df-1
return
end
subroutine namat(x,match,n,nef,work,iwork,xmiss,tol)
#returns match (order) and work(1:nef) is the sorted unique x
implicit double precision(a-h,o-z)
integer match(1),n,nef,iwork(1),index
double precision x(1),xmiss,work(1),tol,xend,xstart
do i=1,n {
	work(i)=x(i)
	iwork(i)=i
	}

call sortdi(work,n,iwork,1,n)
xstart=x(iwork(1))
index=n
xend=x(iwork(n))
while(xend >= xmiss & index > 1){
	index=index-1
	xend=x(iwork(index))
	}

tol=tol*(xend-xstart)
index=1
work(1)=xstart
for(i=1;i<=n;i=i+1){
	while((x(iwork(i))-xstart)<tol){
		match(iwork(i))=index
		i=i+1
		if(i>n)goto 10
		}
	xstart= x(iwork(i))
	index=index+1
	match(iwork(i))=index
	work(index)=xstart
	}
10 if(xstart >= xmiss) 
	{nef=index-1}
	 else {nef=index}
return 
end

subroutine simfit(x,y,w,n,p,dfoff,cost,wp,gcv,coef,s,type,center,work)
#
# computes constant and linear fits, and selects the best using gcv
#
implicit double precision (a-h,o-z)

integer n,p,type
double precision x(n),y(n,p),w(n),cost,dfoff,wp(p),gcv,coef(2,p),
	s(n,p),work(p)
logical center
double precision sx,sy,sumw, dcent,sxx,syy,sxy
dcent=1-dble(center)
#center is F for no centering, else T
#Note: if type enters 1 or 2, no selection is made
sumw=0d0;gcvc=0d0;
do i=1,n {
  sumw=sumw+w(i)
}
if(type!=1){#either 0 or 2 in which case the linear is needed as well
	sx=0.0 ;  sxx=0.0; gcvl=0d0;

	do i=1,n {
  	sx=sx+w(i)*x(i)
	}
	xbar=sx/sumw
	do i=1,n {
	  sxx=sxx+w(i)*(x(i)-xbar)*x(i)
	}
}
do j=1,p{
	sy=0d0;syy=0d0;
	do i=1,n{
		sy=sy+w(i)*y(i,j)
	}
	work(j)=sy/sumw
  	do i=1,n{
		syy=syy+w(i)*(y(i,j)-work(j))*y(i,j)
	}
	gcvc=gcvc+wp(j)*syy
	if(type!=1){ #once again, do for linear as well
		sxy=0.0;
  		do i=1,n{
  			sxy=sxy+w(i)*(x(i)-xbar)*y(i,j)
		}
		coef(2,j)=sxy/sxx
		gcvl=gcvl+wp(j)*(syy -sxy*coef(2,j))
	}
}
if(type==0){
	gcvc =gcvc/ (sumw* (1 - (dfoff*cost + dcent)/sumw)**2 )
	gcvl=gcvl/(sumw* (1 - (dcent + (dfoff +1)* cost)/sumw)**2)
	if(gcvc<=gcvl){
		type=1
		gcv=gcvc
	}
	else{
		type=2
		gcv=gcvl
	}
}
else {

	if(type==1) {gcv=gcvc/(sumw* (1 - (dfoff*cost + dcent)/sumw)**2 )}
	else {gcv=gcvl/(sumw* (1 - (dcent + (dfoff + 1)*cost)/sumw)**2)}
}
	
if(type==1){
	do j=1,p{
		coef(1,j)=work(j)*dcent
		coef(2,j)=0d0
		do i=1,n {s(i,j)=coef(1,j)}
		}
	}
else{
	do j=1,p{
		coef(1,j)=work(j)*dcent-xbar*coef(2,j)
		do i=1,n {s(i,j)=coef(1,j)+coef(2,j)*x(i)}
		}
	}	
return
end

 
subroutine sspl2(x,y,w,n,p,knot,nk,wp,match,nef,
		dfoff,dfmax,cost,lambda,df,gcv,coef,s,type,center,
		xrange,work,tol,ier)
double precision x(n),y(n,p),w(n),knot(nk+4),wp(p),
	 	dfoff,dfmax,cost,lambda,df,gcv,coef(1),s(n,p),
		xrange(2),work(1),tol
#this routine selects from the linear and constant model as well
#see documentation for sspl
#workspace must be (2*p+2)*nefp1 + (p+16)*nk + 2*n
#if type>0 then no selection is performed; the fit is simply computed.
integer         n,p,nk,nef,type,match(n),ier,method
double precision coef1,coef2,cv
logical center
#center is F for no centering, else T
if(type==3){
	method=1
call sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,
		dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,work(1),
		xrange,work(n+1),ier)
	return
}
if(type>0){
       call simfit(x,y,w,n,p,dfoff,cost,wp,gcv,coef,s,
		type,center,work)
	df=dble(type)- dble(center)
return
}
#selection is being performed
method=3
call sspl0(x,y,w,n,p,knot,nk,method,tol,wp,match,nef,center,
		dfoff,dfmax,cost,lambda,df,cv,gcv,coef,s,work(1),
		xrange,work(n+1),ier)
gcv1=gcv
call simfit(x,y,w,n,p,dfoff,cost,wp,gcv,work,work(2*p+1),type,center,
		work((n+2)*p+1))
if(gcv<=gcv1){
#the coef swapping is so as not to destroy the spline coefs if needed
	df=dble(type)- dble(center)
	do j=1,p{
		coef1=work(1+(j-1)*2)
		coef2=work(2+(j-1)*2)
		if(type==1){
			do i=1,n {s(i,j)=coef1}
			}
		else{
			do i=1,n {s(i,j) =coef1+coef2*x(i)}
			}
		coef(1+(j-1)*2)=coef1
		coef(2+(j-1)*2)=coef2
		}
	}
else{
	type=3
	gcv=gcv1
	}
return
end
subroutine psspl2(x,n,p,knot,nk,xrange,coef,coefl,s,order,type)
implicit double precision(a-h,o-z) 
#make predictions from a fitted smoothing spline 
double precision x(n),knot(nk+4),xrange(2),coef(nk,p),coefl(2,1),s(n,p)
integer n,p,nk,order, type
double precision ytemp
switch(type){
	case 1:{
		do j=1,p{
			if(order>=1){ytemp=0d0} else {ytemp=coefl(1,j)}
			do i=1,n {s(i,j)=ytemp}
		}
	}
	case 2:{
		if(order>=1){
			do j=1,p{
				if(order==1){ytemp=coefl(2,j)}
					else {ytemp=0d0}
				do i =1,n {s(i,j)=ytemp}
			}
		}
		else{ 
			do j=1,p{
				do i=1,n {s(i,j)=coefl(1,j)+coefl(2,j)*x(i)}
			}
		}
	}
	case 3: {
		call psspl(x,n,p,knot,nk,xrange,coef,s,order)
	}
}
return
end

		

subroutine psspl(x,n,p,knot,nk,xrange,coef,s,order)
implicit double precision(a-h,o-z) 
#make predictions from a fitted smoothing spline, linear or constant
double precision x(n),knot(nk+4),xrange(2),coef(nk,p),s(n,p)
integer n,p,nk,order
double precision xcs,xmin,xdif, endv(2),ends(2),xends(2),stemp
double precision bvalue
integer where
if(order>2|order<0)return
xdif=xrange(2)-xrange(1)
xmin=xrange(1)
xends(1)=0d0
xends(2)=1d0
do j=1,p{
	if(order==0){
		endv(1)=bvalue(knot,coef(1,j),nk,4,0d0,0)
		endv(2)=bvalue(knot,coef(1,j),nk,4,1d0,0)
		}
	if(order<=1){
		ends(1)=bvalue(knot,coef(1,j),nk,4,0d0,1)
		ends(2)=bvalue(knot,coef(1,j),nk,4,1d0,1)
		}
	do i=1,n{
		xcs=(x(i)-xmin)/xdif
		where=0
		if(xcs<0d0){where=1}
		if(xcs>1d0){where=2}
		if(where>0){#beyond extreme knots
			switch(order){
				case 0:
					stemp=endv(where)+
						(xcs-xends(where))*ends(where)
				case 1:
					stemp=ends(where)
				case 2:
					stemp=0d0
			}
		}
		else {stemp=bvalue(knot,coef(1,j),nk,4,xcs,order)}
		if(order>0){s(i,j)=stemp/(xdif**dble(order))}
		else s(i,j)=stemp
	}
}
return
end
