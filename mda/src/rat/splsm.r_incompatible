subroutine splsm(x,y,w,n,match,nef,spar,dof,smo,s0,cov,ifcov,work,lenw)
#This subroutine performs a smoothing spline fit
#All arguments are either double precision or integer
#INPUT
#
#x	double precision length n ; x variable for smoothing
#y	double precision length n	; y variable for smoothing
#w	double precision length n ; weights for smoothing, > 0
#n	integer length above
#match 	integer length n -- in S language x[i] == sort(unique(x)[match[i]]
#		match is produced by subroutine namat
#nef	number of unique elements in x; so match has values between 1 and nef+1
#		missing data are given the match number nef+1
#spar	double precision smoothing parameter 0<spar<1.5
#		if spar is 0 and dof is 0, gcv is used
#		if spar>0, spar is used
#dof	double precision  equivalent degrees of freedom
#		if dof > 0 and spar = 0 then dof is used to 
#		select amount of smoothing
#		note: dof does not use the constant term
#ifcov	logical if true, the unscaled variance information is computed
#work	double precision workspace of length (10+2*4)*(nef+2)+5*nef+n+15
#
#OUTPUT
#
#x,y,w,n,match,nef  are untouched
#spar	if spar was 0 and dof was 0, then spar is that spar 
#			that minimized gcv
#	if spar was 0 and dof > 0, then spar is that which achieves dof
#dof	the dof of the fitted smooth. Note: even if dof was given
#		as 4, it will be returned as say 3.995 which is what
#		spar produces
#smo	double precision length n the fitted values, with weighted average 0
#s0	double precision weighted mean of y
#cov	double precision length nef the unscaled variance elements for the NONLINEAR
#		and UNIQUE part of smo, in the order of sort(unique(x))
#		cov is lev(i)/w(i) -h(i)/w where h(i) is the hat element from
#		the simple weighted least squares fit. This is passed on
#		to bakfit and used in gamcov
#
# splsm calls (eventually after some memory management dummy calls)
# sbart, the spline routine of Finbarr O'Sullivan, slightly modified
# by Trevor Hastie, 8/2/89	

implicit double precision(a-h,o-z)
integer n,match(n),nef,lenw
double precision x(n),y(n),w(n),spar,dof,smo(n),s0,cov(n),work(lenw)
logical ifcov
# work should be (10+2*ld4)*nk+5*nef+n+15 double precision
# ld4 =4  nk<= nef+2
call splsm1(x,y,w,n,match,nef,spar,dof,smo,s0,cov,ifcov,
#	xin(nef+1),yin(nef+1), win(nef+1),   knot(n+6),
	work(1),   work(nef+2),work(2*nef+3),work(3*nef+4),
	work(3*nef+n+10),lenw)
return
end

subroutine splsm1(x,y,w,n,match,nef,spar,dof,smo,s0,lev,ifcov,
		xin,yin,win,knot,
		work,lenw)
implicit double precision(a-h,o-z)
integer n,match(n),nef,lenw
double precision x(n),y(n),w(n),spar,dof,smo(n),s0,lev(n),work(lenw)
logical ifcov
double precision xin(nef+1),yin(nef+1),win(nef+1),knot(nef+4)
integer nk,ldnk,ld4,k
double precision xmin,xrange
call suff(n,nef,match,x,y,w,xin,yin,win,work(1))
xmin=xin(1)
xrange=xin(nef)-xin(1)
do i=1,nef {xin(i)=(xin(i)-xmin)/xrange}
call sknotl(xin,nef,knot,k)
nk=k-4
ld4=4
ldnk=1 # p21p nd ldnk is not used
call splsm2(x,y,w,n,match,nef,spar,dof,smo,s0,lev,ifcov,
		xin,yin,win,knot,
#       coef(nk),sout(nef+1),  levout(nef+1), xwy(nk),
#	hs0(nk),           hs1(nk),           hs2(nk),
#       hs3(nk),
#	sg0(nk),           sg1(nk),           sg2(nk),
#       sg3(nk),
#	abd(ld4,nk),   p1ip(ld4,nk),
#       p2ip(ldnk,nk)
	work(1), work(nk+1),   work(nk+nef+2),work(nk+2*nef+3),
	work(2*nk+2*nef+3),work(3*nk+2*nef+3),work(4*nk+2*nef+3),
	work(5*nk+2*nef+3),
	work(6*nk+2*nef+3),work(7*nk+2*nef+3),work(8*nk+2*nef+3),
	work(9*nk+2*nef+3),
	work(10*nk+2*nef+3),work((10+ld4)*nk+2*nef+3),
	work((10+2*ld4)*nk+2*nef+3),
	ld4,ldnk,nk)

return
end
subroutine splsm2(x,y,w,n,match,nef,spar,dof,smo,s0,lev,ifcov,
		xin,yin,win,knot,
		coef,sout,levout,xwy,
		hs0,hs1,hs2,hs3,
		sg0,sg1,sg2,sg3,
		abd,p1ip,p2ip,ld4,ldnk,nk)
implicit double precision(a-h,o-z)
integer n,match(n),nef
double precision x(n),y(n),w(n),spar,dof,smo(n),s0,lev(n)
integer nk,ldnk,ld4
logical ifcov
double precision xin(nef+1),yin(nef+1),win(nef+1),knot(nk+4)
double precision coef(nk),sout(nef+1),levout(nef+1),xwy(nk),
	 	hs0(nk),hs1(nk),hs2(nk),hs3(nk),
         	sg0(nk),sg1(nk),sg2(nk),sg3(nk),
	 	abd(ld4,nk),p1ip(ld4,nk),p2ip(ldnk,1)
#    local variables
integer ispar,icrit,isetup,ier
double precision lspar,uspar,tol,penalt,
		sumwin,dofoff,crit,xbar,dsum
double precision wmean
crit=0e0
if(spar==0e0){
	ispar=0
	dofoff=0e0
	if(dof==0e0) 
		{icrit=2 }
	else 
		{dofoff=dof+1e0;icrit=3}
	}
else {ispar=1;dofoff=0e0;icrit=3}
isetup=0
ier=1
penalt=1e0
lspar=1e-10
uspar=1.5
tol=1e-3

call sbart(penalt,dofoff,xin,yin,win,nef,knot,nk,
		  coef,sout,levout,
		  crit,icrit,spar,ispar,lspar,uspar,tol,
		  isetup,
		  xwy,
		  hs0,hs1,hs2,hs3,
         	  sg0,sg1,sg2,sg3,
		  abd,p1ip,p2ip,ld4,ldnk,ier)

#
# if ier ne 0, repeat the above call reducing uspar
#
nit=0
while((ier != 0) & (nit <7)) {
  crit=0e0;spar=0
  if(spar==0e0){
	ispar=0
	dofoff=0e0
	if(dof==0e0) 
		{icrit=2 }
	else 
		{dofoff=dof+1e0;icrit=3}
	}
   else {ispar=1;dofoff=0e0;icrit=3}
  isetup=0
   ier=1
   penalt=1e0
   lspar=1e-10
   tol=1e-3

   nit=nit+1
  uspar=.9+(uspar-.9)*.5
  call sbart(penalt,dofoff,xin,yin,win,nef,knot,nk,
		  coef,sout,levout,
		  crit,icrit,spar,ispar,lspar,uspar,tol,
		  isetup,
		  xwy,
		  hs0,hs1,hs2,hs3,
         	  sg0,sg1,sg2,sg3,
		  abd,p1ip,p2ip,ld4,ldnk,ier)
}



#return
#now clean up 
dof=0e0
sumwin=0e0
do i=1,nef {
	win(i)=win(i)*win(i) #we sqrted them in sbart
	}
sbar=wmean(nef,sout,win)
xbar=wmean(nef,xin,win)
do i=1,nef 
	sumwin=sumwin+win(i)
s0=wmean(n,y,w)
do i=1,nef {lev(i)=(xin(i)-xbar)*sout(i)	}
xsbar=wmean(nef,lev,win)
do i=1,nef {lev(i)=(xin(i)-xbar)**2	}
dsum=wmean(nef,lev,win)
do i=1,nef {
	if(win(i)>0e0) {
		lev(i)=levout(i)/win(i)-1e0/sumwin -lev(i)/(sumwin*dsum)
		}
	else {lev(i)=0e0}
	}
do i=1,nef {dof=dof+lev(i)*win(i)}
dof=dof+1e0
do i=1,nef
	sout(i)=sout(i)-sbar -(xin(i)-xbar)*xsbar/dsum
call unpack(n,nef,match,sout,smo)
return
end


double precision function wmean(n,y,w)
integer n
double precision y(n),w(n),wtot,wsum
wtot=0e0
wsum=0e0
do i=1,n{
	wsum=wsum+y(i)*w(i)
	wtot=wtot+w(i)
}
if(wtot > 0e0) {wmean=wsum/wtot} else {wmean=0e0}
return
end
