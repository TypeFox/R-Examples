subroutine bruto(x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,
		dfmax,cost,lambda,df,coef,type,xrange,
		gcvsel,gcvbak,dfit,maxit,nit,
		eta,resid,thresh,work,trace)
implicit double precision(a-h,o-z) 
double precision x(n,q),y(n,p),w(n),knot(nkmax+4,q),wp(p),
	 	dfmax(q),cost,lambda(q),df(q),coef(nkmax*p,q),xrange(2,q),
		gcvsel(q,maxit(1)),gcvbak(q,maxit(2)),dfit(q,maxit(1)),
		eta(n,p),resid(n,p),thresh,work(1)
integer n,q,p,nkmax,nk(q),match(n,q),nef(q),type(q),maxit(2),
		nit(2)
logical trace, select

#Compute residuals
do j=1,p{
	do i=1,n{resid(i,j)=y(i,j)-eta(i,j)}
}

#bruto backfitting with selection
select=.true.
call bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,
		dfmax,cost,lambda,df,coef,type,xrange,
		gcvsel,dfit,maxit(1),nit(1),eta,resid,thresh*10d0,work,trace)
#regular backfitting
select=.false.
call bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,
		dfmax,cost,lambda,df,coef,type,xrange,
		gcvbak,dfit,maxit(2),nit(2),eta,resid,thresh,work,trace)
do j=1,p{
	do i=1,n{eta(i,j)=y(i,j)-resid(i,j)}
}
return
end
subroutine bakssp(select,x,n,q,y,p,w,knot,nkmax,nk,wp,match,nef,
		dfmax,cost,lambda,df,coef,type,xrange,
		gcv,dfit,maxit,nit,s,resid,thresh,work,trace)
implicit double precision(a-h,o-z) 
double precision x(n,q),y(n,p),w(n),knot(nkmax+4,q),wp(p),
	 	dfmax(q),cost,lambda(q),df(q),coef(nkmax*p,q),
		xrange(2,q),gcv(q,maxit),dfit(q,maxit),
		s(n,p),resid(n,p),thresh,work(1)
integer n,q,p,nkmax,nk(q),match(n,q),nef(q),type(q),maxit,nit,ier,ntype
double precision   dfoff, gcv0, ndf,gcv1,gcvrat,ndfoff,wmean,sbar,rss,tol
logical center, select, trace
center=.true.
tol=1d-3  # this is the convergence criterion for gcv
#remove the mean from the residuals, and compute rss
rss=0d0
do j=1,p{
	sbar=wmean(n,resid(1,j),w)
	do i=1,n{
		resid(i,j)=resid(i,j) - sbar
		work(i)=resid(i,j)**2
		}
rss=rss+wp(j)*wmean(n,work,w)
}
dfoff=0
do k=1,q {dfoff=dfoff+df(k)}
gcv1=rss/((1-(1+dfoff*cost)/n)**2)
gcvrat=1d0
nit=0
while(nit<maxit&gcvrat >thresh ){
	gcv0=gcv1
	nit=nit+1
	do k=1,q{
		gcv(k,nit)=gcv1
		if(!select&type(k)==1)next
# form partial residuals if necessary
		if(type(k)>1){ #get the fitted values
			call psspl2(x(1,k),n,p,knot(1,k),nk(k),xrange(1,k),
				coef(1,k),coef(1,k),s,0,type(k))
			do j=1,p{
				sbar=wmean(n,s(1,j),w)
				do i=1,n{resid(i,j)=resid(i,j) + s(i,j)-sbar}
			}
		}

		ndfoff=dfoff-df(k)
		if(select) {ntype=0}
		else {ntype=type(k)}
		call sspl2(x(1,k),resid,w,n,p,knot(1,k),nk(k),wp,match(1,k),
		nef(k),ndfoff,dfmax(k),cost,lambda(k),ndf,
			gcv1,coef(1,k),s,ntype,center,
			xrange(1,k),work,tol,ier)
		gcv(k,nit)=gcv1
		if(select){
			dfit(k,nit)=ndf
			df(k)=ndf
			dfoff=ndfoff+ndf
			type(k)=ntype
		}
		if(type(k)>1){
			do j=1,p{
				do i=1,n{
					resid(i,j)=resid(i,j) - s(i,j)
					}
			}
		}
	}
gcvrat=dabs(gcv1-gcv0)/gcv0
if(trace){
	call intpr("outer iteration",15,nit,1)
	call dblepr("gcv  ",5,gcv1,1)
	call dblepr("ratio",5,gcvrat,1)
}
}
return
end
subroutine pbruto(x,n,q,ybar,p,knot,nkmax,nk,coef,type,xrange,
		eta,work)
implicit double precision(a-h,o-z) 
double precision x(n,q),ybar(p),knot(nkmax+4,q),coef(nkmax*p,q),
		xrange(2,q),eta(n,p),work(n,p)
integer n,q,p,nkmax,nk(q),type(q)
#initialization
do j=1,p{
	do i =1,n {eta(i,j)=ybar(j)}
}
do k=1,q{
	if(type(k)==1)next
	call psspl2(x(1,k),n,p,knot(1,k),nk(k),xrange(1,k),
			coef(1,k),coef(1,k),work,0,type(k))
	do j=1,p{
		do i=1,n{eta(i,j)=eta(i,j) + work(i,j)}
	}
}
return
end
