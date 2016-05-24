metareg <-
function(y, x, v, threshold=0.00001, maxiter=100, npermut=1000) {

	n=length(y)
	if(!(is.matrix(x))) stop("x must be a matrix")
	p=ncol(x)
	detmat=det(t(x) %*% x)
	
	if (detmat < .00001) stop("Multicolinearity in the moderators\nSOLUTION: Try omitting one moderator\n");
	
	ll=c(rep(0,(p+1)))
	lldif=c(rep(0,p))
	lldifb=lldif
	
	pll=c(rep(0,p))
	pllb=pll
	pp=pll
	pz=pll
	ptd=pll
		
	scalfact=c(rep(1,p))
	
	xx=x
	
	# i=1 yields estimates, i=2...p+1 yields tests
	for (i in 1:(p+1)) {
		
		if (i==1) {x=xx}
		if (i!=1) {x=xx[,-(i-1),drop=FALSE]}
		
		aa=estim(y,x,v,threshold,maxiter)
		
		if (aa[3]==1) {stop("Random effects variance estimation did not converge\nTry increasing maxiter\n")}
 		if (aa[3]==-1) {warning("Variance estimate has been set to zero\n")}
		
		sigmasq=aa[1]; iter=aa[2]
		w=diag(1/(v + sigmasq))
		if(p==1 & i==2) {
			b=matrix(0,1,1) # this covers the case for testing the intercept only model
			ll[i]=sum(log(v+sigmasq))+sum(((y)^2)/(v+sigmasq))
		} else {
			b=(solve(t(x) %*% w %*% x)) %*% t(x) %*% w %*% y
			ll[i]=sum(log(v+sigmasq))+sum(((y-(x %*% b))^2)/(v+sigmasq))
		}
		
		if(i==1) {
			bstore=b
			wstore=w
			se= sqrt(diag(solve(t(x) %*% w %*% x)))
			tval=bstore/se
		}
		
		if(i!=1) {
			
			if(p==1 & i==2) residual=y
			else residual= y- (x %*% b)        # required for permutation test
						
			pp[i-1]=permute(b,y,x,v,residual,xx,tval[i-1],npermut,(i-1),threshold,maxiter)   #output is p value & lb & ub
			lldif[(i-1)]=ll[i]-ll[1] # ll ratio test
			pll[i-1]=pchisq(lldif[(i-1)], 1, lower.tail = FALSE) # chis square tests are always two sided
			scalfact[i-1]=bartlett(wstore,xx,bstore[i-1],(i-1),n,p) # bartlett
			lldifb[(i-1)]=scalfact[i-1]*lldif[i-1] 
			pllb[i-1]=pchisq(lldifb[(i-1)], 1, lower.tail = FALSE)
		}
	}
	
	# results on residual variance estimate    	    
 	w_f= diag(1/v)    
 	b_f= (solve(t(xx) %*% wstore %*% xx)) %*% t(xx) %*% wstore %*% y
 	ll_f=sum(log(v))+sum(((y-(xx %*% b_f))^2)/v)
 	lldif_f=ll_f-ll[1]
 	llpval_f=ifelse(all.equal(lldif_f,0), pchisq(lldif_f, 1, lower.tail = FALSE), 0.5*pchisq(lldif_f, 1, lower.tail = FALSE))
 	varmat=cbind(sigmasq,lldif_f,1,llpval_f) 
	
	## z method: pvalue
	pz=pnorm(abs(tval),lower.tail = FALSE, log.p = FALSE)*2	            # two sided
	
	## t method: pvalue
	ptd=pt(abs(tval), (n-p), lower.tail = FALSE, log.p = FALSE)*2	    # two sided	
	
	outp=list(convergence=aa[3],iter=iter, 
		coefficients=bstore, se=se, tval=tval, 
		pLLR=pll, pBartlett=pllb, pZtest=pz, pttest=ptd, ppermtest=pp, 
		LLR=lldif, bartLLR=lldifb, bartscale=scalfact, dfttest=n-p,
		variance=varmat)
		
	return(outp)
}

