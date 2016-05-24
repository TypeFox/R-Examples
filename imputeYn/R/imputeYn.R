aft.kmweight<-function(Y, delta)
	{
	srt<-order(Y)
	sy<-as.double(Y[srt])
	sdelta<-as.integer(delta[srt])
	n <- length(sdelta)
	if(n != length(sdelta) || n != length(Y))
	stop("dimensions of Y and delta don't match!")
	kmweights <- numeric(n)
	kmweights[1] <- 1/n
	for(i in 2 : n)
	 {
	 kmweights[i] <- kmweights[i-1] * (n-i+2)/(n-i+1) * (((n-i+1)/(n-i+2))^sdelta[i-1])
	 }
	kmwts<-kmweights*sdelta
	if (sdelta[n]==0)
	kmwts[n]<-1-sum(kmwts)
	return(list(kmwts=kmwts))
}



aft.qp<-function(X, Y, delta)
{
	n <- nrow(X) # number of samples
	p <- ncol(X) # number of predictors
	if(n != length(delta) || n != length(Y))
	stop("dimensions of X, Y and delta don't match!")
	eta<-0.01*sqrt(2*log(p))
	ndelta <- n - sum(delta)
	kw <- aft.kmweight(Y,delta)$kmwts
	XW <- apply(as.matrix(X[delta == 1,] * kw[delta == 1]), 2,
	sum) / sum(kw[delta == 1])
	YW <- sum(Y[delta == 1] * kw[delta == 1]) /
	sum(kw[delta == 1])
	for(i in 1:n)
	X[i,] <- X[i,] - XW
	X <- as.matrix(sqrt(kw) * X)
	Y <- sqrt(kw) * (Y - YW)
	Dmat <- t(X[delta==1, ]) %*% X[delta==1, ] + diag(eta, p)
	dvec <- drop(Y)[delta==1] %*% X[delta==1, ]
	Amat <- matrix(0, p, 1) # a psudo constraint just to make solve.QP work
	bvec <- -1
	qpobj <- solve.QP(Dmat,dvec,Amat,bvec)
	beta0 <- YW - XW %*% qpobj$solution
	return(list(beta0=beta0 , beta=qpobj$solution))
}



lss.mod<-function(formula, data, subset, trace=FALSE, mcsize=500, maxiter=10, tolerance=0.001, cov=FALSE, na.action=na.exclude, residue=FALSE)
{
    lss.eres<- function(e,z,delta,eps,sh=F)
	 {
	nobs=length(e)
	ord <- order(e)
	ei <- e[ord]
	zi <- z[ord]
	deltai <- delta[ord]
	tie <- c(diff(ei)>eps,1)
	tie1 <- c(1,diff(ei)>eps)
	dummy <- 1:nobs
	repeats <- diff(c(dummy[c(T,diff(ei)>eps)],nobs+1))
	Ni <- rev(cumsum(rev(zi)))
	di=cumsum(zi*deltai)
	di=di[tie>eps]
	di=c(di[1],diff(di))
	ieb <- 1 - di / Ni[tie1>eps]
	Shat <- cumprod(ieb)
	if(sh)
	{
		return(Shat)	
	}
	Shat <- rep(Shat,repeats)
	edif <- c(diff(ei),0)
	ehat <- rev(cumsum(rev(edif * Shat)))
	ehat[Shat<eps] <- 0
	Shat[Shat<eps] <- 1
	ehat <- ehat/Shat + ei
	eres <- ehat
	eres[dummy[ord]] <- ehat
	eres
       }

     lss.betag<-function(x,y,delta,z)
       {
	row=ncol(x)
	col=ncol(z)		
	betagm<-matrix(0,ncol=col,nrow=row)
	ynew<-1000*(length(y))^2
	dimnum<-dim(x)
	n1<-dimnum[1]
	n2<-dimnum[2]	
	yy0<-rep(y,rep(n1,n1))
	delta1<-rep(delta,rep(n1,n1))
	yy1<-rep(y,n1)
	yy2<-delta1*(yy0-yy1)
	xx0<-matrix(rep(as.vector(x),rep(n1,n1*n2)),nrow=n1*n1)
	xx1<-t(matrix(rep(as.vector(t(x)),n1),nrow=n2))
	xx2<-xx0-xx1	
	for(i in 1:col)
	{
		zz=rep(z[,i],rep(n1,n1))*rep(z[,i],n1)
		xxdif<-xx2*zz*delta1
		xnew<-apply(xxdif,2,sum)
		xnew<-rbind(xxdif)
		yynew<-c(yy2*zz)
		fit <- aft.qp(xnew, yynew, delta1)$beta
		betagm[,i] <- fit
	}
	betagm
	 }

	eps <- .Machine$double.eps^(2/3)
	call <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, sys.parent())
	Terms <- attr(mf, "terms")
	xvars <- as.character(attr(Terms, "variables"))
	yvar <- attr(Terms, "response")
	if((yvar <- attr(Terms, "response")) > 0)
		xvars <- xvars[ - yvar]
	else xlevels <- NULL
	y <- model.extract(mf,"response")
	x <- model.matrix(Terms, mf)
	if(all(x[, 1] == 1))
		x <- x[, -1]			
	if(ncol(as.matrix(y)) != 2)
		stop("Response must be a right-censored survival object!")
	nobs <- nrow(y)
	nvar <- ncol(x)

	fit <- list(converged = FALSE, mcsize=mcsize)
	class(fit) <- c("lss.mod")
	fit$cnames <- dimnames(x)[[2]]
	fit$niter <- 0
	fit$printkm <- residue
	z <- matrix(rexp(nobs*mcsize), ncol=mcsize)
	zdummy <- matrix(rep(1,nobs), ncol=1)
	beta <- lss.betag(x, y[,1], y[,2], zdummy)
	betastar <- lss.betag(x, y[,1], y[,2], z)	
	bbar <- apply(betastar, 1, mean)
	tmp <- betastar - bbar
	if(trace)
		cat("\nbetag: ", format(beta), "\n\n")

	niter=0
	xbar <- apply(x, 2, mean)
	xm <- x - rep(xbar, rep(nobs, nvar))
	xinv <- solve(t(xm) %*% xm)
	xinvstar <- array(1,dim=c(nvar,nvar,mcsize))
	for(i in 1:mcsize)
	{
		xinvstar[,,i] <- solve(t(xm) %*% (xm*z[,i]))
	}
	
	while(niter < maxiter)
	{
		niter <- niter + 1
		betaprev <- beta

		e <- y[,1] - x %*% beta
		eres <- lss.eres(e, zdummy, y[,2], eps)
		yhat <- y[,2] * y[,1] + (1 - y[,2]) * (eres + x %*% beta)
		ybar <- mean(yhat)
		beta <- xinv %*% (t(xm) %*% (yhat - ybar))

		if(trace)
		{
			cat("Iteration: ", niter)
			cat("\n  Beta: ", format(beta), "\n")
		}

		for(i in 1:mcsize)
		{
			e <- y[,1] - x %*% betastar[,i]
			eres <- lss.eres(e, z[,i], y[,2], eps)
			yhat <- y[,2] * y[,1] + (1 - y[,2]) * (eres + x %*% betastar[,i])
			ybar <- mean(yhat)
			betastar[,i] <- xinvstar[,,i] %*% (t(xm*z[,i]) %*% (yhat - ybar))
		}
		
		bb <- abs(beta)
		bb[bb<0.01] <- 0.01
		mm <- max(abs(beta - betaprev) / bb)
		if(mm < tolerance)
		{
			fit$residue <- y[,1] - x %*% beta
			fit$km.residue <- lss.eres(fit$residue, zdummy, y[,2], eps, sh=TRUE)
			fit$converged <- TRUE
			break
		}
	}

	fit$niter <- niter
	fit$lse <- beta
	bbar <- apply(betastar, 1, mean)
	tmp <- betastar - bbar
	fit$tol <- tolerance
	fit$residue
	fit$covmatrix <- tmp %*% t(tmp)/(mcsize - 1)
	fit$sd <- sqrt(diag(fit$covmatrix))
	fit$zvalue <- beta/fit$sd
	fit$pvalue <- (1 - pnorm(abs(fit$zvalue))) * 2
	dimnames(fit$covmatrix) <- list(fit$cnames,fit$cnames)
	fit
}


imputeYn<-function(X, Y, delta, method = "condMean", beta=NULL)
    {
    nn <- length(Y) # number of samples
    srt<-order(Y) 
    sdelta<-as.integer(delta[srt])
    if(sdelta[nn]!=0)
    stop("The largest observation is not censored!")

	if(method=="PDQ"){
	beta<-NULL
	n<-nrow(Y)}
	else
      {
	call <- match.call()
	n <- nrow(X) # number of samples
	p <- ncol(X) # number of predictors
    	if (is.null(beta))
	beta <- c(aft.qp(X,Y,delta)$beta0,aft.qp(X,Y,delta)$beta)
	delta[n]<-1
      }

    condMean<-function(X, Y, delta, beta) 
        {
	N<-length(delta)
	u<-cbind(1,X) %*% beta
	res<-Y-u
	niceorder<-order(as.vector(res), -delta)
	kk<-which( niceorder == N ) #### this should work with tie
	resorder<-res[niceorder]
	dorder<-delta[niceorder]
	dorder[N]<-1
	uorder<-u[niceorder]
	ystar<-Y[niceorder]
	xorder<-as.matrix(X[niceorder, ])
	temp<-WKM(x = resorder, d = dorder, zc = 1:N)
	jifen<-rev(cumsum(rev(resorder * temp$jump)))
	Sresorder<-temp$surv
	if(Sresorder[kk]>0) {
	ystar[kk]<-uorder[kk] + jifen[kk]/(Sresorder[kk])
	}
	return( ystar[kk] )
	  }

    condMedian<-function(X, Y, delta, beta) 
        {
	N<-length(delta)
	u<-cbind(1,X) %*% beta
	res<-Y - u
	delta[N]<-1 ### do we need this?
	niceorder<-order(as.vector(res), -delta)
	kk<-which(niceorder == N )
	resorder<-res[niceorder]
	dorder<-delta[niceorder]
	dorder[N]<-1
	uorder<-u[niceorder]
	ystar<-Y[niceorder]
	xorder<-as.matrix(X[niceorder, ])
	temp<-WKM(x = resorder, d = dorder, zc = 1:N)
	Sresorder<-temp$surv
	if(Sresorder[kk] >0 ) {
	Ptheta<-Sresorder[kk]/2
	pivec<-temp$jump
	pivec[1:(kk-1)]<-0
	posi<-sum( cumsum(pivec) < Ptheta )
	theta<-temp$times[ posi +1]
	ystar[kk]<-uorder[kk] + theta
	}
	return(ystar[kk])
        }

    RcondMean<-function(X, Y, delta) 
         {
	N<-length(delta)
	res<-lss.mod(cbind(Y,delta) ~ X, mcsize=100,trace=FALSE,maxiter=20,tolerance=0.1)$residue
	u<-Y-res
	niceorder<-order(as.vector(res), -delta)
	kk<-which( niceorder == N ) #### this should work with tie
	resorder<-res[niceorder]
	dorder<-delta[niceorder]
	dorder[N]<-1
	uorder<-u[niceorder]
	ystar<-Y[niceorder]
	xorder<-as.matrix(X[niceorder, ])
	temp<-WKM(x = resorder, d = dorder, zc = 1:N)
	jifen<-rev(cumsum(rev(resorder * temp$jump)))
	Sresorder<-temp$surv
	if(Sresorder[kk]>0) {
	ystar[kk]<-uorder[kk] + jifen[kk]/(Sresorder[kk])
	}
	return( ystar[kk] )
        }

    RcondMedian<-function(X, Y, delta) 
        {
	N<-length(delta)
	res<-lss.mod(cbind(Y,delta) ~ X, mcsize=100,trace=FALSE,maxiter=20,tolerance=0.1)$residue
	u<-Y-res
	delta[N]<-1 ### do we need this?
	niceorder<-order(as.vector(res), -delta)
	kk<-which(niceorder == N )
	resorder<-res[niceorder]
	dorder<-delta[niceorder]
	dorder[N]<-1
	uorder<-u[niceorder]
	ystar<-Y[niceorder]
	xorder<-as.matrix(X[niceorder, ])
	temp<-WKM(x = resorder, d = dorder, zc = 1:N)
	Sresorder<-temp$surv
	if(Sresorder[kk] >0 ) {
	Ptheta<-Sresorder[kk]/2
	pivec<-temp$jump
	pivec[1:(kk-1)]<-0
	posi<-sum( cumsum(pivec) < Ptheta )
	theta<-temp$times[ posi +1]
	ystar[kk]<-uorder[kk] + theta
	}
	return(ystar[kk])
	   }

    pred.diff <- function(Y, delta) {
  	sorted<-order(Y) 
  	Y<-as.double(Y[sorted])
  	delta<-as.integer(delta[sorted])
  	status <- delta
  	yy <- Y
  	N <- length(yy)
  	timeorig <- yy
  	order.orig <- 1:N
  	dummystrat <- factor(rep(1, N))
  	sse <- 0
  	n <- 0	
  	nonconv <- FALSE	
  	oldsse <- sse
  	one <- rep(1, N)
  	state <- status
  	state[yy == max(yy)] <- 1
  	m2 <- order(yy)
  	yy.order <- yy[m2]
  	surv.obj <- Surv(yy,state)~1
	KM.hat <- KM.ehat <- survfit(surv.obj,type="kaplan-meier",conf.type = "none", se.fit = FALSE)
  	n.risk <- KM.ehat$n.risk
  	surv <- KM.ehat$surv
  	repeats <- c(diff( - n.risk), n.risk[length(n.risk)])
  	surv <- rep(surv, repeats)
  	w <-  - diff(c(1, surv))
  	m2 <- order(Y)
  	log.y.order <- yy[m2]
  	bla <- NA
  	for(iiii in 1:length(log.y.order)) {
    	bla[iiii] <- sum((log.y.order*w)[-(1:iiii)])
  	}
  	bla2 <- bla/surv
  	bl <- bla2
  	bl[(1:N)[m2]] <- bla2
  	yhat <- bl
  	yy.new <- yy
  	yy.new[state == 0] <- yhat[state == 0]
  	yy.imp <- yhat[state==0]
  	le<-length(yy.imp)
  	med.dev <-sum(abs(yy.imp-median(yy.imp)))/le
  	yy.cens<-yy[state==0]
  	diff<-yy.imp-yy.cens
  	f <- list()
  	f$Ynew<-yy.new
  	f$Ycens<-yy.cens
  	f$imp<-yy.imp
  	f$diff<-diff
  	D<-f$diff[1:(length(f$diff)-0)]
  	Censored.Y<-f$Ycens[1:(length(f$Ycens)-0)]
  	Ymax<-f$Ynew[length(f$Ynew)]
  	w<-1/(Ymax-Censored.Y)
  	reg<-lm(D~Censored.Y, weights=w)
  	diff.pred<-reg$coefficients[1]+reg$coefficients[2]*Ymax
  	Ymax.hat<-Ymax+diff.pred
  	f$diff.pred<-diff.pred
  	f$yyN<-Ymax.hat
  	f
}

	if(method == "condMean"){
	Y[n]<-condMean(X,Y,delta,beta)
	newbeta <- aft.qp(X,Y,delta)$beta}

	else if(method == "condMedian"){
	Y[n]<-condMedian(X,Y,delta,beta)
	newbeta <- aft.qp(X,Y,delta)$beta}

	else if(method == "RcondMean"){
	Y[n]<-RcondMean(X,Y,delta)
	newbeta <- aft.qp(X,Y,delta)$beta}

	else if(method == "RcondMedian"){
	Y[n]<-RcondMedian(X,Y,delta)
	newbeta <- aft.qp(X,Y,delta)$beta}

	else if(method == "PDQ"){
	Y[n]<-pred.diff(Y,delta)$yyN
	newbeta <- NULL}
	else stop("invalid treating method for Y(n)+")
	weighting=list(Yn=Y[n],newdata=cbind(Y=Y,delta=delta),newcoefficients=newbeta,coefficients=beta[-1])
	class(weighting)<-("imputeYn")
	return(weighting)
}


print.imputeYn<-function(x, digits = max(3, getOption("digits") - 3), ...)
    	{
	cat("Y(n):\n")
	print(x$Yn)

	if (length(coef(x))) 
		{
		cat("Coefficients:\n")
		print.default(format(coef(x), digits = digits), print.gap = 2,
		quote = FALSE)
		cat("New Coefficients:\n")
		print.default(format(x$newcoefficients, digits = digits), print.gap = 2,
		quote = FALSE)
		}
	else cat("No coefficients\n")
	cat("\n")
	}


imputeYn.extra<-function(Y, delta, hc.Yn, method = "km.TPQ", trans.sprob=NULL, stime2=NULL, sprob2=NULL, trace=F)
    {
    nn <- length(Y) # number of samples
    srt<-order(Y) 
    sdelta<-as.integer(delta[srt])
    if(sdelta[nn]!=0)
    stop("The largest observation is not censored!")

    it.pred.diff<-function(Y, delta, hc.Yn, trace)
	{
	pred.diff <- function(Y, delta) {
  	sorted<-order(Y) 
  	Y<-as.double(Y[sorted])
  	delta<-as.integer(delta[sorted])
  	status <- delta
  	yy <- Y
  	N <- length(yy)
  	timeorig <- yy
  	order.orig <- 1:N
  	dummystrat <- factor(rep(1, N))
  	sse <- 0
  	n <- 0	
  	nonconv <- FALSE	
  	oldsse <- sse
  	one <- rep(1, N)
  	state <- status
  	state[yy == max(yy)] <- 1
  	m2 <- order(yy)
  	yy.order <- yy[m2]
  	surv.obj <- Surv(yy,state)~1
  	KM.hat <- KM.ehat <- survfit(surv.obj,type="kaplan-	meier",conf.type = "none", se.fit = FALSE)
  	n.risk <- KM.ehat$n.risk
  	surv <- KM.ehat$surv
  	repeats <- c(diff( - n.risk), n.risk[length(n.risk)])
  	surv <- rep(surv, repeats)
  	w <-  - diff(c(1, surv))
  	m2 <- order(Y)
  	log.y.order <- yy[m2]
  	bla <- NA
  	for(iiii in 1:length(log.y.order)) {
    	bla[iiii] <- sum((log.y.order*w)[-(1:iiii)])
  	}
  	bla2 <- bla/surv
  	bl <- bla2
  	bl[(1:N)[m2]] <- bla2
  	yhat <- bl
  	yy.new <- yy
  	yy.new[state == 0] <- yhat[state == 0]
  	yy.imp <- yhat[state==0]
  	le<-length(yy.imp)
  	med.dev <-sum(abs(yy.imp-median(yy.imp)))/le
  	yy.cens<-yy[state==0]
  	diff<-yy.imp-yy.cens
  	f <- list()
  	f$Ynew<-yy.new
  	f$Ycens<-yy.cens
  	f$imp<-yy.imp
  	f$diff<-diff
  	D<-f$diff[1:(length(f$diff)-0)]
  	Censored.Y<-f$Ycens[1:(length(f$Ycens)-0)]
  	Ymax<-f$Ynew[length(f$Ynew)]
  	w<-1/(Ymax-Censored.Y)
  	reg<-lm(D~Censored.Y, weights=w)
  	diff.pred<-reg$coefficients[1]+reg$coefficients[2]*Ymax
  	Ymax.hat<-Ymax+diff.pred
  	f$diff.pred<-diff.pred
  	f$yyN<-Ymax.hat
  	f
    }

     sfit<- survfit(Surv(Y, delta == 1) ~ 1)
	im<-pred.diff(Y, delta)
	D<-im$diff[1:(length(im$diff))]
	Censored.Y<-im$Ycens[1:(length(im$Ycens))]
	max<-hc.Yn[1]
	w<-1/(max-Censored.Y)
	reg<-lm(D~Censored.Y, weights=w)
	im.hc.Yn<-NULL
      for (i in 1:(length(hc.Yn)-1))
        {
	  im.hc.Yn[1]<-hc.Yn[1]+reg$coefficients[1]+reg$coefficients[2]*hc.Yn[1]
	  im.hc.Yn[i+1]<-im.hc.Yn[i]+reg$coefficients[1]+reg$coefficients[2]*im.hc.Yn[i]
	  }	
	sorted<-order(Y) 
	sY<-Y[sorted]
	sdelta<-delta[sorted]
	cobs<-(length(sY)-length(hc.Yn)+1):length(sY)
	sY[cobs]<-im.hc.Yn
	sdelta[cobs]<-1
      sfit2<- survfit(Surv(sY, sdelta == 1) ~ 1)
	if(trace){
	par(mfrow=c(1,2))
	pl<-plot(sfit,col=3,xlab="Survival time", ylab="Survival probability",main="For original data",cex=0.75)
	p2<-plot(sfit2,col=3,xlab="Survival time", ylab="Survival probability",main="For data with imputation",cex=0.75)}      
	return(list(Y=Y[sorted], delta=delta[sorted], newY=sY, newdelta=sdelta, hc.Yn=hc.Yn, im.hc.Yn=im.hc.Yn, trace=trace))
} 
   

km.trend.pred<-function(Y, delta, hc.Yn, trans.sprob=NULL, stime2=NULL, sprob2=NULL, trace)
	{
	k<-nchar(length(hc.Yn))
      sfit<- survfit(Surv(Y, delta == 1) ~ 1)
      stime<-stime.extra<-sfit$time
      sprob<-sprob.extra<-sfit$surv
      if (is.null(stime2))
 	{stime<-stime
      sprob<-sprob}     
	else
      {stime<-stime2
      sprob<-sprob2}
      reg<-lm(stime~sprob)
      sprob.more<-0:(sprob.extra[length(sprob.extra)]*10^k)
	sprob.sample<-sample(sprob.more,length(hc.Yn))/10^k
	if(is.null(trans.sprob))
	sprob.sample<-sprob.sample
	else if (trans.sprob=="exp")
	sprob.sample<-exp(sprob.sample)
	else
	sprob.sample<-sprob.sample^trans.sprob
	imp.hc.Yn<-NULL
        for (i in 1:length(sprob.sample))
         {
        imp.hc.Yn[i]<-reg$coefficients[1]+reg$coefficients[2]*sprob.sample[i]
         }
	sorted<-order(Y) 
	sY<-Y[sorted]
	sdelta<-delta[sorted]
	cobs<-(length(sY)-length(hc.Yn)+1):length(sY)
	sY[cobs]<-sort(imp.hc.Yn)
	sdelta[cobs]<-1
      sfit2<- survfit(Surv(sY, sdelta == 1) ~ 1)
	if(trace){
	par(mfrow=c(1,2))
	pl<-plot(sfit,col=3,xlab="Survival time", ylab="Survival probability",main="For original data",cex=0.75)
	p2<-plot(sfit2,col=3,xlab="Survival time", ylab="Survival probability",main="For data with imputation",cex=0.75)}      
	return(list(Y=Y[sorted], delta=delta[sorted], newY=sY, newdelta=sdelta, hc.Yn=hc.Yn, imp.hc.Yn=sort(imp.hc.Yn), stime=stime.extra,sprob=sprob.extra,
	stime2=stime2, sprob2=sprob2, trans=trans.sprob, trace=trace))
}


	if(method == "it.PDQ"){
     it.PDQ<-it.pred.diff(Y, delta, hc.Yn, trace)
	return(it.PDQ)}

	else if(method == "km.TPQ"){
	km.TPQ<-km.trend.pred(Y, delta, hc.Yn, trans.sprob=NULL, 	stime2=NULL, sprob2=NULL,  trace)
	return(km.TPQ)}

	else stop("invalid additional imputing method for imputing more than 	one Y(n)+ observation")
}
	
data<-function(n, p, r, b1, sig, Cper)
	{
      r.imp<-function(a,r)
		{
	jeroval<-a[]==0
	imputed<-a
	imputed[jeroval]<-r
	return(imputed)
		}
	diagmat<-diag(p)
	cormat<-r.imp(diagmat,r)
	rr<-chol(cormat)
	z<-matrix(runif(n*p),ncol=p)
	x<-z%*%rr
	e<-rnorm(n)
	y<-x%*%b1+e*sig
      repeat{
		repeat{
		c<-runif(n,range(y)[1],range(y)[2]-range(y)[2]*Cper)
		z<-rep(NA,n) 
		d<-rep(NA,n) 
		for(i in 1:n)
			  {
			if (y[i]<c[i]){
			z[i]<-y[i]
			d[i]<-1}
			else {z[i]<-c[i]
			d[i]<-0}
			  }
		sorted<-order(z) 
		sz<-as.double(z[sorted])
		sy<-as.double(y[sorted]) 
		sstat<-as.integer(d[sorted])
		sx<-x[sorted,]  
		if (sstat[n]==0) {break}
			}
	sz<-as.matrix(sz);sstat<-as.matrix(sstat);sx<-as.matrix(sx);sy<-as.matrix(sy)
	Pper<-sum(sstat==0)*100/length(sstat)
	if(Pper==Cper*100+50){break}
		}
	list(y=sz,x=sx,delta=sstat,Pper=Pper)
	}


