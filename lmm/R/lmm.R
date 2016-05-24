###############################################################################
fastmcmc <- function(y,subj,pred,xcol,zcol,prior,seed,vmax,occ,
	start.mode,maxits=100,eps=.0001,iter=1000,start.mcmc,df=4){
        tmp <- table(subj)
	m <- length(tmp)
	nmax <- max(tmp)
	ntot <- length(y)
	pcol <- ncol(pred)
	q <- length(zcol)
	p <- length(xcol)
	g <- as.integer(round(q*(q+1)/2))
	#
	{if(missing(vmax)){
		vmax <- diag(rep(1,nmax))
                occ <- integer(ntot)
		iflag <- as.integer(1)}
	else iflag <- as.integer(0)}
	storage.mode(vmax) <- "double"
	#	
	{if(!missing(start.mode)){
		beta <- start.mode$beta
		sigma2 <- start.mode$sigma2
		xi <- start.mode$psi/start.mode$sigma2
		storage.mode(beta) <- "double"
		storage.mode(xi) <- "double"
		storage.mode(sigma2) <- "double"
		sflag <- as.integer(1)}
	else{
		beta <- numeric(p)
		xi <- matrix(0,q,q)
		sigma2 <- 0
		sflag <- as.integer(0)}}
	abc <- c(prior$a,prior$b,prior$c)
        dinv <- prior$Dinv
	storage.mode(abc) <- "double"
	storage.mode(dinv) <- "double"
	cat("Performing FAST-MODE...")
	now <- proc.time()
        err <- 0
	tmp <- .Fortran("fastmode",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,vmax,
		w=array(0,c(nmax,nmax,m)),
		vinv=array(0,c(nmax,nmax,m)),
		pcol,as.double(pred),q,as.integer(zcol),
		ztvinv=array(0,c(q,nmax,m)),
		ztvinvz=array(0,c(q,q,m)),
		isflag=c(iflag,sflag),err=as.integer(err),
		msg=integer(1),u=array(0,c(q,q,m)),
		iter=integer(1),sigma2=sigma2,p,
		as.integer(xcol),beta=beta,as.double(y),
		delta=rep(0,ntot),xtw=matrix(0,p,nmax),
		xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),wkqq1=matrix(0,q,q),
		wkqq2=matrix(0,q,q),
		xi=xi,wkqnm=array(0,c(q,nmax,m)),
		b=matrix(0,q,m),cvgd=integer(1),
		oxi=matrix(0,q,q),maxits=as.integer(maxits),
		llvec=numeric(as.integer(maxits)),
		eps=as.double(eps),xiecme=matrix(0,q,q),
		g=g,reject=integer(maxits),
		ztvinvx=array(0,c(q,p,m)),
		a=array(0,c(q,q,m)),wkqp=matrix(0,q,p),
		wkg=rep(0,g+1),wkgg=matrix(0,g+1,g+1),
		abc=abc,dinv=dinv,PACKAGE="lmm")
	cat("\n")
	iter.mode <- tmp$iter
        msg <- tmp$msg
	{if(msg==1)
          warning("Supplied V <- i matrix is not positive definite")
	else if(msg==2)
          warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
	else if(msg==3)
          warning("Inadequate information to obtain starting value of psi")
	else if(msg==4)
          warning("Value of psi became non-pos.def. during iterations")
	else if(msg==5)
          warning("t(X)%*%W%*%X became non-pos.def. during iterations")
	else if(msg==10)
	  warning("non-pos.def. xi at input to scoring step")
	else if(msg==11)
	  warning("log-density not concave at one or more scoring steps")}
        llvec <- tmp$llvec[1:iter.mode]
	reject <- tmp$reject[1:iter.mode]
	converged <- tmp$cvgd==as.integer(1)
	if(!converged) warning(paste("FAST-MODE did not converge by",
	   format(iter.mode),"iterations"))
	#
	sigma2 <- tmp$sigma2
	beta <- tmp$beta
	xi <- tmp$xi
	psi <- xi*sigma2
	cov.beta <- tmp$sigma2*tmp$xtwxinv
	#
	mode.list <- list(beta=beta,sigma2=sigma2,psi=psi,
		converged=converged,iter=iter.mode,reject=reject,
		logpost=llvec,cov.beta=cov.beta)
	#
	{if(!missing(start.mcmc)){
		beta <- start.mcmc$beta
		sigma2 <- start.mcmc$sigma2
		xi <- start.mcmc$psi/start.mcmc$sigma2
		storage.mode(beta) <- "double"
		storage.mode(xi) <- "double"
		storage.mode(sigma2) <- "double"}}
	#
	wkgg <- tmp$wkgg
	sig2hat <- tmp$sigma2
	xihat <- tmp$xi
	#
	junk <- .Fortran("rngs",as.integer(seed),PACKAGE="lmm")
	cat("Performing FAST-MCMC...")
	tmp <- .Fortran("fastmcmc",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,vmax,
		w=array(0,c(nmax,nmax,m)),
		vinv=array(0,c(nmax,nmax,m)),
		pcol,as.double(pred),q,as.integer(zcol),
		ztvinv=array(0,c(q,nmax,m)),
		ztvinvz=array(0,c(q,q,m)),
		iflag=iflag,
		msg=integer(1),u=array(0,c(q,q,m)),
		sigma2=sigma2,p,
		as.integer(xcol),beta=beta,as.double(y),
		delta=rep(0,ntot),xtw=matrix(0,p,nmax),
		xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),wkqq1=matrix(0,q,q),
		wkqq2=matrix(0,q,q),
		xi=xi,wkqnm=array(0,c(q,nmax,m)),
		b=matrix(0,q,m),maxits=as.integer(iter),
		abc=abc,dinv=dinv,sqrtu=array(0,c(q,q,m)),
		sigma2s=rep(0,iter),psis=array(0,c(q,q,iter)),
		g=g,wkgg=wkgg,wkgg2=matrix(0,g+1,g+1),wkg=rep(0,g+1),
		sig2hat=sig2hat,xihat=xihat,xigibbs=matrix(0,q,q),
		reject=integer(iter),ratios=rep(0,iter),df=as.double(df),PACKAGE="lmm")
	cat("\n")
	clock <- proc.time()-now
	#
        psi <- tmp$xi*tmp$sigma2
        beta <- tmp$beta
	if(!is.null(dimnames(pred)[[2]])){
           colnames <- dimnames(pred)[[2]]
           names(beta) <- colnames[xcol]
           dimnames(psi) <- list(colnames[zcol],colnames[zcol])}
        #
	list(beta=beta,sigma2=tmp$sigma2,psi=psi,
		sigma2.series=tmp$sigma2s,psi.series=tmp$psis,
		ratios=tmp$ratios,accept=(tmp$reject==0),
		mode.list=mode.list)}
###############################################################################
fastmode <- function(y,subj,pred,xcol,zcol,prior,vmax,occ,start,maxits=100,
	eps=.0001){
        tmp <- table(subj)
	m <- length(tmp)
	nmax <- max(tmp)
	ntot <- length(y)
	pcol <- ncol(pred)
	q <- length(zcol)
	p <- length(xcol)
	g <- as.integer(round(q*(q+1)/2))
	maxits <- as.integer(round(maxits))
	#
	{if(missing(vmax)){
		vmax <- diag(rep(1,nmax))
                occ <- integer(ntot)
		iflag <- as.integer(1)}
	else iflag <- as.integer(0)}
	storage.mode(vmax) <- "double"
	#	
	{if(!missing(start)){
		beta <- start$beta
		sigma2 <- start$sigma2
		xi <- start$psi/start$sigma2
		storage.mode(beta) <- "double"
		storage.mode(xi) <- "double"
		storage.mode(sigma2) <- "double"
		sflag <- as.integer(1)}
	else{
		beta <- numeric(p)
		xi <- matrix(0,q,q)
		sigma2 <- 0
		sflag <- as.integer(0)}}
	abc <- c(prior$a,prior$b,prior$c)
        dinv <- prior$Dinv
	storage.mode(abc) <- "double"
	storage.mode(dinv) <- "double"
	cat("Performing FAST-MODE...")
	now <- proc.time()
	tmp <- .Fortran("fastmode",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,vmax,
		w=array(0,c(nmax,nmax,m)),
		vinv=array(0,c(nmax,nmax,m)),
		#
		pcol,as.double(pred),q,as.integer(zcol),
		ztvinv=array(0,c(q,nmax,m)),
		ztvinvz=array(0,c(q,q,m)),
		isflag=c(iflag,sflag),
		msg=integer(1),u=array(0,c(q,q,m)),
		iter=integer(1),
		#
		sigma2=sigma2,p,
		as.integer(xcol),beta=beta,as.double(y),
		delta=rep(0,ntot),xtw=matrix(0,p,nmax),
		xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),
		#
		wkqq1=matrix(0,q,q),wkqq2=matrix(0,q,q),
		xi=xi,wkqnm=array(0,c(q,nmax,m)),
		b=matrix(0,q,m),cvgd=integer(1),
		oxi=matrix(0,q,q),maxits=as.integer(maxits),
		llvec=numeric(as.integer(maxits)),
		eps=as.double(eps),
		#
		xiecme=matrix(0,q,q),
		g=g,reject=integer(maxits),
		ztvinvx=array(0,c(q,p,m)),
		a=array(0,c(q,q,m)),wkqp=matrix(0,q,p),
		wkg=rep(0,g+1),wkgg=matrix(0,g+1,g+1),
		abc=abc,dinv=dinv,PACKAGE="lmm")
	clock <- proc.time()-now
	cat("\n")
	iter <- tmp$iter
	cov.beta <- tmp$sigma2*tmp$xtwxinv
        msg <- tmp$msg
	{if(msg==1)
          warning("Supplied V <- i matrix is not positive definite")
	else if(msg==2)
          warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
	else if(msg==3)
          warning("Inadequate information to obtain starting value of psi")
	else if(msg==4)
          warning("Value of psi became non-pos.def. during iterations")
	else if(msg==5)
          warning("t(X)%*%W%*%X became non-pos.def. during iterations")
	else if(msg==10)
	  warning("non-pos.def. xi at input to scoring step")
	else if(msg==11)
	  warning("log-density not concave at one or more scoring steps")}
        llvec <- tmp$llvec[1:iter]
	reject <- tmp$reject[1:iter]
	converged <- tmp$cvgd==as.integer(1)
	if(!converged) warning(paste("did not converge by",
	   format(iter),"iterations"))
	#
        psi <- tmp$xi*tmp$sigma2
        beta <- tmp$beta
	if(!is.null(dimnames(pred)[[2]])){
           colnames <- dimnames(pred)[[2]]
           names(beta) <- colnames[xcol]
           dimnames(psi) <- list(colnames[zcol],colnames[zcol])}
        #
	list(beta=beta,sigma2=tmp$sigma2,psi=psi,
		converged=converged,iter=iter,reject=reject,
		logpost=llvec,cov.beta=cov.beta)}
###############################################################################
fastrml <- function(y,subj,pred,xcol,zcol,vmax,occ,start,maxits=50,
	eps=.0001){
        tmp <- table(subj)
	m <- length(tmp)
	nmax <- max(tmp)
	ntot <- length(y)
	pcol <- ncol(pred)
	q <- length(zcol)
	p <- length(xcol)
	g <- as.integer(round(q*(q+1)/2))
	#
	{if(missing(vmax)){
		vmax <- diag(rep(1,nmax))
                occ <- integer(ntot)
		iflag <- as.integer(1)}
	else iflag <- as.integer(0)}
	storage.mode(vmax) <- "double"
	#	
	{if(!missing(start)){
		beta <- start$beta
		sigma2 <- start$sigma2
		xi <- start$psi/start$sigma2
		storage.mode(beta) <- "double"
		storage.mode(xi) <- "double"
		storage.mode(sigma2) <- "double"
		sflag <- as.integer(1)}
	else{
		beta <- numeric(p)
		xi <- matrix(0,q,q)
		sigma2 <- 0
		sflag <- as.integer(0)}}
	cat("Performing FAST-RML...")
	now <- proc.time()
        err <- 0
	tmp <- .Fortran("fastrml",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,vmax,
		w=array(0,c(nmax,nmax,m)),
		vinv=array(0,c(nmax,nmax,m)),
		pcol,as.double(pred),q,as.integer(zcol),
		ztvinv=array(0,c(q,nmax,m)),
		ztvinvz=array(0,c(q,q,m)),
		iflag=iflag,err=as.integer(err),
		msg=integer(1),u=array(0,c(q,q,m)),
		iter=integer(1),sflag,sigma2=sigma2,p,
		as.integer(xcol),beta=beta,as.double(y),
		delta=rep(0,ntot),xtw=matrix(0,p,nmax),
		xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),wkqq1=matrix(0,q,q),
		wkqq2=matrix(0,q,q),
		xi=xi,wkqnm=array(0,c(q,nmax,m)),
		b=matrix(0,q,m),cvgd=integer(1),
		oxi=matrix(0,q,q),maxits=as.integer(maxits),
		llvec=numeric(as.integer(maxits)),
		eps=as.double(eps),xiecme=matrix(0,q,q),
		g=g,reject=integer(maxits),
		ztvinvx=array(0,c(q,p,m)),
		a=array(0,c(q,q,m)),wkqp=matrix(0,q,p),
		wkg=rep(0,g+1),wkgg=matrix(0,g+1,g+1),PACKAGE="lmm")
	clock <- proc.time()-now
	cat("\n")
	iter <- tmp$iter
        msg <- tmp$msg
	{if(msg==1)
          warning("Supplied V <- i matrix is not positive definite")
	else if(msg==2)
          warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
	else if(msg==3)
          warning("Inadequate information to obtain starting value of psi")
	else if(msg==4)
          warning("Value of psi became non-pos.def. during iterations")
	else if(msg==5)
          warning("t(X)%*%W%*%X became non-pos.def. during iterations")
	else if(msg==10)
	  warning("non-pos.def. xi at input to scoring step")
	else if(msg==11)
	  warning("loglikelihood not concave at one or more scoring steps")}
        llvec <- tmp$llvec[1:iter]
	reject <- tmp$reject[1:iter]
	converged <- tmp$cvgd==as.integer(1)
	if(!converged) warning(paste("did not converge by",
	   format(iter),"iterations"))
	b <- tmp$b
	u <- tmp$u
	a <- tmp$a
	xtwxinv <- tmp$xtwxinv
	sigma2 <- tmp$sigma2
	ztvinvx <- tmp$ztvinvx
	wkgg <- tmp$wkgg
	xi <- tmp$xi
	cov.beta <- tmp$sigma2*tmp$xtwxinv
	b.hat <- tmp$b
	cov.b <- tmp$sigma2*tmp$u
	tmp1 <- .Fortran("cebayes",m,p,q,b,u,a,xtwxinv,sigma2,ztvinvx,g,wkgg,
		wkqp=matrix(0,q,p),
		wkpg=matrix(0,p,g),wkqg=matrix(0,q,g),
		wkp=rep(0,p),varbeta=matrix(0,p,p),
		varb=array(0,c(q,q,m)),
		covbbeta=array(0,c(q,p,m)),
		wkpg2=matrix(0,p,g+1),wkqg2=matrix(0,q,g+1),
		xi=xi,wkqq1=matrix(0,q,q),ntot,err=integer(1),PACKAGE="lmm")
	{if(tmp1$err==0){
		cov.beta.new <- tmp1$varbeta
		cov.b.new <- tmp1$varb
		cov.b.beta.new <- tmp1$covbbeta}
	else{
		warning("loglikelihood not concave at solution")
		cov.b.beta.new <- NULL
		cov.b.new <- NULL
		cov.beta.new <- NULL}}
	#
        psi <- tmp$xi*tmp$sigma2
        beta <- tmp$beta
	if(!is.null(dimnames(pred)[[2]])){
           colnames <- dimnames(pred)[[2]]
           names(beta) <- colnames[xcol]
           dimnames(psi) <- list(colnames[zcol],colnames[zcol])}
        #
	list(beta=beta,sigma2=tmp$sigma2,psi=psi,
		converged=converged,iter=iter,reject=reject,
		loglik=llvec,cov.beta=cov.beta,b.hat=b.hat,cov.b=cov.b,
		cov.beta.new=cov.beta.new,cov.b.new=cov.b.new,
		cov.b.beta.new=cov.b.beta.new)}
###############################################################################
ecmerml <- function(y,subj,pred,xcol,zcol,vmax,occ,start,maxits=1000,
	eps=.0001){
        tmp <- table(subj)
	m <- length(tmp)
	nmax <- max(tmp)
	ntot <- length(y)
	pcol <- ncol(pred)
	q <- length(zcol)
	p <- length(xcol)
	#
	{if(missing(vmax)){
		vmax <- diag(rep(1,nmax))
                occ <- integer(ntot)
		iflag <- as.integer(1)}
	else iflag <- as.integer(0)}
	storage.mode(vmax) <- "double"
	#	
	{if(!missing(start)){
		beta <- start$beta
		sigma2 <- start$sigma2
		xi <- start$psi/start$sigma2
		storage.mode(beta) <- "double"
		storage.mode(xi) <- "double"
		storage.mode(sigma2) <- "double"
		sflag <- as.integer(1)}
	else{
		beta <- numeric(p)
		xi <- matrix(0,q,q)
		sigma2 <- 0
		sflag <- as.integer(0)}}
	cat("Performing ECME for RML estimation...")
	now <- proc.time()
        err <- 0
	tmp <- .Fortran("ecmerml",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,vmax,
		w=array(0,c(nmax,nmax,m)),
		vinv=array(0,c(nmax,nmax,m)),
		pcol,as.double(pred),q,as.integer(zcol),
		ztvinv=array(0,c(q,nmax,m)),
		ztvinvz=array(0,c(q,q,m)),
		iflag=iflag,err=as.integer(err),
		msg=integer(1),u=array(0,c(q,q,m)),
		iter=integer(1),sflag,sigma2=sigma2,p,
		as.integer(xcol),beta=beta,as.double(y),
		delta=rep(0,ntot),xtw=matrix(0,p,nmax),
		xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),wkqq1=matrix(0,q,q),
		wkqq2=matrix(0,q,q),
		xi=xi,wkqnm=array(0,c(q,nmax,m)),
		b=matrix(0,q,m),cvgd=integer(1),obeta=rep(0,p),
		oxi=matrix(0,q,q),maxits=as.integer(maxits),
		llvec=numeric(as.integer(maxits)),
		eps=as.double(eps),
		ztvinvx=array(0,c(q,p,m)),a=array(0,c(q,q,m)),
		wkqp=matrix(0,q,p),PACKAGE="lmm")
	clock <- proc.time()-now
	cat("\n")
	iter <- tmp$iter
        msg <- tmp$msg
	{if(msg==1)
          warning("Supplied V <- i matrix is not positive definite")
	else if(msg==2)
          warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
	else if(msg==3)
          warning("Inadequate information to obtain starting value of psi")
	else if(msg==4)
          warning("Value of psi became non-pos.def. during iterations")
	else if(msg==5)
          warning("t(X)%*%W%*%X became non-pos.def. during iterations")}
        llvec <- tmp$llvec[1:iter]
	converged <- tmp$cvgd==as.integer(1)
	cov.beta <- tmp$xtwxinv*tmp$sigma2
	b.hat <- tmp$b
	cov.b <- tmp$u*tmp$sigma2
	#
        psi <- tmp$xi*tmp$sigma2
        beta <- tmp$beta
	if(!is.null(dimnames(pred)[[2]])){
           colnames <- dimnames(pred)[[2]]
           names(beta) <- colnames[xcol]
           dimnames(psi) <- list(colnames[zcol],colnames[zcol])}
        #
	list(beta=beta,sigma2=tmp$sigma2,psi=psi,
		converged=converged,iter=iter,
		loglik=llvec,cov.beta=cov.beta,b.hat=b.hat,cov.b=cov.b)}
###############################################################################
mgibbs <- function(y,subj,pred,xcol,zcol,prior,seed,vmax,occ,start,
   iter=1000){
        tmp <- table(subj)
	m <- length(tmp)
	nmax <- max(tmp)
	ntot <- length(y)
	pcol <- ncol(pred)
	q <- length(zcol)
	p <- length(xcol)
	#
	{if(missing(vmax)){
		vmax <- diag(rep(1,nmax))
                occ <- integer(ntot)
		iflag <- as.integer(1)}
	else iflag <- as.integer(0)}
	storage.mode(vmax) <- "double"
	#	
	{if(!missing(start)){
		sigma2 <- start$sigma2
		xi <- start$psi/start$sigma2
		storage.mode(xi) <- "double"
		storage.mode(sigma2) <- "double"
		sflag <- as.integer(1)}
	else{
		xi <- matrix(0,q,q)
		sigma2 <- 0
		sflag <- as.integer(0)}}
	beta <- numeric(p)
	abc <- c(prior$a,prior$b,prior$c)
	storage.mode(abc) <- "double"
	dinv <- prior$Dinv
	storage.mode(dinv) <- "double"
	junk <- .Fortran("rngs",as.integer(seed),PACKAGE="lmm")
	cat("Performing modified Gibbs...")
	now <- proc.time()
        err <- 0
	tmp <- .Fortran("mgibbs",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,vmax,
		w=array(0,c(nmax,nmax,m)),
		vinv=array(0,c(nmax,nmax,m)),
		pcol,as.double(pred),q,as.integer(zcol),
		ztvinv=array(0,c(q,nmax,m)),
		ztvinvz=array(0,c(q,q,m)),
		isflag=c(iflag,sflag),err=as.integer(err),
		msg=integer(1),u=array(0,c(q,q,m)),
		sigma2=sigma2,p,
		as.integer(xcol),beta=beta,as.double(y),
		delta=rep(0,ntot),xtw=matrix(0,p,nmax),
		xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),wkqq1=matrix(0,q,q),
		wkqq2=matrix(0,q,q),
		xi=xi,wkqnm=array(0,c(q,nmax,m)),
		b=matrix(0,q,m),maxits=as.integer(iter),
		abc=abc,dinv=dinv,sqrtu=array(0,c(q,q,m)),
		sigma2s=rep(0,iter),psis=array(0,c(q,q,iter)),PACKAGE="lmm")
	clock <- proc.time()-now
	cat("\n")
        msg <- tmp$msg
	{if(msg==1)
          warning("Supplied V <- i matrix is not positive definite")
	else if(msg==2)
          warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
	else if(msg==3)
          warning("Inadequate information to obtain starting value of psi")
	else if(msg==4)
          warning("Value of psi became non-pos.def. during iterations")
	else if(msg==5)
          warning("t(X)%*%W%*%X became non-pos.def. during iterations")}
	#
        psi <- tmp$xi*tmp$sigma2
        beta <- tmp$beta
	if(!is.null(dimnames(pred)[[2]])){
           colnames <- dimnames(pred)[[2]]
           names(beta) <- colnames[xcol]
           dimnames(psi) <- list(colnames[zcol],colnames[zcol])}
        #
	list(beta=beta,sigma2=tmp$sigma2,psi=psi,
		sigma2.series=tmp$sigma2s,psi.series=tmp$psis)}
###############################################################################
fastml <- function(y,subj,pred,xcol,zcol,vmax,occ,start,maxits=50,
	eps=.0001){
        tmp <- table(subj)
	m <- length(tmp)
	nmax <- max(tmp)
	ntot <- length(y)
	pcol <- ncol(pred)
	q <- length(zcol)
	p <- length(xcol)
	g <- as.integer(round(q*(q+1)/2))
	#
	{if(missing(vmax)){
		vmax <- diag(rep(1,nmax))
                occ <- integer(ntot)
		iflag <- as.integer(1)}
	else iflag <- as.integer(0)}
	storage.mode(vmax) <- "double"
	#	
	{if(!missing(start)){
		beta <- start$beta
		sigma2 <- start$sigma2
		xi <- start$psi/start$sigma2
		storage.mode(beta) <- "double"
		storage.mode(xi) <- "double"
		storage.mode(sigma2) <- "double"
		sflag <- as.integer(1)}
	else{
		beta <- numeric(p)
		xi <- matrix(0,q,q)
		sigma2 <- 0
		sflag <- as.integer(0)}}
	cat("Performing FAST-ML...")
	now <- proc.time()
        err <- 0
	tmp <- .Fortran("fastml",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,vmax,
		w=array(0,c(nmax,nmax,m)),
		vinv=array(0,c(nmax,nmax,m)),
		pcol,as.double(pred),q,as.integer(zcol),
		ztvinv=array(0,c(q,nmax,m)),
		ztvinvz=array(0,c(q,q,m)),
		iflag=iflag,err=as.integer(err),
		msg=integer(1),u=array(0,c(q,q,m)),
		iter=integer(1),sflag,sigma2=sigma2,p,
		as.integer(xcol),beta=beta,as.double(y),
		delta=rep(0,ntot),xtw=matrix(0,p,nmax),
		xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),wkqq1=matrix(0,q,q),
		wkqq2=matrix(0,q,q),
		xi=xi,wkqnm=array(0,c(q,nmax,m)),
		b=matrix(0,q,m),cvgd=integer(1),
		oxi=matrix(0,q,q),maxits=as.integer(maxits),
		llvec=numeric(as.integer(maxits)),
		eps=as.double(eps),xiecme=matrix(0,q,q),
		g=g,reject=integer(maxits),
		wkg=rep(0,g+1),wkgg=matrix(0,g+1,g+1),PACKAGE="lmm")
	clock <- proc.time()-now
	cat("\n")
	iter <- tmp$iter
        msg <- tmp$msg
	{if(msg==1)
          warning("Supplied V <- i matrix is not positive definite")
	else if(msg==2)
          warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
	else if(msg==3)
          warning("Inadequate information to obtain starting value of psi")
	else if(msg==4)
          warning("Value of psi became non-pos.def. during iterations")
	else if(msg==5)
          warning("t(X)%*%W%*%X became non-pos.def. during iterations")
	else if(msg==10)
	  warning("non-pos.def. xi at input to scoring step")
	else if(msg==11)
	  warning("loglikelihood not concave at one or more scoring steps")}
        llvec <- tmp$llvec[1:iter]
	reject <- tmp$reject[1:iter]
	converged <- tmp$cvgd==as.integer(1)
	if(!converged) warning(paste("did not converge by",
	   format(iter),"iterations"))
	cov.beta <- tmp$xtwxinv*tmp$sigma2
	b.hat <- tmp$b
	cov.b <- tmp$sigma2*tmp$u
	#
        psi <- tmp$xi*tmp$sigma2
        beta <- tmp$beta
	if(!is.null(dimnames(pred)[[2]])){
           colnames <- dimnames(pred)[[2]]
           names(beta) <- colnames[xcol]
           dimnames(psi) <- list(colnames[zcol],colnames[zcol])}
        #
	list(beta=beta,sigma2=tmp$sigma2,psi=psi,
		converged=converged,iter=iter,reject=reject,
		loglik=llvec,cov.beta=cov.beta,b.hat=b.hat,cov.b=cov.b)}
###############################################################################
ecmeml <- function(y,subj,pred,xcol,zcol,vmax,occ,start,maxits=1000,
	eps=.0001){
        tmp <- table(subj)
	m <- length(tmp)
	nmax <- max(tmp)
	ntot <- length(y)
	pcol <- ncol(pred)
	q <- length(zcol)
	p <- length(xcol)
	#
	{if(missing(vmax)){
		vmax <- diag(rep(1,nmax))
                occ <- integer(ntot)
		iflag <- as.integer(1)}
	else iflag <- as.integer(0)}
	storage.mode(vmax) <- "double"
	#	
	{if(!missing(start)){
		beta <- start$beta
		sigma2 <- start$sigma2
		xi <- start$psi/start$sigma2
		storage.mode(beta) <- "double"
		storage.mode(xi) <- "double"
		storage.mode(sigma2) <- "double"
		sflag <- as.integer(1)}
	else{
		beta <- numeric(p)
		xi <- matrix(0,q,q)
		sigma2 <- 0
		sflag <- as.integer(0)}}
	cat("Performing ECME...")
	now <- proc.time()
        err <- 0
	tmp <- .Fortran("ecmeml",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,vmax,
		w=array(0,c(nmax,nmax,m)),
		vinv=array(0,c(nmax,nmax,m)),
		pcol,as.double(pred),q,as.integer(zcol),
		ztvinv=array(0,c(q,nmax,m)),
		ztvinvz=array(0,c(q,q,m)),
		iflag=iflag,err=as.integer(err),
		msg=integer(1),u=array(0,c(q,q,m)),
		iter=integer(1),sflag,sigma2=sigma2,p,
		as.integer(xcol),beta=beta,as.double(y),
		delta=rep(0,ntot),xtw=matrix(0,p,nmax),
		xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),wkqq1=matrix(0,q,q),
		wkqq2=matrix(0,q,q),
		xi=xi,wkqnm=array(0,c(q,nmax,m)),
		b=matrix(0,q,m),cvgd=integer(1),obeta=rep(0,p),
		oxi=matrix(0,q,q),maxits=as.integer(maxits),
		llvec=numeric(as.integer(maxits)),
		eps=as.double(eps),PACKAGE="lmm")
	clock <- proc.time()-now
	cat("\n")
	iter <- tmp$iter
        msg <- tmp$msg
	{if(msg==1)
          warning("Supplied V <- i matrix is not positive definite")
	else if(msg==2)
          warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
	else if(msg==3)
          warning("Inadequate information to obtain starting value of psi")
	else if(msg==4)
          warning("Value of psi became non-pos.def. during iterations")
	else if(msg==5)
          warning("t(X)%*%W%*%X became non-pos.def. during iterations")}
        llvec <- tmp$llvec[1:iter]
	converged <- tmp$cvgd==as.integer(1)
	cov.beta <- tmp$xtwxinv*tmp$sigma2
	b.hat <- tmp$b
	cov.b <- tmp$u*tmp$sigma2
	#
        psi <- tmp$xi*tmp$sigma2
        beta <- tmp$beta
	if(!is.null(dimnames(pred)[[2]])){
           colnames <- dimnames(pred)[[2]]
           names(beta) <- colnames[xcol]
           dimnames(psi) <- list(colnames[zcol],colnames[zcol])}
        #
	list(beta=beta,sigma2=tmp$sigma2,psi=psi,
		converged=converged,iter=iter,
		loglik=llvec,cov.beta=cov.beta,b.hat=b.hat,cov.b=cov.b)}
###############################################################################
