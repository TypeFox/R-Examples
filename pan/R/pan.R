###############################################################################
# A new version of pan that allows psi to be block-diagonal.
# This looks the same as the original version, except that psi is
# now an array of dimension (q x q x r). Also, the prior
# hyperparameters for psi, c and Dinv, have a new form;  c is now
# a vector of length r, and Dinv is now a (q x q x r) array.
pan.bd <- function(y,subj,pred,xcol,zcol,prior,seed,iter=1,start){
	if(any(is.na(pred)))
		stop("missing values in pred not allowed")
        # change y and pred to matrices, if necessary
	if(is.vector(y)) y <- matrix(y,ncol=1)
	if(is.vector(pred)) pred <- matrix(pred,ncol=1)
	# change prior$Binv and Dinv to matrices if necessary
	if(length(prior$Binv)==1) prior$Binv <- matrix(prior$Binv,1,1)
	if(length(prior$Dinv)==1) prior$Dinv <- matrix(prior$Binv,1,1)
	m <- as.integer(length(table(subj)))
	ntot <- as.integer(nrow(y))
	r <- as.integer(ncol(y))
	p <- length(xcol)
	q <- length(zcol)
	zcol <- as.integer(zcol)
	xcol <- as.integer(xcol)
	pcol <- as.integer(ncol(pred))
	#
	{if(missing(start)){
		beta <- matrix(0,p,r)
		sigma <- matrix(0,r,r)
		psi <- array(0,c(q,q,r))
		eps <- matrix(0,ntot,r)
		sflag <- as.integer(0)}
	 else{
		beta <- start$beta
		sigma <- start$sigma
		psi <- start$psi
		eps <- start$y
		sflag <- as.integer(1)
		storage.mode(eps) <- "double"
		storage.mode(beta) <- "double"
		storage.mode(sigma) <- "double"
		storage.mode(psi) <- "double"}}
	# This part packs the hyperparameters a, Binv, c, Dinv into a vector
	# of length 1 + r*(r+1)/2 + r + r*q*(q+1)/2
	nhyp <- as.integer(1 + r*(r+1)/2 + r + r*q*(q+1)/2)
	hyp <- rep(0,nhyp)
	hyp[1] <- prior$a
	k <- 1
	for(i in 1:r){
		for(j in i:r){
			k <- k+1
			hyp[k] <- prior$Binv[i,j]}}
	for(i in 1:r){
		k <- k+1
		hyp[k] <- prior$c[i]}
	for(l in 1:r){
		for(i in 1:q){
			for(j in i:q){
				k <- k+1
				hyp[k] <- prior$Dinv[i,j,l]}}}
	storage.mode(hyp) <- "double"
	# create rmat, npatt and patt to keep track of missingness patterns
	rmat <- 1-1*is.na(y)
	storage.mode(rmat) <- "integer"
	revcpatt <- rep("",ntot)
        for(i in 1:r) revcpatt <- paste(as.character(rmat[,i]),revcpatt,sep="")
	nulpat0 <- ""
        nulpat2 <- ""
        for(i in 1:r){
           nulpat0 <- paste(nulpat0,"0",sep="")
           nulpat2 <- paste(nulpat2,"2",sep="")}
        revcpatt[revcpatt==nulpat0] <- nulpat2
	tmp <- rev(table(revcpatt))
	npatt <- length(tmp)
	if(any(revcpatt==nulpat2)) npatt <- npatt-1
	w <- !duplicated(revcpatt)
	upatt <- revcpatt[w]
	rmat <- rmat[w,]
        if(r==1) rmat <- matrix(rmat,ncol=1)
	w <- rev(order(upatt))
	upatt <- upatt[w]
	rmat <- rmat[w,]
        if(r==1) rmat <- matrix(rmat,ncol=1)
	if(any(upatt==nulpat2)){
		rmat <- rmat[-1,]
		upatt <- upatt[-1]}
        patt <- integer(ntot)
	patt[revcpatt==nulpat2] <- 0
	for(i in 1:npatt) patt[revcpatt==upatt[i]] <- i
	storage.mode(npatt) <- "integer"
	storage.mode(rmat) <- "integer"
	storage.mode(patt) <- "integer"
	iposn <- as.integer(1:ntot)
	w <- order(patt)
	iposn <- iposn[w]
	pstfin <- matrix(0,npatt,2)
	{if(any(patt==0)){
		st <- tmp[1]+1
		for(i in 2:(npatt+1)){
			pstfin[i-1,1] <- st
			pstfin[i-1,2] <- st+tmp[i]-1
			st <- st+tmp[i]}}
	else{
		st <- 1
		for(i in 1:npatt){
			pstfin[i,1] <- st
			pstfin[i,2] <- st+tmp[i]-1
			st <- st+tmp[i]}}}
	storage.mode(pstfin) <- "integer"
	#
	storage.mode(y) <- "double"
	y[is.na(y)] <- -999.99
	storage.mode(pred) <- "double"
	#
	.Fortran("rngs",as.integer(seed),PACKAGE="pan")
	#
	tmp <- .Fortran("mgibbsbd",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),pcol,pred,
		q,zcol,ztz=array(0,c(q,q,m)),patt,nstar=integer(1),
		r,y=y,p,xcol,npatt,rmat,sflag,
		beta=beta,sigma=sigma,psi=psi,
		b=array(0,c(q,r,m)),xtxinv=matrix(0,p,p),wkpp=matrix(0,p,p),
		wkpr=matrix(0,p,r),eps=eps,
		wkrqrq1=matrix(0,r*q,r*q),wkrqrq2=matrix(0,r*q,r*q),
		sig=array(0,c(r*q,r*q,m)),wkrr1=matrix(0,r,r),
		wkrr2=matrix(0,r,r),
		iter=as.integer(iter),wkqr1=matrix(0,q,r),wkqr2=matrix(0,q,r),
		wkqrv=numeric(q*r),nhyp=nhyp,hyp=hyp,delta=matrix(0,ntot,r),
		iposn=iposn,pstfin=pstfin,
		betas=array(0,c(p,r,iter)),sigmas=array(0,c(r,r,iter)),
		psis=array(0,c(q,q,r,iter)),
		wkqq1=matrix(0,q,q),wkqq2=matrix(0,q,q),PACKAGE="pan")
	# restore y to being a vector if necessary
	if(r==1) tmp$y <- as.vector(tmp$y)
	#
	list(beta=tmp$betas,sigma=tmp$sigmas,psi=tmp$psis,
		y=tmp$y,last=list(beta=tmp$beta,sigma=tmp$sigma,
		psi=tmp$psi,y=tmp$y))}
###############################################################################
pan <- function(y,subj,pred,xcol,zcol,prior,seed,iter=1,start){
	if(any(is.na(pred)))
		stop("missing values in pred not allowed")
        # change y and pred to matrices, if necessary
	if(is.vector(y)) y <- matrix(y,ncol=1)
	if(is.vector(pred)) pred <- matrix(pred,ncol=1)
	# change prior$Binv and Dinv to matrices if necessary
	if(length(prior$Binv)==1) prior$Binv <- matrix(prior$Binv,1,1)
	if(length(prior$Dinv)==1) prior$Dinv <- matrix(prior$Binv,1,1)
	m <- as.integer(length(table(subj)))
	ntot <- as.integer(nrow(y))
	r <- as.integer(ncol(y))
	p <- length(xcol)
	q <- length(zcol)
	zcol <- as.integer(zcol)
	xcol <- as.integer(xcol)
	pcol <- as.integer(ncol(pred))
	#
	{if(missing(start)){
		beta <- matrix(0,p,r)
		sigma <- matrix(0,r,r)
		psi <- matrix(0,q*r,q*r)
		eps <- matrix(0,ntot,r)
		sflag <- as.integer(0)}
	 else{
		beta <- start$beta
		sigma <- start$sigma
		psi <- start$psi
		eps <- start$y
		sflag <- as.integer(1)
		storage.mode(eps) <- "double"
		storage.mode(beta) <- "double"
		storage.mode(sigma) <- "double"
		storage.mode(psi) <- "double"}}
	#
	nhyp <- as.integer(2+(r*(r+1)/2)+(r*q*(r*q+1)/2))
	hyp <- rep(0,nhyp)
	hyp[1] <- prior$a
	k <- 1
	for(i in 1:r){
		for(j in i:r){
			k <- k+1
			hyp[k] <- prior$Binv[i,j]}}
	k <- k+1
	hyp[k] <- prior$c
	for(i in 1:(r*q)){
		for(j in i:(r*q)){
			k <- k+1
			hyp[k] <- prior$Dinv[i,j]}}
	storage.mode(hyp) <- "double"
	# create rmat, npatt and patt to keep track of missingness patterns
	rmat <- 1-1*is.na(y)
	storage.mode(rmat) <- "integer"
	revcpatt <- rep("",ntot)
        for(i in 1:r) revcpatt <- paste(as.character(rmat[,i]),revcpatt,sep="")
	nulpat0 <- ""
        nulpat2 <- ""
        for(i in 1:r){
           nulpat0 <- paste(nulpat0,"0",sep="")
           nulpat2 <- paste(nulpat2,"2",sep="")}
        revcpatt[revcpatt==nulpat0] <- nulpat2
	tmp <- rev(table(revcpatt))
	npatt <- length(tmp)
	if(any(revcpatt==nulpat2)) npatt <- npatt-1
	w <- !duplicated(revcpatt)
	upatt <- revcpatt[w]
	rmat <- rmat[w,]
        if(r==1) rmat <- matrix(rmat,ncol=1)
	w <- rev(order(upatt))
	upatt <- upatt[w]
	rmat <- rmat[w,]
        if(r==1) rmat <- matrix(rmat,ncol=1)
	if(any(upatt==nulpat2)){
		rmat <- rmat[-1,]
		upatt <- upatt[-1]}
        patt <- integer(ntot)
	patt[revcpatt==nulpat2] <- 0
	for(i in 1:npatt) patt[revcpatt==upatt[i]] <- i
	storage.mode(npatt) <- "integer"
	storage.mode(rmat) <- "integer"
	storage.mode(patt) <- "integer"
	iposn <- as.integer(1:ntot)
	w <- order(patt)
	iposn <- iposn[w]
	pstfin <- matrix(0,npatt,2)
	{if(any(patt==0)){
		st <- tmp[1]+1
		for(i in 2:(npatt+1)){
			pstfin[i-1,1] <- st
			pstfin[i-1,2] <- st+tmp[i]-1
			st <- st+tmp[i]}}
	else{
		st <- 1
		for(i in 1:npatt){
			pstfin[i,1] <- st
			pstfin[i,2] <- st+tmp[i]-1
			st <- st+tmp[i]}}}
	storage.mode(pstfin) <- "integer"
	#
	storage.mode(y) <- "double"
	y[is.na(y)] <- -999.99
	storage.mode(pred) <- "double"
	#
	.Fortran("rngs",as.integer(seed),PACKAGE="pan")
	#
	tmp <- .Fortran("mgibbs",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),pcol,pred,
		q,zcol,ztz=array(0,c(q,q,m)),patt,nstar=integer(1),
		r,y=y,p,xcol,npatt,rmat,sflag,
		beta=beta,sigma=sigma,psi=psi,
		b=array(0,c(q,r,m)),xtxinv=matrix(0,p,p),wkpp=matrix(0,p,p),
		wkpr=matrix(0,p,r),eps=eps,
		wkrqrq1=matrix(0,r*q,r*q),wkrqrq2=matrix(0,r*q,r*q),
		sig=array(0,c(r*q,r*q,m)),wkrr1=matrix(0,r,r),
		wkrr2=matrix(0,r,r),
		iter=as.integer(iter),wkqr1=matrix(0,q,r),wkqr2=matrix(0,q,r),
		wkqrv=numeric(q*r),nhyp,hyp,delta=matrix(0,ntot,r),
		iposn=iposn,pstfin=pstfin,
		betas=array(0,c(p,r,iter)),sigmas=array(0,c(r,r,iter)),
		psis=array(0,c(r*q,r*q,iter)),PACKAGE="pan")
	# restore y to being a vector if necessary
	if(r==1) tmp$y <- as.vector(tmp$y)
	#
	list(beta=tmp$betas,sigma=tmp$sigmas,psi=tmp$psis,
		y=tmp$y,last=list(beta=tmp$beta,sigma=tmp$sigma,
		psi=tmp$psi,y=tmp$y))}
###############################################################################
ecme <- function(y,subj,occ,pred,xcol,zcol=NULL,vmax,start,maxits=1000,
	eps=.0001,random.effects=F){
	# Calculates ML estimates for the univariate linear mixed model
	# using the ECME algorithm described by Schafer (1997).
	m <- length(table(subj))
	ntot <- length(y)
	nmax <- as.integer(max(occ))
	pcol <- ncol(pred)
	q <- length(zcol)
	p <- length(xcol)
	#
	{if(missing(vmax)){
		vmax <- diag(rep(1,nmax))
		iflag <- as.integer(1)}
	else iflag <- as.integer(0)}
	storage.mode(vmax) <- "double"
	#	
	if(is.null(zcol)){
	q <- as.integer(1)
	tmp <- .Fortran("prelim",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,as.double(vmax),
		vh=array(0,c(nmax,nmax,m)),vi=array(0,c(nmax,nmax,m)),pcol,
		as.double(pred),q,
		integer(1),ztv=array(0,c(q,nmax,m)),
		sig=array(0,c(q,q,m)),iflag,PACKAGE="pan")
		ist <- tmp$ist
		ifin <- tmp$ifin
		vh <- tmp$vh
		vi <- tmp$vi
		ztv <- tmp$ztv
		sig0 <- tmp$sig
	tmp <- .Fortran("nopsi",ntot,m,ist,ifin,as.integer(occ),nmax,
		vi,vh,pcol,as.double(pred),q,integer(1),ztv,
		sig0,iflag,sig=array(0,c(q,q,m)),psi=matrix(0,q,q),
		sigma2=numeric(1),
		p,as.integer(xcol),beta=numeric(p),matrix(0,q,q),
		matrix(0,q,q),
		matrix(0,q,q),as.double(y),delta=numeric(ntot),
		b=matrix(0,q,m),
		wk=array(0,c(q,nmax,m)),w=array(0,c(nmax,nmax,m)),
		xtw=matrix(0,p,nmax),
	        xtwx=matrix(0,p,p),xtwy=numeric(p),xtwxinv=matrix(0,p,p),
		ll=numeric(1),PACKAGE="pan")
	beta <- tmp$beta
	sigma2 <- tmp$sigma2
	#llvec <- tmp$ll  I'm not sure that this loglikelihood is right,
	#              so I'll omit it for now.
	llvec <- NULL
	cov.beta <- tmp$xtwxinv
	converged <- T
	iter <- as.integer(1)
	psi <- NULL
	bhat <- NULL	
	cov.b <- NULL
	}
	#
	else{
	{if(!missing(start)){
		beta <- start$beta
		sigma2 <- start$sigma2
		psi <- start$psi
		storage.mode(beta) <- "double"
		storage.mode(psi) <- "double"
		storage.mode(sigma2) <- "double"
		sflag <- as.integer(1)}
	else{
		beta <- numeric(p)
		psi <- matrix(0,q,q)
		sigma2 <- 0
		sflag <- as.integer(0)}}
	llvec <- numeric(as.integer(maxits))
	cat("Performing ECME...")
	tmp <- .Fortran("ecme3",ntot,as.integer(subj),m,ist=integer(m),
		ifin=integer(m),as.integer(occ),nmax,
		vi=array(0,c(nmax,nmax,m)),
		vh=array(0,c(nmax,nmax,m)),
		pcol,as.double(pred),q,as.integer(zcol),
		ztv=array(0,c(q,nmax,m)),
		sig0=array(0,c(q,q,m)),iflag=iflag,
		sig=array(0,c(q,q,m)),psi=psi,
		sigma2=sigma2,
		p,as.integer(xcol),beta=beta,
		matrix(0,q,q),matrix(0,q,q),
		matrix(0,q,q),as.double(y),delta=rep(0,ntot),
		b=matrix(0,q,m),wk=array(0,c(q,nmax,m)),
		w=array(0,c(nmax,nmax,m)),xtw=matrix(0,p,nmax),
	        xtwx=matrix(0,p,p),xtwy=numeric(p),
		xtwxinv=matrix(0,p,p),llvec=llvec,vmax,
		sflag=sflag,eps=as.double(eps),
		obeta=rep(0,p),opsi=matrix(0,q,q),
		maxits=as.integer(maxits),
		iter=integer(1),cvgd=integer(1),PACKAGE="pan")
	cat("\n")
	iter <- tmp$iter
        llvec <- tmp$llvec[1:iter]
	converged <- tmp$cvgd==as.integer(1)
	cov.beta <- tmp$xtwxinv*tmp$sigma2
	cov.beta <- cov.beta+t(cov.beta)-diag(cov.beta)
	{if(random.effects){
		bhat <- tmp$b
		cov.b <- tmp$sig
		cov.b <- .Fortran("bdiag",q,m,cov.b=cov.b,PACKAGE="pan")$cov.b}
	else{
		bhat <- NULL
		cov.b <- NULL}}
	}
	#
	list(beta=tmp$beta,sigma2=tmp$sigma2,psi=tmp$psi,
		converged=converged,iter=iter,
		loglik=llvec,cov.beta=cov.beta,bhat=bhat,cov.b=cov.b)}
