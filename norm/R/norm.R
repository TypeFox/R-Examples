#***********************************************************************
# Changes NA's to single precision missing value code
.na.to.snglcode <- function(x,mvcode){
      x[is.na(x)] <- as.double(mvcode)
  x}
#***********************************************************************
# Changes missing value code to NA
.code.to.na <- function(x,mvcode){
      x[x==mvcode] <- NA
  x}
#***********************************************************************
#  Perform preliminary manipulations on matrix of continuous data.  
#  Rows are sorted by missing data pattern.
prelim.norm <- function(x){
# get dimensions of x
  if(is.vector(x)) x <- matrix(x,length(x),1)
  n <- nrow(x); p <- ncol(x); storage.mode(x) <- "double"
# find missingness patterns
  r <- 1*is.na(x)
  nmis <- as.integer(apply(r,2,sum))
  names(nmis) <- dimnames(x)[[2]]
# index the missing data patterns
  mdp <- as.integer((r%*%(2^((1:ncol(x))-1)))+1)
# do row sort
  ro <- order(mdp)
  x <- matrix(x[ro,],n,p)
  mdp <- mdp[ro]
  r <- matrix(r[ro,],n,p)
  ro <- order(ro)
# compress missing data patterns
  mdpst <- as.integer(seq(along=mdp)[!duplicated(mdp)])
  mdp <- unique(mdp); npatt <- length(mdpst)
# create r-matrix for display purposes
  r <- 1-r; r <- matrix(r[mdpst,],npatt,p)
  if(npatt==1) tmp <- format(n)
  if(npatt>1)  tmp <- format(c(mdpst[2:npatt],n+1)-mdpst)
  dimnames(r) <- list(tmp,dimnames(x)[[2]])
  storage.mode(r) <- "integer"
# center and scale the columns of x
  if(sum(is.na(x))<length(x)){
    mvcode <- as.double(max(x[!is.na(x)])+1000)
    x <- .na.to.snglcode(x,mvcode)
    tmp <- .Fortran("ctrsc",x,n,p,numeric(p),numeric(p),mvcode,PACKAGE="norm")
    x <- tmp[[1]]; xbar <- tmp[[4]]; sdv <- tmp[[5]]
    x <- .code.to.na(x,mvcode)}
  if(sum(is.na(x))==length(x)){
    xbar <- rep(0,p); sdv <- rep(1,p)}  
# form matrix of packed storage indices
  d <- as.integer((2+3*p+p^2)/2)
  psi <- .Fortran("mkpsi",p,matrix(as.integer(0),p+1,p+1),PACKAGE="norm")[[2]]
# other bookkeeping quantities
  if(npatt>1) nmdp <- as.integer(c(mdpst[-1],n+1)-mdpst)
  if(npatt==1) nmdp <- n
  sj <- .Fortran("sjn",p,npatt,r,integer(p),PACKAGE="norm")[[4]]
  nmon <- .Fortran("nmons",p,npatt,r,nmdp,sj,integer(p),PACKAGE="norm")[[6]]
  last <- .Fortran("lasts",p,npatt,sj,integer(npatt),PACKAGE="norm")[[4]]
  tmp <- .Fortran("layers",p,sj,integer(p),integer(1),PACKAGE="norm")
  layer <- tmp[[3]]; nlayer <- tmp[[4]]
# return list
  list(x=x,n=n,p=p,r=r,nmis=nmis,ro=ro,mdpst=mdpst,
    nmdp=nmdp,npatt=npatt,xbar=xbar,sdv=sdv,d=d,psi=psi,sj=sj,
    nmon=nmon,last=last,layer=layer,nlayer=nlayer)}
#***********************************************************************
# Retrieves means and covariances from theta. If corr=FALSE , returns
# a list containing a vector of means and a covariance matrix. If
# corr=TRUE , returns a list containing a vector of means, a vector of
# standard deviations, and a correlation matrix.
getparam.norm <- function(s,theta,corr=FALSE ){
  mu <- theta[s$psi[1,2:(s$p+1)]]*s$sdv + s$xbar
  names(mu) <- dimnames(s$x)[[2]]
  sigma <- theta[s$psi[2:(s$p+1),2:(s$p+1)]]
  sigma <- matrix(sigma,s$p,s$p)
  tmp <- matrix(s$sdv,s$p,s$p)
  sigma <- sigma*tmp*t(tmp)
  dimnames(sigma) <- list(names(mu),names(mu))
  if(corr){
    sdv <- sqrt(diag(sigma)); names(sdv) <- names(mu)
    tmp <- matrix(sdv,s$p,s$p)
    r <- sigma/(tmp*t(tmp)); dimnames(r) <- list(names(mu),names(mu))
    result <- list(mu=mu,sdv=sdv,r=r)}
  else result <- list(mu=mu,sigma=sigma)
  result}
#***********************************************************************
# Makes a theta vector out of a list of specified parameters.
makeparam.norm <- function(s,thetalist){
  result <- numeric(s$d); result[1] <- -1
  xbar <- s$xbar;sdv <- s$sdv
  mu <- (thetalist[[1]]-xbar)/sdv
  result[2:(s$p+1)] <- mu
  if(length(thetalist)==3){
    tmp <- matrix(thetalist[[2]],s$p,s$p)
    sigma <- thetalist[[3]]*tmp*t(tmp)}
  else sigma <- thetalist[[2]]
  tmp <- matrix(sdv,s$p,s$p)
  sigma <- sigma/(tmp*t(tmp))
  tmp <- as.vector(s$psi[2:(s$p+1),2:(s$p+1)])
  result[tmp] <- as.vector(sigma)
  result}
#***********************************************************************
# Finds posterior mode of theta under the multivariate
# normal model. If no prior is specified, finds the mle.
em.norm <- function(s,start,showits=TRUE ,maxits=1000,criterion=.0001,
     prior){
  s$x <- .na.to.snglcode(s$x,999)
  if(missing(start)){
    start <- .Fortran("stvaln",s$d,numeric(s$d),s$p,s$psi,PACKAGE="norm")[[2]]}
  if(missing(prior)){
    mle <- as.integer(1)
    tau <- numeric(1); m <- numeric(1); mu0 <- numeric(s$p);
    lambdainv <- matrix(0,s$p,s$p)}
  if(!(missing(prior))){
    mle <- as.integer(0)
    tau <- as.numeric(prior[[1]]); m <- as.numeric(prior[[2]])
    mu0 <- as.numeric(prior[[3]])
    lambdainv <- as.numeric(prior[[4]])}
  tmp <- as.integer(numeric(s$p))
  tobs <- .Fortran("tobsn",s$d,numeric(s$d),s$p,s$psi,s$n,s$x,s$npatt,
    s$r,s$mdpst,s$nmdp,tmp,PACKAGE="norm")[[2]]
# iterate to mle
  it <- 0; converged <- FALSE
  if(showits) cat(paste("Iterations of EM:","\n"))
  while((!converged)&(it<maxits)){
  old <- start
  start <- .Fortran("emn",s$d,old,start,tobs,s$p,s$psi,s$n,
    s$x,s$npatt,s$r,s$mdpst,s$nmdp,tmp,tmp,numeric(s$p),
    mle,tau,m,mu0,lambdainv,PACKAGE="norm")[[3]]
# print iteration number
  it <- it+1; if(showits) cat(paste(format(it),"...",sep=""))
  converged <- max(abs(old-start))<=criterion}
  if(showits)cat("\n")
  start}       
#***********************************************************************
# Calculates log observed-data posterior at theta
logpost.norm <- function(s,theta,prior){
  s$x <- .na.to.snglcode(s$x,999)
  l1 <- .Fortran("lobsn",s$d,theta,numeric(s$d),s$p,s$psi,s$n,s$x,
    s$npatt,s$r,s$mdpst,s$nmdp,as.integer(numeric(s$p)),numeric(s$p),
    0,PACKAGE="norm")[[14]]
  if(!(missing(prior))){
  tau <- as.numeric(prior[[1]]); m <- as.numeric(prior[[2]])
  mu0 <- as.numeric(prior[[3]]); lambdainv <- as.numeric(prior[[4]])}
  if(missing(prior)){
    tau <- as.numeric(0); m <- as.numeric(-1); mu0 <- numeric(s$p)
    lambdainv <- matrix(0,s$p,s$p)}
  l2 <- .Fortran("lprin",s$d,theta,s$p,s$psi,numeric(s$p),tau,m,mu0,
    lambdainv,0,PACKAGE="norm")[[10]]
  l1+l2}
#***********************************************************************
# Calculates observed-data loglikelihood at theta
loglik.norm <- function(s,theta){
  s$x <- .na.to.snglcode(s$x,999)
  .Fortran("lobsn",s$d,theta,numeric(s$d),s$p,s$psi,s$n,s$x,s$npatt,
    s$r,s$mdpst,s$nmdp,as.integer(numeric(s$p)),numeric(s$p),0,PACKAGE="norm")[[14]]}
#***********************************************************************
# Simulate a value of theta from a normal-inverted Wishart
# distribution
ninvwish <- function(s,params){
    tau <- as.numeric(params[[1]]); m <- as.numeric(params[[2]])
    mu0 <- params[[3]]; lambdainv <- params[[4]]
    tmpi <- as.integer(numeric(s$p)); tmpr <- as.double((numeric(s$p)))
    pri <- numeric(s$d); pri[2:(s$p+1)] <- mu0
    tmp <- as.vector(s$psi[2:(s$p+1),2:(s$p+1)])
    pri[tmp] <- as.vector(lambdainv)
    .Fortran("ninvwn",s$d,pri,tau,m,s$p,s$psi,
       numeric((s$p)^2),tmpr,numeric(s$d),tmpi,PACKAGE="norm")[[2]]}
#***********************************************************************
# Data augmentation for the multivariate normal.
# Produces a new draw of theta from its posterior distribution. 
da.norm <- function(s,start,prior,steps=1,showits=FALSE ,return.ymis=FALSE ){
  s$x <- .na.to.snglcode(s$x,999)
  tmpi <- as.integer(numeric(s$p)); tmpr <- as.double((numeric(s$p)))
  tobs <- .Fortran("tobsn",s$d,numeric(s$d),s$p,s$psi,s$n,s$x,s$npatt,
    s$r,s$mdpst,s$nmdp,tmpi,PACKAGE="norm")[[2]]
  if(missing(prior)){
    tau <- 0; m <- -1; mu0 <- numeric(s$p);lambdainv <- matrix(0,s$p,s$p)}
  if(!(missing(prior))){
    tau <- as.numeric(prior[[1]]); m <- as.numeric(prior[[2]])
    mu0 <- prior[[3]]; lambdainv <- prior[[4]]}
  pri <- numeric(s$d); pri[2:(s$p+1)] <- mu0
  tmp <- as.vector(s$psi[2:(s$p+1),2:(s$p+1)])
  pri[tmp] <- as.vector(lambdainv)
  if(showits) cat(paste("Steps of Data Augmentation:","\n"))
  for(i in 1:steps){
    if(showits) cat(paste(format(i),"...",sep=""))
    tmp <- .Fortran("is1n",s$d,start,start,tobs,s$p,s$psi,s$n,s$x,
      s$npatt,s$r,s$mdpst,s$nmdp,tmpi,tmpi,tmpr,start,PACKAGE="norm")
    start <- tmp[[3]]
    start <- .Fortran("ps1n",s$d,start,m,tau,pri,s$p,s$psi,s$n,
      matrix(0,s$p,s$p),tmpr,start,tmpi,PACKAGE="norm")[[5]]}
  if(showits)cat("\n")
  if(return.ymis){
    ymis <- tmp[[8]]*matrix(s$sdv,s$n,s$p,TRUE)+matrix(s$xbar,s$n,s$p,TRUE)
    ymis <- ymis[s$ro,];ymis <- ymis[s$x[s$ro,]==999]
    start <- list(parameter=start,ymis=ymis)}
  start}
#***********************************************************************
# Generates a single imputed dataset under theta
imp.norm <- function(s,theta,x){
  s$x <- .na.to.snglcode(s$x,999)
  tmpi <- as.integer(numeric(s$p)); tmpr <- as.double((numeric(s$p)))
  tobs <- .Fortran("tobsn",s$d,numeric(s$d),s$p,s$psi,s$n,s$x,s$npatt,
    s$r,s$mdpst,s$nmdp,tmpi,PACKAGE="norm")[[2]]
  s$x <- .Fortran("is1n",s$d,theta,theta,tobs,s$p,s$psi,s$n,s$x,
      s$npatt,s$r,s$mdpst,s$nmdp,tmpi,tmpi,tmpr,theta,PACKAGE="norm")[[8]]
  s$x <- s$x*matrix(s$sdv,s$n,s$p,TRUE)+matrix(s$xbar,s$n,s$p,TRUE)
  s$x <- s$x[s$ro,]
  if(!missing(x))x[is.na(x)] <- s$x[is.na(x)]
  else{x <- s$x; storage.mode(x) <- "double"}
  x}
#***********************************************************************
# Monotone data augmentation for the multivariate normal.
# Produces a new draw of theta from its posterior distribution. 
mda.norm <- function(s,theta,steps=1,showits=FALSE ){
  s$x <- .na.to.snglcode(s$x,999)
  tobs <- .Fortran("tobsmn",s$p,s$psi,s$n,s$x,s$npatt,s$r,s$mdpst,s$nmdp,
    s$last,integer(s$p),s$sj,s$layer,s$nlayer,s$d,
    matrix(0,s$nlayer,s$d),PACKAGE="norm")[[15]]
  if(showits) cat(paste("Steps of Monotone Data Augmentation:",
    "\n")) 
  for(i in 1:steps){
    if(showits) cat(paste(format(i),"...",sep=""))
    s$x <- .Fortran("is2n",s$d,theta,s$p,s$psi,s$n,s$x,s$npatt,
      s$r,s$mdpst,s$nmdp,s$sj,s$last,integer(s$p),integer(s$p),
      double(s$p),theta,PACKAGE="norm")[[6]]
    theta <- .Fortran("ps2n",s$p,s$psi,s$n,s$x,s$npatt,s$r,s$mdpst,
      s$nmdp,integer(s$p),integer(s$p),s$nmon,s$sj,s$nlayer,s$d,
      tobs,numeric(s$d),numeric(s$d),numeric(s$p+1),numeric(s$d),PACKAGE="norm")[[19]]}
  if(showits)cat("\n")
  theta}
#***********************************************************************
# Initializes random number generator seed. Argument should be a 
# positive integer
rngseed <- function(seed){
  seed <- as.integer(seed)
  if(seed<=0)stop("Seed must be a positive integer")
  tmp <- .Fortran("rngs",seed,PACKAGE="norm")
  invisible()}
#***********************************************************************
# multiple imputation inference
mi.inference <- function(est,std.err,confidence=.95){
  qstar <- est[[1]]
  for(i in 2:length(est)){qstar <- cbind(qstar,est[[i]])}
  qbar <- apply(qstar,1,mean)
  u <- std.err[[1]]
  for(i in 2:length(std.err)){u <- cbind(u,std.err[[i]])}
  dimnames(u)[[1]] <- dimnames(qstar)[[1]]
  u <- u^2
  ubar <- apply(u,1,mean)
  bm <- apply(qstar,1,var)
  m <- dim(qstar)[2]
  tm <- ubar+((1+(1/m))*bm)
  rem <- (1+(1/m))*bm/ubar
  nu <- (m-1)*(1+(1/rem))**2
  alpha <- 1-(1-confidence)/2
  low <- qbar-qt(alpha,nu)*sqrt(tm)
  up <- qbar+qt(alpha,nu)*sqrt(tm)
  pval <- 2*(1-pt(abs(qbar/sqrt(tm)),nu))
  fminf <- (rem+2/(nu+3))/(rem+1)
  result <- list(est=qbar,std.err=sqrt(tm),df=nu,signif=pval,lower=low,
  upper=up,r=rem,fminf=fminf)
  result}
#***********************************************************************

# Added function:

# ".First.lib" <-
# function(lib, pkg) {
#    library.dynam("norm", pkg, lib) }
