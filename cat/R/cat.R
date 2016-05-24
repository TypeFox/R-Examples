#************************************************************************
#  Perform preliminary manipulations on matrix of categorical data.  
#  Rows are sorted by missing data pattern and by the values of the
#  observed variables within each pattern. Data must be coded as
#  positive integers.
prelim.cat <- function(x,counts,levs){
# get dimensions of x
  if(is.vector(x)) x <- matrix(x,length(x),1)
  p <- ncol(x); n <- nrow(x)
####################################################
# Necessary to reverse the order of 'missing tests #
# since in R once counts has been assigned to it   #
# is no longer missing            EFH 17 June 2003 #
####################################################
  if(!missing(counts)){grouped <- TRUE}
  if(missing(counts)){
     counts <- rep(1,nrow(x))
     grouped <- FALSE}
##if(!missing(counts)){grouped <- T}
# get dimensions of contingency table
####################################################
# Not necessary to reverse order of 'missing tests #
# since levs is not been assigned to here          #
####################################################
  if(missing(levs)){
    d <- as.integer(numeric(p))
    for(j in 1:p){d[j] <- as.integer(max(x[!is.na(x[,j]),j]))}}
  if(!(missing(levs))){d <- as.integer(levs)}
# get response indicators
  r <- 1*is.na(x)
  nmis <- as.integer(apply(r,2,sum))
  names(nmis) <- dimnames(x)[[2]]
# index the missing data patterns
  mdp <- as.integer((r%*%(2^((1:ncol(x))-1)))+1)
# calculate the known part of the cell number for rows of x
  cumd <- as.integer(cumprod(d))
  ncells <- cumd[p]; jmp <- as.integer(cumd/d)
  r <- 1-r
  tmp <- x;  tmp[is.na(tmp)] <- as.integer(1)	
  mobs0 <- as.integer(1+((tmp-1)*r)%*%(cumd/d))
# do row sort
  ro <- order(mdp,mobs0)
############################################################
# The following line changed: as.matrix has no params in R #
# EFH 16/06/2003                                           #
############################################################
##x <- as.matrix(x[ro,],n,p)
  x <- matrix(x[ro,],n,p)
  mdp <- mdp[ro]; mobs0 <- mobs0[ro]
############################################################
# The following line changed: as.matrix has no params in R #
# EFH 17/06/2003                                           #
############################################################
##r <- as.matrix(r[ro,],n,p)
  r <- matrix(r[ro,],n,p)
  counts <- counts[ro]
  ro <- order(ro)
  storage.mode(x) <- "integer"
# compress missing data patterns
  mdpst <- as.integer(seq(along=mdp)[!duplicated(mdp)])
  mdp <- unique(mdp); npatt <- length(mdp)
  if(npatt>1) nmdp <- as.integer(c(mdpst[-1],n+1)-mdpst)
  if(npatt==1) nmdp <- n
# group the data within missing data patterns
  mdpgrp <- numeric(npatt)
  for(j in 1:npatt){
    mdpgrp[j] <- length(table(mobs0[mdpst[j]:(mdpst[j]+nmdp[j]-1)]))}
  mdpgrp <- as.integer(mdpgrp)
  ngrp <- as.integer(sum(mdpgrp))
  if(ngrp==1) mdpgst <- as.integer(1)
  if(ngrp>1) mdpgst <- as.integer(cumsum(c(1,mdpgrp[-npatt])))
  mobs <- numeric(ngrp);mobsst <- numeric(ngrp);nmobs <- numeric(ngrp)
  for(i in 1:npatt){
    w <- mdpst[i]:(mdpst[i]+nmdp[i]-1)
    g <- mdpgst[i]:(mdpgst[i]+mdpgrp[i]-1)
    mobs[g] <- unique(mobs0[w])
    mobsst[g] <- w[!duplicated(mobs0[w])]
    for(j in g){nmobs[j] <- sum((counts[w])[mobs0[w]==mobs[j]])}}
  storage.mode(mobs) <- "integer"
  storage.mode(mobsst) <- "integer"
  storage.mode(nmobs) <- "integer"
# create r-matrix for display purposes
  r <- matrix(r[mdpst,],npatt,p)
  if(!grouped){
    if(npatt==1) tmp <- format(n)
    if(npatt>1)  tmp <- format(c(mdpst[2:npatt],n+1)-mdpst)}
  if(grouped){
    tmp <- numeric(npatt)
    for(i in 1:npatt){
      w <- mdpgst[i]:(mdpgst[i]+mdpgrp[i]-1)
      tmp[i] <- sum(nmobs[w])
      tmp <- format(tmp)}}
  dimnames(r) <- list(tmp,dimnames(x)[[2]])
  storage.mode(r) <- "integer"
# find s <- j's for monotone algorithm
  sj <- .Fortran("sjn",p,npatt,r,integer(p),PACKAGE="cat")[[4]]  
# return list
  list(x=x,n=n,p=p,d=d,jmp=jmp,r=r,nmis=nmis,ro=ro,mdpgrp=mdpgrp,
    mdpgst=mdpgst,mobs=mobs,mobsst=mobsst,nmobs=nmobs,ncells=ncells,
    ngrp=ngrp,npatt=npatt,grouped=grouped,sj=sj)}
#************************************************************************
# Default starting value is a uniform table, and default prior is
# uniform (yielding the MLE). The elements of the starting value need
# not sum to one. Structural zeroes are indicated by zeros in the
# starting value and NAs in the prior.
em.cat <- function(s,start,prior=1,showits=TRUE,maxits=1000,eps=.0001){
  if(length(prior)==1) prior <- rep(prior,s$ncells)
  w <- !is.na(prior)
  if(missing(start)){
    start <- rep(1,s$ncells)
    start[!w] <- 0}
  else{
    if(any(start[w]==0)){
        warning("Starting value on the boundary")}
    if(any(!w)){
      if(any(start[!w]!=0)){
      stop("Starting value has nonzero elements for structural zeroes")}}}
  prior <- as.double(prior); start <- as.double(start); it <- 0; convgd <- FALSE
  if(showits) cat(paste("Iterations of EM:","\n"))
  while((!convgd)&(it<maxits)){
    old <- start
    start <- .Fortran("estepc",s$ncells,start,start,s$npatt,s$p,s$r,
      s$mdpgrp,s$ngrp,s$mobs,s$nmobs,s$d,s$jmp,as.integer(0),
      integer(s$p),integer(s$p),PACKAGE="cat")
    if(start[[13]]==1)
      stop("Bad parameter value: assigns zero prob. to an observed event")
    else start <- start[[3]]
    start[w] <- (start[w]+prior[w]-1)
    start[w] <- start[w]/sum(start[w])
    if(any(start<0))
      stop("Estimate outside the parameter space. Check prior.")
    start[w][start[w]<(.0000001/sum(w))] <- 0
    it <- it+1; if(showits) cat(paste(format(it),"...",sep=""))
  convgd <- all(abs(old-start)<=(eps*abs(old)))}
  if(showits) cat("\n")
  array(start,s$d)}
#************************************************************************
# Data augmentation. Produces a new draw of theta. 
# The default prior is Dirichlet with all hyperparameters=.5 (Jeffreys).
# Unlike em.cat, you must supply a starting value. The elements of the
# starting value need not sum to one. Structural zeros are indicated
# by zeros in the starting value and NAs in the prior.
da.cat <- function(s,start,prior=.5,steps=1,showits=FALSE){
  if(length(prior)==1) prior <- rep(prior,s$ncells)
  w <- !is.na(prior)
  if(any(start[w]==0)){
    warning("Starting value on the boundary")}
  if(any(!w)){
    if(any(start[!w]!=0)){
    stop("Starting value has nonzero elements for structural zeroes")}}
  start <- as.double(start); prior <- as.double(prior)
  if(showits) cat(paste("Steps of Data Augmentation:","\n"))
  for(i in 1:steps){
    if(showits) cat(paste(format(i),"...",sep=""))
    start <- .Fortran("istepc",s$ncells,start,start,s$npatt,s$p,s$r,
      s$mdpgrp,s$ngrp,s$mobs,s$nmobs,s$d,s$jmp,as.integer(0),s$d,s$d,PACKAGE="cat")
    if(start[[13]]==1)
      stop("Bad parameter value: assigns zero prob. to an observed event")
    else start <- start[[3]]
    start[w] <- start[w]+prior[w]
    start[!w] <- -999
    start <- .Fortran("pstep1c",s$ncells,start,start,as.integer(0),PACKAGE="cat")
    if(start[[4]]==1) stop("Improper posterior - check prior.")
    start <- start[[3]]}
  if(showits)cat("\n")
  array(start,s$d)}
#************************************************************************
# generates a single imputed dataset under theta
##############################################
# First need to define function slice.index  #
# per Peter Dalgaard to Fernando Tussell on  #
# 30 March 2000             EFH 17 June 2003 #
# Recent R versions now have this in base,   #
# so now wrapped in if(!exists(...))         #
#                      EFH 30 September 2003 #
##############################################
# if(!exists("slice.index")){
#  slice.index <- function(x,i){
#    k <- dim(x)[i]; sweep(x,i,1:k,function(x,y)y)
#    }
# }
# Above was old def;  replaced by enhanced   #
# def below (R-1.18.1) EFH 30 September 2003 #
##############################################
# New code for slice.index:                  #
if(!exists("slice.index")){
  function (x, MARGIN) 
  {
      d <- dim(x)
      if (is.null(d)) 
          d <- length(x)
      n <- length(d)
      if ((length(MARGIN) > 1) || (MARGIN < 1) || (MARGIN > n)) 
          stop("incorrect value for MARGIN")
      if (any(d == 0)) 
          return(array(integer(0), d))
      y <- rep(rep(seq(1:d[MARGIN]), prod(d[seq(length = MARGIN - 
          1)]) * rep(1, d[MARGIN])), prod(d[seq(from = MARGIN + 
          1, length = n - MARGIN)]))
      dim(y) <- d
      y
  }
}
# End of change        EFH 30 September 2003 #
##############################################
imp.cat <- function(s,theta){
  theta <- as.double(theta/sum(theta))
  if(!(s$grouped)){
    x <- s$x
    x[is.na(x)] <- as.integer(0)
    x <- .Fortran("impc",s$n,s$p,x,s$ncells,theta,s$npatt,s$r,s$mdpgst,
      s$mdpgrp,s$ngrp,s$mobs,s$nmobs,s$mobsst,s$d,s$jmp,s$d,s$d,PACKAGE="cat")[[3]]
    result <- x[s$ro,]}
  if(s$grouped){
    tab <- .Fortran("istepc",s$ncells,theta,theta,s$npatt,s$p,s$r,
      s$mdpgrp,s$ngrp,s$mobs,s$nmobs,s$d,s$jmp,as.integer(0),s$d,s$d,PACKAGE="cat")[[3]]
    tab <- as.integer(round(tab))
    tmp <- array(integer(1),s$d)
    levmat <- matrix(integer(1),s$ncells,s$p)
    for(j in 1:s$p) levmat[,j] <- slice.index(tmp,j)
    result <- list(x=levmat[tab!=0,],counts=tab[tab!=0])}
  result}
#************************************************************************
# Calculates observed-data log-posterior density at theta. If no prior
# is specified, gives observed-data loglikelihood. Structural zeros
# should be indicated by NAs in the prior and zeros in theta.
logpost.cat <- function(s,theta,prior){
  tmp <- .Fortran("llc",s$ncells,as.double(theta),s$npatt,s$p,s$r,
    s$mdpgrp,s$ngrp,s$mobs,s$nmobs,s$d,s$jmp,s$d,s$d,
    numeric(1),as.integer(0),PACKAGE="cat")
  if(tmp[[15]]==1){
    stop("Bad parameter: has zero probability for a nonzero cell")}
  lp <- tmp[[14]]
  if(!missing(prior)){
    lp <- lp+sum((prior[!is.na(prior)]-1)*log(theta[!is.na(prior)]))}
  lp}
#************************************************************************
# Default starting value is a uniform table, and default prior is
# uniform (yielding the MLE). The elements of the starting value need
# not sum to one. Structural zeroes are indicated by a zeroes in the
# starting value and NAs in the prior.
ecm.cat <- function(s,margins,start,prior=1,showits=TRUE,maxits=1000,eps=.0001){
  margins <- as.integer(margins)
  if(length(prior)==1) prior <- rep(prior,s$ncells)
  w <- !is.na(prior)
  if(missing(start)){
    start <- rep(1,s$ncells)
    start[!w] <- 0}
  else{
    if(any(start[w]==0)){
        warning("Starting value on the boundary")}
    if(any(!w)){
      if(any(start[!w]!=0)){
      stop("Starting value has nonzero elements for structural zeroes")}}}
  prior <- as.double(prior); start <- as.double(start); it <- 0; convgd <- FALSE
  eps1 <- .0000001*sum(s$nmobs)/sum(w)
  if(showits) cat(paste("Iterations of ECM:","\n"))
  d <- integer(s$p)
  while((!convgd)&(it<maxits)){
    old <- start
    start <- .Fortran("estepc",s$ncells,start,start,s$npatt,s$p,s$r,
      s$mdpgrp,s$ngrp,s$mobs,s$nmobs,s$d,s$jmp,as.integer(0),
      integer(s$p),integer(s$p),PACKAGE="cat")
    if(start[[13]]==1)
      stop("Bad parameter value: assigns zero prob. to an observed event")
    else start <- start[[3]]
#   do cm-step
    start[w] <- start[w]+prior[w]-1
    if(any(start[w]<0))
      stop("Estimate outside the parameter space. Check prior.")
    start <- .Fortran("ipf",s$ncells,as.double(start),as.double(old),
      length(margins),as.integer(margins),s$p,s$d,s$jmp,
      d,d,d,eps1,PACKAGE="cat")[[3]]
    it <- it+1; if(showits) cat(paste(format(it),"...",sep=""))
    convgd <- all(abs(old-start)<=(eps*abs(old)))}
  if(showits) cat("\n")
  start[w] <- start[w]/sum(start[w])
  array(start,s$d)}
#************************************************************************
# Initializes random number generator seed. Argument should be a 
# positive integer
rngseed <- function(seed){
  seed <- as.integer(seed)
  if(seed<=0)stop("Seed must be a positive integer")
  tmp <- .Fortran("rngs",seed,PACKAGE="cat")
  invisible()}
#************************************************************************
# Monotone data augmentation. Produces a new draw of theta. At each 
# iteration, the missing data Ymis* are filled in to make a 
# monotone pattern. The default prior is Dirichlet with all prior 
# counts = .5. A structural zero is indicated by a starting value
# of zero and a prior hyperparameter of NA.
mda.cat <- function(s,start,steps=1,prior=.5,showits=FALSE){
  start <- as.double(start/sum(start))
  if(length(prior)==1) prior <- rep(prior,s$ncells)
  if((!missing(prior))&(length(prior)==1)){
    prior <- as.double(rep(prior,s$ncells))}
  if((!missing(prior))&any(is.na(prior))){
   if(any(start[is.na(prior)]!=0)){
   stop("Bad starting value: nonzero elements for structural zeroes")}}
  if(showits) cat(paste("Steps of Monotone Data Augmentation:","\n"))
  for(i in 1:steps){
    if(showits) cat(paste(format(i),"...",sep=""))
    start <- .Fortran("mgstepc",s$ncells,start,start,s$npatt,s$p,s$r,
      s$mdpgrp,s$ngrp,s$mobs,s$nmobs,s$d,s$jmp,s$d,s$d,prior,
      prior,as.integer(0),s$sj,PACKAGE="cat")
    if(start[[17]]==1)
      stop("Bad parameter value: assigns zero prob. to an observed event")
    if(start[[17]]==2) stop("Improper posterior distribution")
    start <- start[[3]]
    start <- start/sum(start)}
  if(showits)cat("\n")
  array(start,s$d)}
#************************************************************************
# Iterative proportional fitting in double precision.
# Starting value should have zeroes in positions corresponding to
# structural zeros and uniform values elsewhere. The default
# starting value is a uniform table. Returns a fitted array whose
# sufficient configurations agree with those of table.
ipf <- function(table,margins,start,eps=.0001,maxits=50,showits=TRUE){
  a <- attributes(table)
  d <- dim(table); n <- length(table)
  j <- as.integer(cumprod(d)/d)
  if(!missing(start)){if(length(start)!=n)
    stop("Starting value has incorrect length")}
  if(missing(start)) start <- rep(1,n)
  convgd <- FALSE; it <- 0
  if(showits) cat(paste("Cycles of IPF:","\n",sep=""))
  eps1 <- .0000001*mean(table)
  while((!convgd)&(it<maxits)){
    old <- start
    it <- it+1
    if(showits) cat(paste(format(it),"...",sep=""))
    start <- .Fortran("ipf",n,as.double(table),as.double(old),
      length(margins),as.integer(margins),length(d),d,j,
      d,d,d,eps1,PACKAGE="cat")[[3]]
    convgd <- all(abs(old-start)<=(eps*abs(old)))}
  if(showits) cat("\n")
  attributes(start) <- a
  start}
#************************************************************************
# Data augmentation/Bayesian IPF algorithm. Produces a new draw of
# theta. At each step, the missing data are imputed and one cycle
# of Bayesian IPF is performed.  The default prior is the Jeffreys
# Dirichlet prior (all hyperparameters = .5). Structural zeros
# are indicated by zeros in start and NAs in the prior.
dabipf <- function(s,margins,start,steps=1,prior=.5,showits=FALSE){
  if(length(prior)==1) prior <- rep(prior,s$ncells)
  w <- !is.na(prior)
  if(any(start[w]==0)){
    stop("Starting value on the boundary")}
  if(any(!w)){
    if(any(start[!w]!=0)){
    stop("Starting value has nonzero elements for structural zeroes")}}
  start <- as.double(start); prior <- as.double(prior)
  prior[!w] <- -999
  if(showits)
    cat(paste("Cycles of Data Augmentation/Bayesian IPF:","\n"))
  for(i in 1:steps){
    if(showits) cat(paste(format(i),"...",sep=""))
    old <- start
    start <- .Fortran("istepc",s$ncells,start,start,s$npatt,s$p,s$r,
      s$mdpgrp,s$ngrp,s$mobs,s$nmobs,s$d,s$jmp,as.integer(0),s$d,s$d,PACKAGE="cat")
    if(start[[13]]==1)
      stop("Bad parameter value: assigns zero prob. to an observed event")
    else start <- start[[3]]
    start <- .Fortran("bipf",s$ncells,start,old,prior,
        length(margins),as.integer(margins),s$p,s$d,
        s$jmp,s$d,s$d,s$d,as.integer(0),PACKAGE="cat")
    if(start[[13]]==1) stop("Improper posterior - Check prior.")
    start <- start[[3]]}
  if(showits)cat("\n")
  start <- start/sum(start)  
  array(start,s$d)}
#************************************************************************
# Bayesian IPF algorithm. The default prior is the Jeffreys
# Dirichlet prior (all hyperparameters = .5). Structural zeros
# are indicated by zeros in start and NAs in the prior. 
# The starting value must lie in the interior of the parameter space
# for the loglinear model; the default is an array with zeros
# corresponding to the structural zeros and uniform values elsewhere.
bipf <- function(table,margins,prior=.5,start,steps=1,showits=FALSE){
  a <- attributes(table)
  d <- dim(table); n <- length(table)
  j <- as.integer(cumprod(d)/d)
  if(length(prior)==1)
##FT:    prior <- rep(prior,s$ncells)
##FT:    changed below to,              20.7.2003
    prior <- rep(prior,n)
##FT:    end of changes
  w <- !is.na(prior)
  if(!missing(start)){
     if(length(start)!=n)
         stop("Starting value has incorrect length")
     if(any(start[w]==0))
         stop("Starting value on the boundary")
     if(any(!w)){
         if(any(start[!w]!=0)){
     stop("Starting value has nonzero elements for structural zeroes")}}}
  if(missing(start)){
     start <- rep(1,n)
     start[!w] <- 0
     start[w] <- start[w]/sum(start[w])}
  prior[!w] <- -999
  if(showits) cat(paste("Cycles of Bayesian IPF:","\n",sep=""))
  for(it in 1:steps){
    if(showits) cat(paste(format(it),"...",sep=""))
    tmp <- .Fortran("bipf",n,as.double(table),as.double(start),
      as.double(prior),length(margins),as.integer(margins),
      length(d),d,j,d,d,d,as.integer(0),PACKAGE="cat")
    if(tmp[[13]]==1) stop("Improper posterior - Check prior.")
    start <- tmp[[3]]}
  if(showits) cat("\n")
  attributes(start) <- a
  start}
#************************************************************************
# multiple imputation inference
mi.inference <- function(est,std.err,confidence=.95){
  qstar <- est[[1]]
  for(i in 2:length(est)){qstar <- cbind(qstar,est[[i]])}
  qbar <- apply(qstar,1,mean)
  u <- std.err[[1]]
  for(i in 2:length(std.err)){u <- cbind(u,std.err[[i]])}
########################################################################
##FT: 24.7.2003  Added next line, otherwise aborts when invoked with scalar
##               estimands
  if (!is.null(dimnames(qstar)[[1]]))
##  End of changes #####################################################
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
