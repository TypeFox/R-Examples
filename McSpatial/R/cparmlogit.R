cparmlogit <- function(form,nonpar,window=.25,bandwidth=0,kern="tcub",distance="Mahal",target=NULL,data=NULL) {

  xmat <- model.frame(form,data=data)
  lvect <- levels(factor(xmat[,1]))
  nchoice = length(lvect)
  K = nchoice-1
  n = nrow(xmat)
  nk = ncol(xmat)
  y <- array(0,dim=n)
  for (j in seq(2,nchoice)) {
    y <- ifelse(factor(xmat[,1])==lvect[j],j-1,y)
  }
  xmat[,1] <- 1
  zmat <- model.frame(nonpar,data=data)
  nz = ncol(zmat)
  y <- as.numeric(y)

  if (nz==1) {vzmat <- var(zmat) }
  if (nz==2) {
    vzmat <- cov(zmat) 
    if (distance=="Euclid"|distance=="E") {vzmat <- diag(diag(vzmat)) }
  }

  if (kern=="rect")  { wgt <- function(psi) {ifelse(abs(psi)>=0,1,0) } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { 1-psi^2 } }
  if (kern=="bisq")  { wgt <- function(psi) { (1-psi^2)^2 } }
  if (kern=="tcub")  { wgt <- function(psi) { (1 - abs(psi)^3)^3 } }
  if (kern=="trwt")  { wgt <- function(psi) { (1 - psi^2)^3 } }
  if (kern=="gauss") { wgt <- function(psi) { exp(-((2.5*psi)^2)/2) } }

  if (bandwidth>0) {window = 0}

  if (nz==1&window>0)    {fit <- locfit(~lp(zmat[,1],nn=window,deg=1),kern=kern) }
  if (nz==2&window>0)    {fit <- locfit(~lp(zmat[,1],zmat[,2],nn=window,deg=1),kern=kern) }
  if (nz==1&bandwidth>0) {fit <- locfit(~lp(zmat[,1],h=2*bandwidth,deg=1),kern=kern) }
  if (nz==2&bandwidth>0) {fit <- locfit(~lp(zmat[,1],zmat[,2],h=2*bandwidth,deg=1),kern=kern) }
 
  if (identical(target,NULL)){
    target <- maketarget(nonpar,window=window,bandwidth=bandwidth,kern="tcub",data=data)$target
  }
  alldata = FALSE
  if (identical(target,"alldata")){
    target <- xmat
    alldata = TRUE
  }
  if (bandwidth>0){window = 0}
  target <- as.matrix(target)
  nt = nrow(target)


  if (distance=="Latlong"|distance=="L") {
    tvect <- attr(terms(nonpar),"term.labels")
    if (substr(tvect[1],1,2)=="la"|substr(tvect[1],1,2)=="La"|substr(tvect[1],1,2)=="LA") {
      la  <- 2*pi*zmat[,1]/360 
      la1 <- 2*pi*target[,1]/360
    }
    if (substr(tvect[2],1,2)=="la"|substr(tvect[2],1,2)=="La"|substr(tvect[2],1,2)=="LA") {
      la  <- 2*pi*zmat[,2]/360 
      la1 <- 2*pi*target[,2]/360 
    }
    if (substr(tvect[1],1,2)=="lo"|substr(tvect[1],1,2)=="Lo"|substr(tvect[1],1,2)=="LO") {
      lo  <- 2*pi*zmat[,1]/360 
      lo1 <- 2*pi*target[,1]/360
    }
    if (substr(tvect[2],1,2)=="lo"|substr(tvect[2],1,2)=="Lo"|substr(tvect[2],1,2)=="LO") {
      lo  <- 2*pi*zmat[,2]/360 
      lo1 <- 2*pi*target[,2]/360
    }
  }

  bstart <- array(0,dim=c(nk,K)) 
  for (k in seq(1:K)) {
    p0 <- y==k
    fit <- glm(p0[y==0|y==k]~as.matrix(xmat)[y==0|y==k,-1],family=binomial(link="logit"))
    bstart[,k] <- fit$coef
  }
  bstart <- c(bstart)

  pmat <- array(0,dim=c(n,K))

  logl <- function(bmat) {
    b <- array(bmat,dim=c(nk,K))
    pmat <- array(0,dim=c(n1,K))

    for (j in seq(1:K)) {
      pmat[,j] <- exp(as.matrix(xmat1)%*%b[,j])
    }
    p0 <- 1/(1+rowSums(pmat))
    pmat <- p0*as.data.frame(pmat)

    gmat <- array(0,dim=c(nk,K))
    for (j in seq(1,K)) {
      u <- ifelse(y[samp]==j,1,0)-pmat[,j]
      gmat[,j] <- t(as.matrix(xmat1))%*%as.array(w*u)
    }
    gmat <- c(gmat)

    pmat <- cbind(p0, pmat)
    pmat <- w*(log(pmat))
    for (j in seq(0,K)) {
      pmat[y[samp]!=j,j+1] <- 0
    }

    out <- -sum(pmat)
    attr(out,"gradient") <- -gmat
    return(out)
    
  }

  xcoef.target <- array(0,dim=c(nt,nk*K))
  xcoef.target.se <- array(0,dim=c(nt,nk*K))



  for (i in seq(1:nt)) {
    if (distance!="Latlong"&distance!="L") {dist <- sqrt(mahalanobis(zmat, target[i,], vzmat)) }
    if (distance=="Latlong"|distance=="L") {
      dist <- pmin(sin(la)*sin(la1[i]) + cos(la)*cos(la1[i])*cos(lo1[i]-lo),  1)
      dist <- acos(dist)*3958
    }
    if (window>0) {h = quantile(dist,window) }
    if (bandwidth>0) {h = bandwidth}
    samp <- dist<=h
    if (kern=="gauss") {samp <- dist<=max(dist)}

    xmat1 <- as.matrix(xmat[samp,]) 
    w <- wgt(dist[samp]/h)
    n1 = length(w)

    nlfit <- nlm(logl,bstart,iterlim=10000)
    xcoef.target[i,] <- nlfit$estimate
    b <- array(nlfit$estimate,dim=c(nk,K))

    pmat <- array(0,dim=c(n1,K))
    for (j in seq(1:K)) {
       pmat[,j] <- exp(xmat1%*%b[,j])
    }
    p0 <- 1/(1+rowSums(pmat))
    pmat <- p0*as.data.frame(pmat)
    umat <- array(0,dim=c(n1,K))
    for (j in seq(1:K)) {
      umat[,j] <- ifelse(y[samp]==j,1,0) - pmat[,j]
    }
    vmat1 <- array(0,dim=c(nk*K,nk*K)) 
    vmat <- array(0,dim=c(nk*k,nk*K))
    for (ii in seq(1,K)) {
    for (j in seq(1,ii)) {
      if (ii==j) {
        vmat1[(1+(ii-1)*nk):(ii*nk),(1+(j-1)*nk):(j*nk)] <- crossprod(as.matrix((sqrt(w*pmat[,j]*(1-pmat[,j]))*xmat1)))
         vmat[(1+(ii-1)*nk):(ii*nk),(1+(j-1)*nk):(j*nk)] <- crossprod(as.matrix((w*umat[,ii])*xmat1))
      }
      if (ii!=j) {
        vmat1[(1+(ii-1)*nk):(ii*nk),(1+(j -1)*nk):( j*nk)] <- -crossprod(as.matrix(sqrt(w*pmat[,ii]*pmat[,j])*xmat1))
        vmat1[(1+ (j-1)*nk):( j*nk),(1+(ii-1)*nk):(ii*nk)] <- vmat1[(1+(ii-1)*nk):(ii*nk),(1+(j-1)*nk):(j*nk)]
          vmat[(1+(ii-1)*nk):(ii*nk),(1+(j-1)*nk):(j*nk)] <- crossprod(as.matrix(w*umat[,ii]*xmat1),as.matrix(w*umat[,j]*xmat1))
          vmat[(1+(j-1)*nk):(j*nk),(1+(ii-1)*nk):(ii*nk)] <- vmat[(1+(ii-1)*nk):(ii*nk),(1+(j-1)*nk):(j*nk)]
      }
    }
    }

  vmat1 <- solve(vmat1)
  vmat <- vmat1%*%vmat%*%vmat1
  vmat <- sqrt(diag(vmat))
 
  xcoef.target.se[i,] <- vmat

 }

  xcoef <- array(0,dim=c(n,nk*K))
  xcoef.se <- array(0,dim=c(n,nk*K))
  for (j in seq(1:ncol(xcoef))) {
    xcoef[,j] <- smooth12(target,xcoef.target[,j],zmat)
    xcoef.se[,j] <- smooth12(target,xcoef.target.se[,j],zmat)
  }

  pmat <- array(0,dim=c(n,K))
  for (j in seq(1:K)) {
     b <- xcoef[,(1+(j-1)*nk):(j*nk)]
     pmat[,j] <- exp(rowSums(xmat*b))
  }
  p0 <- 1/(1+rowSums(pmat))
  pmat <- p0*as.data.frame(pmat)
  pmat <- cbind(p0,pmat)

  p0 <- log(p0)
  for (j in seq(1:K)) {
    p0 <- ifelse(y==j,log(pmat[,j+1]),p0)
  }
  lnl = sum(p0)

  rname <- colnames(xmat)
  rname[1] <- "Intercept"
  colnames(xcoef.target) <- rep(rname,K)
  colnames(xcoef.target.se) <- rep(rname,K)
  colnames(xcoef) <- rep(rname,K)
  colnames(xcoef.se) <- rep(rname,K)
  colnames(pmat) <- paste("p",seq(0,K),sep="")

  out <- list(target,xcoef.target, xcoef.target.se,xcoef,xcoef.se, pmat, lnl)
  names(out) <- c("target","xcoef.target", "xcoef.target.se","xcoef","xcoef.se","pmat","lnl")
  return(out)

}






