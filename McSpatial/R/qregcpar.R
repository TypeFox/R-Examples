qregcpar <- function(form,nonpar,taumat=c(.10,.25,.50,.75,.90),window=.25,bandwidth=0,kern="tcub",
  distance="Mahal",target=NULL,data=NULL) {

  mat <- model.frame(form,data=data)
  y <- mat[,1]
  n = length(y)
  xmat <- model.matrix(form,data=data)
  zmat <- model.frame(nonpar,data=data)
  nz = ncol(zmat)
  nk = ncol(xmat)
  ntau = length(taumat)
  xname <- colnames(xmat)

  if (nz==1) {vzmat <- var(zmat) }
  if (nz==2) {
    vzmat <- cov(zmat) 
    if (distance=="Euclid"|distance=="E") {vzmat <- diag(diag(vzmat)) }
  }

  if (kern=="rect")  { wgt <- function(psi) {1 } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { 1-psi^2 } }
  if (kern=="bisq")  { wgt <- function(psi) { (1-psi^2)^2 } }
  if (kern=="tcub")  { wgt <- function(psi) { (1 - abs(psi)^3)^3 } }
  if (kern=="trwt")  { wgt <- function(psi) { (1 - psi^2)^3 } }
  if (kern=="gauss") { wgt <- function(psi) { exp(-((2.5*psi)^2)/2) } }

  if (bandwidth>0) {window = 0}

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

  xcoef.target <-    array(0,dim=c(nt,ntau,nk))
  xcoef.target.se <- array(0,dim=c(nt,ntau,nk))
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
    k <- wgt(dist[samp]/h)

    for (j in seq(1:ntau)) {
      fit <- summary(rq(y[samp]~xmat[samp,-1], weights=k, tau=taumat[j]), covariance=TRUE)
      xcoef.target[i,j,] <- fit$coef[1:nk]
      xcoef.target.se[i,j,] <- diag(fit$cov)
    }
  }

  if (alldata==TRUE) {
    xcoef <- xcoef.target
    xcoef.se <- xcoef.target.se
  }

  if (alldata==FALSE) {
    xcoef <-    array(0,dim=c(n,ntau,nk))
    xcoef.se <- array(0,dim=c(n,ntau,nk))
    for (itau in seq(1:ntau)) {
    for (j in seq(1:nk)) {
      xcoef[,itau,j] <- smooth12(target,xcoef.target[,itau,j],zmat)
      xcoef.se[,itau,j] <- smooth12(target,xcoef.target.se[,itau,j],zmat)
    }
    }
  }

  yhat <- array(0,dim=c(n,ntau))
  for (itau in seq(1:ntau)) {
    yhat[,itau] <- rowSums(xmat*xcoef[,itau,])
  }
    
  out <- list(target,xcoef.target,xcoef.target.se,xcoef,xcoef.se,yhat)
  names(out) <- c("target","xcoef.target","xcoef.target.se","xcoef","xcoef.se","yhat")
  return(out)   
  detach(data)
}


