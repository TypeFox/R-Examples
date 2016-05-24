cparlwr <- function(form,nonpar,window=.25,bandwidth=0,kern="tcub",distance="Mahal",targetobs=NULL,data=NULL) {

  mat <- model.frame(form,data=data)
  y <- mat[,1]
  n = length(y)
  xmat <- model.matrix(form,data=data)
  zmat <- model.frame(nonpar,data=data)
  nz = ncol(zmat)
  nk = ncol(xmat)

  xname <- colnames(xmat)
  nocons = xname[1]!="(Intercept)"

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

  alldata = FALSE
  ttype = "provided"
  if (identical(targetobs,NULL)){ttype = "locfit"}
  if (identical(targetobs,"alldata")){ttype = "alldata"}

  if (ttype=="provided"){
    target <- zmat[targetobs,]
    xvect <- xmat[targetobs,]
  }
  if (ttype=="locfit"){
    fit <- maketarget(nonpar, window=window, bandwidth=bandwidth, kern=kern, actualobs=TRUE, data=data)
    target <- fit$target
    xvect <- xmat[fit$obs,]
    obs <- fit$obs
  }
  if (ttype=="alldata"){
    target <- zmat
    xvect <- xmat
    obs <- seq(1,n)
  }
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

  xcoef.target <- array(0,dim=c(nt,nk))
  df1target <- array(0,dim=nt)
  df2target <- array(0,dim=nt)
  xcoef.target.se <- array(0,dim=c(nt,nk))

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

    k <- wgt(dist[samp]/h)
    xmat2 <- k*xmat1
    xx <- solve(crossprod(xmat1,xmat2))
    xmat1 <- xx%*%t(xmat2)
    xcoef.target[i,] <- xmat1%*%y[samp]
    vmat <- tcrossprod(xmat1)

    xmat1 <- xvect[i,]%*%xmat1
    df2target[i] = sum(xmat1^2)
    dist[samp] <- xmat1
    dist[!samp] <- 0
    df1target[i] = dist[obs[i]]

    xcoef.target.se[i,] <- sqrt(diag(vmat))
  }

  xcoef <- array(0,dim=c(n,nk))
  xcoef.se <- array(0,dim=c(n,nk))

  if (alldata==FALSE) {
    for (j in seq(1:nk)) {
      xcoef[,j] <- smooth12(target,xcoef.target[,j],zmat)
      xcoef.se[,j] <- smooth12(target,xcoef.target.se[,j],zmat)
     }
    infl <- smooth12(target,df1target,zmat)
    df1 = sum(infl)
    df2 <- sum(smooth12(target,df2target,zmat))
  }

  if (alldata==TRUE) {
    xcoef <- xcoef.target
    xcoef.se <- xcoef.target.se
    infl <- df1target
    df1 = sum(infl)
    df2 = sum(df2target)
  }

#  yhat <- diag(tcrossprod(xmat,xcoef))   
  yhat <- xmat[,1]*xcoef[,1]
  if (nk>1) {
    for (j in seq(2,nk)) {
      yhat <- yhat + xmat[,j]*xcoef[,j]
    }
  }

  rss = sum((y-yhat)^2)
  sig2 = rss/(n-2*df1 + df2)
  cv = mean(((y-yhat)/(1-infl))^2)
  gcv = n*rss/((n-df1)^2)

  xcoef.target.se <- sqrt(sig2)*xcoef.target.se
  xcoef.se <- sqrt(sig2)*xcoef.se

  colnames(xcoef.target) <- xname
  colnames(xcoef.target.se) <- xname
  colnames(xcoef) <- xname
  colnames(xcoef.se) <- xname
  ytarget <- yhat[obs]

  out <-      list(target,  ytarget,  xcoef.target,  xcoef.target.se,  yhat,  xcoef,  xcoef.se,  df1,  df2,  sig2,  cv,  gcv,  infl)
  names(out) <- c("target","ytarget","xcoef.target","xcoef.target.se","yhat","xcoef","xcoef.se","df1","df2","sig2","cv","gcv","infl")
  return(out)   
  detach(data)
}


