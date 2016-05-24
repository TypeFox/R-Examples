lwr <- function(form,window=.25,bandwidth=0,kern="tcub",distance="Mahal",target=NULL,data=NULL) {

  form <- as.formula(form,env=data)

  mat <- model.frame(form,data=data)
  y <- mat[,1]
  n = length(y)
  xmat <- model.matrix(form,data=data)
  xmat <- as.matrix(xmat[,-1])
  nk = ncol(xmat)
  if (nk==1) {vxmat <- var(xmat) }
  if (nk==2) {
    vxmat <- cov(xmat) 
    if (distance=="Euclid"|distance=="E") {vxmat <- diag(diag(vxmat)) }
  }

  if (kern=="rect")  { wgt <- function(psi) {ifelse(abs(psi)>=0,1,0) } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { 1-psi^2 } }
  if (kern=="bisq")  { wgt <- function(psi) { (1-psi^2)^2 } }
  if (kern=="tcub")  { wgt <- function(psi) { (1 - abs(psi)^3)^3 } }
  if (kern=="trwt")  { wgt <- function(psi) { (1 - psi^2)^3 } }
  if (kern=="gauss") { wgt <- function(psi) { exp(-((2.5*psi)^2)/2) } }

  if (identical(target,NULL)){
    target <- maketarget(form,window=window,bandwidth=bandwidth,kern="tcub",data=data)$target
  }
  alldata = FALSE
  if (identical(target,"alldata")){
    target <- xmat
    alldata = TRUE
  }
  if (bandwidth>0){window = 0}
  target <- as.matrix(target)
  nt = nrow(target)
  
  ytarget     <- array(0,dim=nt)
  ytarget.se  <- array(0,dim=nt)
  dtarget1    <- array(0,dim=nt)
  dtarget1.se <- array(0,dim=nt)
  dtarget2    <- array(0,dim=nt)
  dtarget2.se <- array(0,dim=nt)
  df1target   <- array(0,dim=nt)
  df2target   <- array(0,dim=nt)

  if (distance=="Latlong"|distance=="L") {
    tvect <- attr(terms(form),"term.labels")
    if (substr(tvect[1],1,2)=="la"|substr(tvect[1],1,2)=="La"|substr(tvect[1],1,2)=="LA") {
      la  <- 2*pi*xmat[,1]/360 
      la1 <- 2*pi*target[,1]/360
    }
    if (substr(tvect[2],1,2)=="la"|substr(tvect[2],1,2)=="La"|substr(tvect[2],1,2)=="LA") {
      la  <- 2*pi*xmat[,2]/360 
      la1 <- 2*pi*target[,2]/360 
    }
    if (substr(tvect[1],1,2)=="lo"|substr(tvect[1],1,2)=="Lo"|substr(tvect[1],1,2)=="LO") {
      lo  <- 2*pi*xmat[,1]/360 
      lo1 <- 2*pi*target[,1]/360
    }
    if (substr(tvect[2],1,2)=="lo"|substr(tvect[2],1,2)=="Lo"|substr(tvect[2],1,2)=="LO") {
      lo  <- 2*pi*xmat[,2]/360 
      lo1 <- 2*pi*target[,2]/360
    }
  }

  for (i in seq(1:nt)) {
    if (distance!="Latlong"&distance!="L")  dist <- sqrt(mahalanobis(xmat, target[i,], vxmat))
    if (distance=="Latlong"|distance=="L") {
      dist <- pmin(sin(la)*sin(la1[i]) + cos(la)*cos(la1[i])*cos(lo1[i]-lo),  1)
      dist <- acos(dist)*3958
    }
    if (window>0) {h = quantile(dist,window) }
    if (bandwidth>0) {h = bandwidth}
    samp <- dist<=h
    if (kern=="gauss") {samp <- dist<=max(dist)}
    
    xmat1 <- cbind(1,xmat[samp,1]-target[i,1])
    if (nk==2) {xmat1 <- cbind(xmat1,(xmat[samp,2]-target[i,2])) }
    k <- wgt(dist[samp]/h)
    xmat2 <- k*xmat1
    xx <- solve(crossprod(xmat1,xmat2))
    xmat1 <- xx%*%t(xmat2)
    bmat <- xmat1%*%y[samp]
    ytarget[i] = bmat[1,1]
    dtarget1[i] = bmat[2,1]
    if (nk==2) {dtarget2[i] = bmat[3,1] }
    df1target[i] = xx[1,1]
    df2target[i] = sum(xmat1[1,]^2)
    vmat <- tcrossprod(xmat1)
    ytarget.se[i] = sqrt(vmat[1,1])
    dtarget1.se[i] = sqrt(vmat[2,2])
    if (nk==2) {dtarget2.se[i] = sqrt(vmat[3,3]) }
  }

  if (alldata==FALSE) {
    yhat <- smooth12(target,ytarget,xmat)
    dhat1 <- smooth12(target,dtarget1,xmat)
    infl <- smooth12(target,df1target,xmat)
    df1 <- sum(infl)
    df2 <- sum(smooth12(target,df2target,xmat))
    yhat.se <- smooth12(target,ytarget.se,xmat)
    dhat1.se <- smooth12(target,dtarget1.se,xmat)
    dhat2 <- NULL
    dhat2.se <- NULL
    if (nk==2){
      dhat2 <- smooth12(target,dtarget2,xmat)
      dhat2.se <- smooth12(target,dtarget2.se,xmat)
    }
  }
    
  if (alldata==TRUE) {
    yhat <- ytarget
    dhat1 <- dtarget1
    dhat2 <- dtarget2
    yhat.se  <- ytarget.se
    dhat1.se <- dtarget1.se
    dhat2.se <- dtarget2.se
    infl <- df1target
    df1 = sum(infl)
    df2 = sum(df2target)
  }

  rss = sum((y-yhat)^2)
  sig2 = rss/(n-2*df1 + df2)
  cv = mean(((y-yhat)/(1-infl))^2)
  gcv = n*rss/((n-df1)^2)
  ytarget.se  <- sqrt(sig2)*ytarget.se
  dtarget1.se <- sqrt(sig2)*dtarget1.se
  dtarget2.se <- sqrt(sig2)*dtarget2.se
  yhat.se  <- sqrt(sig2)*yhat.se
  dhat1.se <- sqrt(sig2)*dhat1.se
  dhat2.se <- sqrt(sig2)*dhat2.se

  out <- list(target,ytarget,dtarget1,dtarget2,ytarget.se,dtarget1.se,dtarget2.se,
                    yhat,dhat1,dhat2,yhat.se,dhat1.se,dhat2.se,df1,df2,sig2,cv,gcv,infl)
  names(out) <- c("target","ytarget","dtarget1","dtarget2","ytarget.se","dtarget1.se","dtarget2.se",
                  "yhat","dhat1","dhat2","yhat.se","dhat1.se","dhat2.se","df1","df2","sig2","cv","gcv","infl")
  return(out)    
}

