# cross-validation
robustgam.crossval <- function(X, y, family, p=3, K=30, c=1.345, sp, show.msg=FALSE, count.lim=200, w.count.lim=50, smooth.basis="tp", wx=FALSE, ngroup=length(y), seed=12345){
  if (family$family=="poisson") expect <- expect.poisson
  else if (family$family=="binomial") expect <- expect.binomial
  else if (family$family=="gaussian") expect <- expect.gaussian

  if (family$family=="poisson"){
     w.fun1 <- w.poisson
  } else if (family$family=="binomial"){
     w.fun1 <- w.binomial
  }

  if (wx){
    w.fun <- function(m.initial,c,X){
      w <- w.fun1(m.initial, c)
      mcd.est <- covMcd(X)
      return(w*sqrt(1/(1+mcd.est$mah)))
    }
  } else {
    w.fun <- function(m.initial,c,X){
      return(w.fun1(m.initial, c))
    }
  }

  #dividing the groups
  set.seed(seed)
  n <- length(y)
  leave.out <- trunc(n/ngroup)
  o <- sample(1:n)
  groups <- vector("list", ngroup)
  for (j in (1:(ngroup-1))){
    jj <- (1+(j-1)*leave.out)
    groups[[j]]<-(o[jj:(jj+leave.out-1)])
  }
  groups[[ngroup]] <- o[(1+(ngroup-1)*leave.out):n]

  #
  u <- vector("list", ngroup)
  cv.fit <- rep(NA, n)
  X <- matrix(X,nrow=length(y))
  data <- data.frame(data.frame(X),data.frame(y))
  for(j in (1:ngroup)){
    u[[j]] <- robustgam(X=X[-groups[[j]],], y=y[-groups[[j]]], family=family, p=p, K=K, c=c, sp=sp, show.msg=show.msg, count.lim=count.lim, w.count.lim=w.count.lim,smooth.basis=smooth.basis, wx=wx)
    cv.fit[groups[[j]]] <- pred.robustgam(u[[j]],data[groups[[j]],],type="response")$predict.values
  }

  #
  int1 <- array(dim=n)
  int2 <- array(dim=n)
  if (family$family=="poisson") {
    cv.fit <- 1e-6*(cv.fit<(1e-6)) + cv.fit*(cv.fit>=(1e-6))
  }
  if (family$family=="binomial") {
    cv.fit <- 1e-6*(cv.fit<(1e-6)) + (1-1e-6)*(cv.fit>(1-1e-6)) + cv.fit*(cv.fit>=(1e-6))*(cv.fit<=(1-1e-6))
  }

  w <- w.fun(cv.fit, c, X)

  for (i in (1:n)){
    int1[i] <- integrate(function(t, c1, y1){Huber.deriv((y1-t)/sqrt(family$variance(t)),c1)/sqrt(family$variance(t))}, lower=y[i], upper=cv.fit[i],c1=c, y1=y[i], subdivisions=500, stop.on.error=FALSE)$value
    int2[i] <- integrate(function(t, c1){expect(t,c1,sqrt(family$variance(t)))/sqrt(family$variance(t))}, lower=y[i], upper=cv.fit[i],c1=c, subdivisions=500, stop.on.error=F)$value
  }
  Qm <- -(sum(int1*w)-sum(int2*w))

  return(list(Qm=Qm, cv.fit=cv.fit, sp=sp, ngroup=ngroup, groups=groups, u=u))
}

robustgam.CV <- function(X, y, family, p=3, K=30, c=1.345, show.msg=FALSE, count.lim=200, w.count.lim=50, smooth.basis="tp", wx=FALSE, sp.min=1e-7, sp.max=1e-3, len=50, show.msg.2=TRUE, ngroup=length(y), seed=12345){
  X <- matrix(X,nrow=length(y))
  nxs <- ncol(X)
  if (length(sp.min)==1){
    sp.min <- rep(sp.min,nxs)
  }
  if (length(sp.max)==1){
    sp.max <- rep(sp.max,nxs)
  }
  if (length(len)==1){
    len <- rep(len,nxs)
  }

  sp <- list()
  for (j in (1:nxs)){
    sp[[j]] <- exp(seq(log(sp.min[j]), log(sp.max[j]), len=len[j]))
  }

  criteria <- array(dim=len)
  for (i in (1:prod(len))){
    tempind <- arrayInd(i,len)
    tempsp <- sapply(1:nxs,function(jj,sp,ind){sp[[jj]][ind[jj]]},sp=sp,ind=as.vector(tempind))
    temp <- robustgam.crossval(X=X, y=y, family=family, p=p, K=K, c=c, sp=tempsp, show.msg=show.msg, count.lim=count.lim, w.count.lim=w.count.lim, smooth.basis=smooth.basis, wx=wx, ngroup=ngroup, seed=seed)
    criteria[tempind] <- temp$Qm
    if (show.msg.2) {cat("#### ", as.vector(tempind), ": sp=",tempsp," finished. ####\n")}
  }
  optim.index <- which.min(criteria)
  optim.index2 <- arrayInd(optim.index,len)
  optim.sp <- sapply(1:nxs,function(jj,sp,ind){sp[[jj]][ind[jj]]},sp=sp,ind=as.vector(optim.index2))
  optim.fit <- robustgam(X=X, y=y, family=family, p=p, K=K, c=c, sp=optim.sp, show.msg=show.msg, count.lim=count.lim, w.count.lim=w.count.lim, smooth.basis=smooth.basis, wx=wx)
  return(list(fitted.values=optim.fit$fitted.values, initial.fitted=optim.fit$initial.fitted, beta=optim.fit$beta, optim.index=optim.index, optim.index2=optim.index2, optim.criterion=criteria[optim.index2], optim.sp=optim.sp, criteria=criteria, sp=sp, optim.fit=optim.fit))
}
