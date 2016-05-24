w.poisson <- function(m, c){
    k1 <- m-c*sqrt(m)
    k2 <- m+c*sqrt(m)
    1/(ppois(k2,m)-ppois(ceiling(k1-1),m))
}

w.binomial <- function(m, c){
    k1 <- m-c*sqrt(m*(1-m))
    k2 <- m+c*sqrt(m*(1-m))
    1/(pbinom(k2,1,m)-pbinom(ceiling(k1-1),1,m))
}

robustgam <- function(X, y, family, p=3, K=30, c=1.345, sp=-1, show.msg=FALSE, count.lim=200, w.count.lim=50, smooth.basis="tp", wx=FALSE){
  X <- matrix(X,nrow=length(y))
  nxs <- ncol(X)
  if (length(sp)==1){
    sp <- rep(sp,nxs)
  }
  if ((length(sp)!=nxs)||(!prod(sp>=0))){
    stop('Please specify smoothing parameter!\n')
  }

  n <- length(y)
  data <- data.frame(data.frame(X),data.frame(y))
  # select splines
  basis <- list()
  for (j in (1:nxs)){
    if (smooth.basis=="tr"){
      stemp <- s(X[,j], bs="tr", m=p, k=p+K+1)
      stemp$term <- names(data)[j]; stemp$label <- paste(c("s(",names(data)[j],")"),collapse="")
      basis[[j]]  <- smooth.construct.tr.smooth.spec(stemp,data,NULL)
    } else if (smooth.basis=="tp") {
      stemp <- s(X[,j], bs="tp", m=p)
      stemp$term <- names(data)[j]; stemp$label <- paste(c("s(",names(data)[j],")"),collapse="")
      basis[[j]]  <- smooth.construct.tp.smooth.spec(stemp,data,NULL)
    } else if (smooth.basis=="cr") {
      stemp <- s(X[,j], bs="cr", k=K)
      stemp$term <- names(data)[j]; stemp$label <- paste(c("s(",names(data)[j],")"),collapse="")
      basis[[j]]  <- smooth.construct.cr.smooth.spec(stemp,data,NULL)
    } else if (smooth.basis=="ps") {
      stemp <- s(X[,j], bs="ps", m=c(2,2), k=2+K+2)
      stemp$term <- names(data)[j]; stemp$label <- paste(c("s(",names(data)[j],")"),collapse="")
      basis[[j]]  <- smooth.construct.ps.smooth.spec(stemp,data,NULL)
      basis[[j]]$df <- dim(basis[[j]]$X)[2]
    }
  }

  # impose the identifiablity constraint
  dfs <- sapply(basis,function(b){b$df})
  B <- basis[[1]]$X
  sD <- matrix(0, nrow=dfs[1]+sum(dfs[-1]-1), ncol=dfs[1]+sum(dfs[-1]-1))
  sD[1:dfs[1],1:dfs[1]] <- sp[1]*basis[[1]]$S[[1]]
  Z <- list()
  Z[[1]] <- NULL
  if (nxs>1){
    tempindex <- dfs[1]
    for (j in (2:nxs)){
      Z[[j]] <- qr.Q(qr(t(basis[[j]]$X)%*%rep(1,n)), TRUE)[,-1]
      B <- cbind(B,basis[[j]]$X%*%Z[[j]])
      sD[(tempindex+1):(tempindex+dfs[j]-1),(tempindex+1):(tempindex+dfs[j]-1)] <- sp[j]*(t(Z[[j]])%*%basis[[j]]$S[[1]]%*%Z[[j]])
      tempindex <- tempindex + dfs[j] - 1
    }
  }
  rS <- mat.sqrt(sD)

  # choose the fisher consistency correction
  if (family$family=="poisson"){
   expect <- expect.poisson
   } else if (family$family=="binomial"){
     expect <- expect.binomial
   } else if (family$family=="gaussian"){
     expect <- expect.gaussian
  }

  # choose weights
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

  # initial estimate
  fit <- fit.gam.sp1(y, B, rS, family)
  beta.old <- as.vector(fit$coefficients)
  m.initial <- fit$fitted.values

  main.fit <- rgam.fast.main(Ry=y, RX=X, Rfamily=family, Rc=c, RB=B, RrS=rS, Rexpect=expect, Rw_fun=w.fun, Rm_initial=m.initial, Rbeta_old=beta.old, Rcount_lim=count.lim, Rw_count_lim=w.count.lim,Rdisplay=show.msg)

  #if (main.fit$converge!=1) {
  #  cat("The algorithm did not converage!\n")
  #}

  # reconstructing the beta in original basis representation
  beta.fit <- main.fit$beta
  beta <- beta.fit[1:dfs[1]]
  if (nxs>1){
    tempindex <- dfs[1]
    for (j in (2:nxs)){
      beta <- c(beta,Z[[j]]%*%beta.fit[(tempindex+1):(tempindex+dfs[j]-1)])
      tempindex <- tempindex + dfs[j]
     }
  }
  return(list(fitted.values=as.vector(main.fit$fitted.values),initial.fitted=m.initial,beta=beta,B=B,sD=sD,basis=basis,converage=(main.fit$converge==1), w=as.vector(main.fit$w),family=family,wx=wx, beta.fit=beta.fit))
}

pred.robustgam <- function(fit, data, type="response"){
  # type only affects fitted
  nxs <- length(fit$basis)
  n <- nrow(data)
  dfs <- sapply(fit$basis,function(b){b$df})
  predict.comp <- matrix(ncol=nxs,nrow=n)
  tempindex <- 0
  for (j in (1:nxs)){
    if (attr(fit$basis[[j]],"class")=="tr.smooth"){
      predict.comp[,j] <- as.vector(as.matrix(Predict.matrix.tr.smooth(fit$basis[[j]], data))%*%fit$beta[(tempindex+1):(tempindex+dfs[j])])
    } else {
      predict.comp[,j] <- as.vector(as.matrix(Predict.matrix(fit$basis[[j]], data))%*%fit$beta[(tempindex+1):(tempindex+dfs[j])])
    }
    tempindex <- tempindex + dfs[j]
  }
  fitted <- apply(predict.comp,1,sum)
  if (type=="link") {
    out <- fitted
  } else if (type=="response") {
    out <- fit$family$linkinv(fitted)
  } else {
    stop("no such option for type\n")
  }
  return(list(predict.comp=predict.comp, predict.values=out))
}

cplot.robustgam <- function(fit, ranges, len=100){
  if ("optim.fit"%in%names(fit)){
    fit <- fit$optim.fit
  }
  ranges <- matrix(ranges, nrow=2)
  nxs <- length(fit$basis)
  X <- matrix(nrow=len, ncol=nxs)
  for (j in (1:nxs)){
    X[,j] <- seq(ranges[1,j], ranges[2,j], len=len)
  }
  dat1 <- data.frame(X=X)
  if (nxs>1){
    names(dat1) <- paste("X", 1:nxs, sep="")
  }
  pre.res <- pred.robustgam(fit=fit, data=dat1, type="response")
  pre.comp <- pre.res$predict.comp
  pre.comp[,1] <- pre.comp[,1]/mean(pre.comp[,1])
  for (j in (1:nxs)){
    plot(X[,j], pre.comp[,j], xlab=paste(c("x",j), collapse=""), ylab=paste(c("f",j), collapse=""), type="l", main=paste(c("Component ",j),collapse=""))
    readline(prompt = "Pause. Press <Enter> to continue...")
  }
  return(invisible())
}
