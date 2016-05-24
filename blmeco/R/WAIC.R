WAIC <- function(mod, bsim=NA, nsim=100){

  mclass <- class(mod)[1]
  if(!is.element(mclass, c("lm", "glm", "lmerMod", "glmerMod"))) stop("mod must be an lm, glm or mer object")

  if(is.element(mclass, c("lm", "glm"))) y <- mod$model[,1] else y <- mod@frame[,1]
  n <- length(y)
  lppdj <- numeric(n)
  lEp <- numeric(n)  # log of expected posterior density
  Elp <- numeric(n)  # expectation of log posterior
  Vs <- numeric(n) # sample variance
  
  if(n > sum(!is.na(y))) stop("y must not contain missing values")
  
  if(class(bsim)=="logical") bsim <- sim(mod, n.sim=nsim) else{
    if(is.element(mclass, c("lm", "glm")))  nsim <- nrow(bsim@coef) else nsim <- nrow(bsim@fixef)
  }
  
  if(is.element(mclass, c("lm", "glm"))) family <- mod$family$family
  if(mclass=="lmerMod") family <- "gaussian"
  if(mclass=="glmerMod") family <- mod@resp$family$family
  if(is.null(family)) family <- "gaussian"
  if(!is.element(family, c("gaussian", "poisson", "binomial"))) stop("family must be an gaussian, binomial or poisson")
  
  
  # lppd, pwaic1, pwaic2
  predorig <- matrix(ncol=nsim, nrow=n) 
  linpred <- matrix(ncol=nsim, nrow=n)
  Xmat <- model.matrix(mod)
  if(is.element(mclass, c("lm", "glm"))) link <- mod$family$link 
  if(mclass=="lmerMod") link <- "identity"
  if(mclass=="glmerMod") link <- mod@resp$family$link
  if(is.null(link)|link=="identity") ilink <- "identity"
  if(link=="log") ilink <- "exp"
  if(link=="logit") ilink <- "plogis"
  
  for(r in 1:nsim){
    if(is.element(mclass, c("lm", "glm"))) linpred[,r] <- Xmat%*%bsim@coef[r,] else {
      bsimi <- matrix(bsim@fixef[r,], nrow=n, ncol=ncol(bsim@fixef), byrow=TRUE)
      colnames(bsimi) <- names(fixef(mod))
      nranefs <- length(ranef(mod))
      for(k in 1:nranefs){
        npars <- dim(bsim@ranef[[k]])[3]
        groupvar <- mod@flist[[k]]
        for(m in 1:npars){
          deltai <- bsim@ranef[[k]][r,,m]
          parname <- dimnames(bsim@ranef[[k]])[[3]][m]
          bsimi[,parname] <- bsimi[,parname] + deltai[match(groupvar, levels(groupvar))]
        }
      }
      for(i in 1:n){
        linpred[i,r] <- Xmat[i,]%*%bsimi[i,]
      }
    } # close ranef case 
  } # close r
  
  if(family=="gaussian"|family=="Gaussian"){
    for(r in 1:nsim) {
      predorig[,r] <- eval(parse(text=ilink))(linpred[,r])
    }
    for(j in 1:n) {
      lppdj[j] <- log((1/nsim) * sum(dnorm(y[j], predorig[j,], bsim@sigma)))
      lEp[j] <- log((1/nsim) * sum(dnorm(y[j], predorig[j,], bsim@sigma)))
      Elp[j] <- (1/nsim) * sum(log(dnorm(y[j], predorig[j,], bsim@sigma)))
      Vs[j] <- var(log(dnorm(y[j], predorig[j,], bsim@sigma)))
    }
  }
  if(family=="binomial"|family=="Binomial"){
    if(length(dim(y))>1) N <- apply(y, 1, sum) else N <- rep(1, n)
    if(length(dim(y))>1) y <- y[,1]
    for(r in 1:nsim) predorig[,r] <- eval(parse(text=ilink))(linpred[,r])
    for(j in 1:n) {
      lppdj[j] <- log((1/nsim) * sum(dbinom(y[j], prob=predorig[j,], size=N[j])))
      lEp[j] <- log((1/nsim) * sum(dbinom(y[j], prob=predorig[j,], size=N[j])))
      Elp[j] <- (1/nsim) * sum(log(dbinom(y[j], prob=predorig[j,], size=N[j])))
      Vs[j] <- var(log(dbinom(y[j], prob=predorig[j,], size=N[j])))      
    }
  }
  if(family=="poisson"|family=="Poisson"){
    for(r in 1:nsim) predorig[,r] <- eval(parse(text=ilink))(linpred[,r])
    for(j in 1:n) {
      lppdj[j] <- log((1/nsim) * sum(dpois(y[j], predorig[j,])))
      lEp[j] <- log((1/nsim) * sum(dpois(y[j], predorig[j,])))
      Elp[j] <- (1/nsim) * sum(log(dpois(y[j], predorig[j,])))  
      Vs[j] <- var(log(dpois(y[j], predorig[j,])))         
    }
  }  

  lppd <- sum(lppdj)
  pwaic1 <- 2*sum(lEp-Elp)
  pwaic2 <- sum(Vs)
  elppdwaic1 <- lppd - pwaic1
  elppdwaic2 <- lppd - pwaic2
  
  WAIC1 <- -2*elppdwaic1
  WAIC2 <- -2*elppdwaic2
  
  return(list(lppd=lppd, pwaic1=pwaic1, pwaic2=pwaic2, WAIC1=WAIC1, WAIC2=WAIC2))

}
