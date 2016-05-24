loo.cv <- function(mod, nsim=100, bias.corr=FALSE){
  # Bayesian leave-one-out cross-validation according to Gelman et al. 2014
  # pp 175+
  # works for lm and glm-objects
  if(!is.element(class(mod)[1], c("lm", "glm"))) stop("mod must be an lm or glm object since leave-one-out cross-validation is only reliable when the observations are independent.")
  y <- mod$model[,1]
  n <- length(y)
  if(n > sum(!is.na(y))) stop("y must not contain missing values")
  
  logpostminusi.yi <- numeric(n)
  lppdji <- matrix(ncol=n, nrow=n)
  lppdj <- numeric(n)
  
  family <- mod$family$family
  if(is.null(family)) family <- "gaussian"
  link <- mod$family$link
  if(is.null(link)|link=="identity") ilink <- "identity"
  if(link=="log") ilink <- "exp"
  if(link=="logit") ilink <- "plogis"
  
  for(i in 1:n){
    dataminusi <- mod$model[-i,]
    if(class(dataminusi)!="data.frame"){
      dataminusi <- data.frame(dataminusi)
      names(dataminusi) <- names(mod$model)
    }
    modminusi <- update(mod, data=dataminusi)
    bsimminusi <- sim(modminusi, n.sim=nsim)
    predi <- numeric(nsim)
    predj <- numeric(nsim)
    Xmat <- model.matrix(mod)[i,]
    
    for(r in 1:nsim) predi[r] <-  eval(parse(text=ilink))(Xmat%*%bsimminusi@coef[r,])  #prediction for i
    

    
    if(family=="gaussian"|family=="Gaussian"){
      logpostminusi.yi[i] <- log(sum(dnorm(y[i], predi, bsimminusi@sigma))/nsim)
      if(bias.corr){
        for(j in 1:n){ 
          Xmatj <- model.matrix(mod)[j,]
          for(r in 1:nsim) predj[r] <-  eval(parse(text=ilink))(Xmatj%*%bsimminusi@coef[r,])  #prediction for i
          lppdji[i,j] <- log(sum(dnorm(y[j], predj, bsimminusi@sigma))/nsim)
        }
      }
    }
    
    if(family=="binomial"|family=="Binomial"){
      if(length(dim(y))>1) N <- apply(y, 1, sum) else N <- rep(1, n)
      if(length(dim(y))>1) y <- y[,1]
      logpostminusi.yi[i] <- log(sum(dbinom(y[i], prob=predi, size=N[i]))/nsim) 
      if(bias.corr){
        for(j in 1:n){ 
          Xmatj <- model.matrix(mod)[j,]
          for(r in 1:nsim) predj[r] <-  eval(parse(text=ilink))(Xmatj%*%bsimminusi@coef[r,])  #prediction for i
          lppdji[i,j] <- log(sum(dbinom(y[j], prob=predj, size=N[j]))/nsim)
        }
      }
      
    }
    if(family=="poisson"|family=="Poisson"){
      logpostminusi.yi[i] <- log(sum(dpois(y[i], predi))/nsim)
      if(bias.corr){
        for(j in 1:n){ 
          Xmatj <- model.matrix(mod)[j,]
          for(r in 1:nsim) predj[r] <-  eval(parse(text=ilink))(Xmatj%*%bsimminusi@coef[r,])  #prediction for i
          lppdji[i,j] <- log(sum(dpois(y[j], predj))/nsim)
        }
      }
      
    }
  }
  lppd_loo.cv <- sum(logpostminusi.yi)
  
  

  lppd_cloo.cv <- NULL # bias correctec cv

  bsim <- sim(mod, n.sim=nsim)
  predorig <- matrix(ncol=nsim, nrow=n)
    
  if(family=="gaussian"|family=="Gaussian"){
  for(r in 1:nsim) predorig[,r] <- eval(parse(text=ilink))(model.matrix(mod)%*%bsim@coef[r,])
    for(j in 1:n) {
      lppdj[j] <- log((1/nsim) * sum(dnorm(y[j], predorig[j,], bsim@sigma)))
    }
  }
  if(family=="binomial"|family=="Binomial"){
    for(r in 1:nsim) predorig[,r] <- eval(parse(text=ilink))(model.matrix(mod)%*%bsim@coef[r,])
    for(j in 1:n) {
      lppdj[j] <- log((1/nsim) * sum(dbinom(y[j], prob=predorig[j,], size=N[j])))
    }
  }
  if(family=="poisson"|family=="Poisson"){
    for(r in 1:nsim) predorig[,r] <- eval(parse(text=ilink))(model.matrix(mod)%*%bsim@coef[r,])
    for(j in 1:n) {
      lppdj[j] <- log((1/nsim) * sum(dpois(y[j], predorig[j,])))
    }
  }  
  if(bias.corr) lppdminusi_quer <- sum(apply(lppdji, 1, sum))/n
  lppd <- sum(lppdj)
  if(bias.corr) b <-lppd-lppdminusi_quer 
  if(bias.corr) lppd_cloo.cv <- lppd_loo.cv+b


  itscale <- -2*ifelse(bias.corr, lppd_cloo.cv, lppd_loo.cv)
  peff <- lppd - ifelse(bias.corr, lppd_cloo.cv, lppd_loo.cv) # estimate of effective number of parameters
  return(list(LOO.CV=lppd_loo.cv, bias.corrected.LOO.CV=lppd_cloo.cv, 
            minus2times_lppd=itscale,est.peff=peff))  
}