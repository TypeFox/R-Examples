
# estimation function for RRlog
RRlog.fit <- function(model, x, y, n.response, p, start, group, setPar2=-1, maxit=1000){
  logLik<- NA;
  coef<- rep(NA,ncol(x));
  iter<- NA;
  hessian<- matrix(NA,ncol=ncol(x),nrow=ncol(x))
  convergence <- NA;
  
  # parameter boundaries
  low <- rep(-Inf,ncol(x))
  up <- rep(Inf,ncol(x))
  if(is2group(model) & setPar2==-1){
    low <- c(low, 0)
    up <- c(up, 1)
  }
  
  try({est <- optim(par=start,fn=RRlog.loglik,
                    method="L-BFGS-B",
                    lower = low, upper = up,
                    control=list(fnscale=-1, maxit=maxit),hessian=T,  
                    X=x, y=y, prand=p, group=group, 
                    model=model, n.response=n.response, setPar2=setPar2)
  logLik=est$value;
  coef=est$par;
  iter=est$counts;
  hessian=est$hessian
  convergence <- est$convergence
  },silent=T)
  
  nams <- colnames(x)
  if(is2group(model)){
    nams <- c(nams, ifelse(model=="SLD", "t", ifelse(model == "UQTunknown", "pi", "gamma")))
  }
  #   print(grad)
  #   print(grad(func=RRlog.SLD.ll,x=est$par,cov=x,y=y,prand=p,group=group,setT=setT))
  res <- list(model=model,
              pString=paste("p = ",paste0(round(p,3),collapse=",")),
              coefficients=coef,
              logLik=logLik,param=nams,
              hessian=hessian,iter=iter, convergence=convergence)
  return(res)
}


# general loglik for RRlog
RRlog.loglik <- function(param, X, y, model, prand, group, n.response, setPar2=-1){
  # adjust t-parameter for SLD (or other 2-group models)
  if(is2group(model)){
    if (setPar2 != -1){
      beta <- param
      par2 <- setPar2
    }else{
      m <- length(param)
      beta <- param[1:(m-1)]
      par2 <- param[m]
    }
  }else{
    beta <- param
    par2 <- NULL
  }
  
  # linear prediction of latent states
  pi <- 1/(1+exp(-X%*%beta))
  
  # find correct randomization probability (maybe: speed up by appling only to unique groups)
  if(length(unique(group)) == 1){
    p <- getPW(model=model, p=prand, group=group[1], par2=par2)[2,]
    p1 <- pi*p["1"] + (1-pi)*p["0"]
  }else{
    # p <- t(sapply(group, function(gg) getPW(model=model, p=prand, group=gg, par2=par2)[2,]))
    p <- t(sapply(unique(group), function(gg) getPW(model=model, p=prand, group=gg, par2=par2)[2,]))
    idx <- match(group, unique(group))
    p1 <- p[idx,"1"]*pi + p[idx,"0"]*(1-pi)
  }
  # predicted probability of response=1
  loglik <- sum(dbinom(y, n.response, p1, log = TRUE))
  return(loglik)
}