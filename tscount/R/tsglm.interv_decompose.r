tsglm.interv_decompose <- function(ts, model, xreg, link, paramvec, isolate){
#Decomposes a time series in the effect of one or several covariates and an INGARCH model with the remaining covariates
  #isolate: Integer vector. Giving the indices of the covariates whose effect should be separated from the time series.
##############################

  ##############
  #Checks and preparations:
  n <- length(ts)
  p <- length(model$past_obs)
  P <- seq(along=numeric(p)) #sequence 1:p if p>0 and NULL otherwise
  p_max <- max(model$past_obs, 0)
  q <- length(model$past_mean)
  Q <- seq(along=numeric(q)) #sequence 1:q if q>0 and NULL otherwise
  q_max <- max(model$past_mean, 0)
  Q_max <- seq(along=numeric(q_max))
  r <- max(ncol(xreg), 0)
  R <- seq(along=numeric(r)) #sequence 1:r if r>0 and NULL otherwise
  param <- list( #transform parameter vector to a list
    intercept=paramvec[1],
    past_obs=paramvec[1+P],
    past_mean=paramvec[1+p+Q],
    xreg=paramvec[1+p+q+R]
  )    
  stopifnot(
    r>0,
    1<=isolate & isolate<=r
  )
  elim <- 1:r %in% isolate
  ##############
  
  ##############
  #Initialisation:
  #The code for this initialisation is taken from the functions ingarch.condmean and loglin.condmean. Do not modify it here without modifying it there.
  if(link == "identity"){
    #Initialisation by stationary solution (and its partial derivatives):
    denom <- 1-sum(param$past_obs)-sum(param$past_mean)    
    kappa_stationary <- param$intercept/denom
    kappa <- c(rep(kappa_stationary, q_max), numeric(n))  
    z <- c(as.integer(rep(round(kappa_stationary), p_max)), ts)
  }
  if(link == "log"){
    #Initialisation by arbritary zero:
    kappa <- c(rep(0, q_max), numeric(n))
    z <- c(as.integer(rep(0, p_max)), ts)
  }
  mu <- numeric(q_max+n)
  contamination <- integer(p_max+n)
  cleaned <- z
  X <- matrix(0, nrow=q_max+n, ncol=r) #regressors are initalised by zero
  X[q_max+(1:n), ] <- xreg  
  ##############  
  
  ##############
  #Recursion:
  if(link == "identity"){
    for(t in 1:n){
      kappa[t+q_max] <- param$intercept + sum(param$past_obs*z[(t-model$past_obs)+p_max]) + sum(param$past_mean*kappa[(t-model$past_mean)+q_max]) + sum(param$xreg*X[t+q_max, ]) - sum(param$past_mean*colSums(model$external*param$xreg*t(X[(t-model$past_mean)+q_max, , drop=FALSE])))
      mu[t+q_max] <- sum(param$past_obs*contamination[(t-model$past_obs)+p_max]) + sum(param$past_mean*mu[(t-model$past_mean)+q_max]) + sum(param$xreg[elim]*X[t+q_max, elim]) - sum(param$past_mean*colSums(model$external[elim]*param$xreg[elim]*t(X[(t-model$past_mean)+q_max, elim, drop=FALSE])))
      contamination[t+p_max] <- min(round(mu[t+q_max]/kappa[t+q_max]*z[t+p_max]),z[t+p_max]) #contamination may not be larger than the original observation   
    }
    cleaned <- z-contamination
  }
  if(link == "log"){
    for(t in 1:n){
      mu[t+q_max] <- sum(param$past_obs*log(1+contamination[(t-model$past_obs)+p_max]/(cleaned[(t-model$past_obs)+p_max]+1))) + sum(param$past_mean*mu[(t-model$past_mean)+q_max]) + sum(param$xreg[elim]*X[t+q_max, elim]) - sum(param$past_mean*colSums(model$external[elim]*param$xreg[elim]*t(X[(t-model$past_mean)+q_max, elim, drop=FALSE])))
      cleaned[t+p_max] <- round(z[t+p_max]/exp(mu[t+q_max]))
      contamination[t+p_max] <- z[t+p_max]-cleaned[t+p_max]
    }    
  }  
     
  result <- list(ts=ts, cleaned=cleaned[p_max+(1:n)], contamination=contamination[p_max+(1:n)])
  return(result)
} 
