start.fit <- function(ts, model, xreg, start.control, linkfunc){
  ##############
  #Checks and preparations: 
  n <- length(ts)
  p <- length(model$past_obs)
  P <- seq(along=numeric(p)) #sequence 1:p if p>0 and NULL otherwise
  p_max <- max(model$past_obs, 0)
  q <- length(model$past_mean)
  Q <- seq(along=numeric(q)) #sequence 1:q if q>0 and NULL otherwise
  q_max <- max(model$past_mean, 0)
  r <- ncol(xreg)
  R <- seq(along=numeric(r)) #sequence 1:r if r>0 and NULL otherwise
  if(!(start.control$method %in% c("fixed", "iid", "GLM", "MM", "CSS", "ML", "CSS-ML"))) stop("Argument 'start.control$method' needs to be one of 'fixed', 'iid', 'GLM', 'MM',\n'CSS', 'ML' or 'CSS-ML'.")
  
  if(linkfunc=="identity") trafo <- function(x) x
  if(linkfunc=="log") trafo <- function(x) if(!is.null(x)) log(x+1) else NULL  
  param_start <- list(intercept=NULL, past_obs=NULL, past_mean=NULL, xreg=NULL)
  if(start.control$method == "fixed"){ #fixed values, use given ones where available
    param_start$intercept <- if(!is.null(start.control$intercept)) start.control$intercept else 1
    param_start$past_obs <- if(!is.null(start.control$past_obs)) start.control$past_obs else rep(0, p)
    param_start$past_mean <- if(!is.null(start.control$past_mean)) start.control$past_mean else rep(0, q)
    param_start$xreg <- if(!is.null(start.control$xreg)) start.control$xreg else rep(0, r)
  }else{
    # # # # # # #
    #Which observations to use for starting estimation?
        if(is.null(start.control$use)) start.control$use <- n
        if(length(start.control$use)<1 | length(start.control$use)>2) stop("Argument 'start.control$use' must be of length 1 or 2")
        if(length(start.control$use)==1){
          if(start.control$use==Inf) start.control$use <- n
          if(start.control$use<p+q+1) stop(paste("Too few observations for start estimation, argument 'start.control$use' must be\ngreater than p+q+1=", p+q+1, sep=""))
          if(start.control$use>n){ start.control$use <- n; warning(paste("Argument 'start.control$use' is out of range and set to the largest possible value n=", n, sep="")) }
          start_use <- 1:start.control$use
        }else{
          if(start.control$use[2]-start.control$use[1]<=p+q+1) stop(paste("Too few observations for start estimation, for argument 'start.control$use' the\ndifference start.control$use[2]-start.control$use[1] must be greater than p+q+1=", p+q+1, sep=""))
          if(start.control$use[2]>n | start.control$use[1]<1) stop(paste("Argument 'start.control$use' is out of range, start.control$use[1] must be\ngreater than 1 and start.control$use[2] lower than n=", n, sep=""))
          start_use <- start.control$use[1]:start.control$use[2]
        }
    ts_start <- ts[start_use]
    # # # # # # #
  }
  if(start.control$method == "iid"){
    param_start$intercept <- intercept <- mean(trafo(ts_start))
    param_start$past_obs <- rep(0, p) 
    param_start$past_mean <- rep(0, q)
    param_start$xreg <- rep(0, r)
  }
  if(start.control$method == "GLM"){
    max_delay <- max(model$past_obs, 0)
    delayed_ts <- function(del, timser=ts_start, max_del=max_delay) c(timser[((1+max_del):length(timser))-del])
    dataset <- data.frame(timser=delayed_ts(del=0), trafo(sapply(model$past_obs, delayed_ts)), if(r>0){apply(xreg[start_use, , drop=FALSE], 2, delayed_ts, del=0)}else{matrix(nrow=length(start_use)-max_delay, ncol=0)})
    startingvalues <- c(mean(trafo(ts_start)), rep(0, ncol(dataset)-1))
    glm_fit <- suppressWarnings(glm(timser ~ ., family=poisson(link=linkfunc), data=dataset, start=startingvalues)$coefficients)
    param_start$intercept <- intercept <- glm_fit[1]
    param_start$past_obs <- glm_fit[1+P] 
    param_start$past_mean <- rep(0, q)
    param_start$xreg <- glm_fit[1+p+R]
  }
  if(start.control$method %in% c("MM", "CSS", "ML", "CSS-ML")){ #approaches via an ARMA representation of the process, which differ only in the method to fit the ARMA process
    ts_start <- trafo(ts_start)
    k <- max(p_max, q_max)
    K <- seq(along=numeric(k)) #sequence 1:k if k>0 and NULL otherwise    
    if(start.control$method == "MM"){ #moment estimator via ARMA(1,1) representation, assume parameters for higher order to be zero
      if(k > 0){ #non-trivial case for q>0 and p>0
        momest <- momest_arma11(ts_start)
        ma <- c(momest["ma1"], rep(0,k-1)) #set higher order parameters to zero
        ar <- c(momest["ar1"], rep(0,k-1)) #see above
        intercept <- momest["intercept"]
        if(intercept < 0){
          warning("The start estimation of the intercept has been negative with\n argument 'method=\"MM\"'. Instead 'method=\"iid\"' is used.")
          ar[1] <- ma[1] <- 0
          intercept <- mean(ts_start)
        }  
      }else{
        ar <- ma <- NULL
        intercept <- mean(ts_start)
      }  
      regressors <- if(!is.null(start.control$xreg)) start.control$xreg else rep(0, r)  
    }
    if(start.control$method %in% c("CSS", "ML", "CSS-ML")){ #least squares or maximum likelihood estimator via ARMA(k,k) representation
      arma_fit <- as.numeric(suppressWarnings(arima(ts_start, order=c(k,0,k), xreg=xreg[start_use,], transform.pars=TRUE, method=start.control$method, optim.method=start.control$optim.method, optim.control=start.control$optim.control)$coef)) #Supress warning messages, which occur quite frequently and are not very relevant to the user, as this is only a start estimation. However, the interested user can find detailed information on this optimisation in the output.                     
      ma <- arma_fit[k+K]
      ar <- arma_fit[K]
      intercept <- arma_fit[k+k+1]
      regressors <- arma_fit[k+k+1+R] 
    }
    param_start$past_obs <- ar[model$past_obs]+ma[model$past_obs]
    param_start$past_mean <- -ma[model$past_mean] 
    param_start$intercept <- intercept*(1-sum(param_start$past_obs)-sum(param_start$past_mean))
    param_start$xreg <- regressors
  }
  result <- param_start
  return(result)
}
