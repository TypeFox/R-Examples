#######################################################
##             compound Poisson GLM                  ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

cpglm <- function(formula, link = "log", data, weights, offset, 
                  subset, na.action = NULL, contrasts = NULL, 
                  control = list(), chunksize = 0, 
                  optimizer = "nlminb", ...) {

  call <- match.call()  
  if (missing(data)) 
    data <- environment(formula) 
  # use bigglm for big data sets 
  if (chunksize) {
    ans <- cpglm.profile.bigglm(call, data)
  } else {
    fr <- cpglm.mf(call, contrasts)
    link.power <- make.link.power(link)
    control <- do.call("cplm.control", control)
    ans <- cpglm.profile(fr, link.power, control, optimizer)  
  }
  ans@formula <- formula 
  ans@call <- call
  ans@na.action <- na.action
  ans@contrasts <- contrasts
  return(ans)
}


# function to implement the profile likelihood approach 
# Use cpglm.profile.bigglm for big data sets  
cpglm.profile <- function(fr, link.power = 0, control = list(), 
                          optimizer = "nlminb"){
  
  # compute the profile loglikelihood: parm = c(log(phi), p)
  llik_profile <- function(parm){
    fit2 <- cpglm.fit(fr, p = parm[2], link.power) 
    0.5 * dtweedie.nlogl(fr$Y, fit2$fitted.values, exp(parm[1]) / fr$wts, parm[2])
  }
      
  # generate starting values for p and phi
  init <- cpglm.init(fr, link.power)
  parm <- c(log(init$phi), init$p)

  # optimize the profile loglikelihood (cplm_optim is a wrapper)
  opt_ans <- cplm_optim(parm, llik_profile, gr = NULL, 
                   lower = c(-Inf, control$bound.p[1]),
                   upper = c(Inf, control$bound.p[2]),
                   control = control, optimizer = optimizer)  
  if (opt_ans$convergence) warning(opt_ans$message)
  
  # mle estimates
  p.max <- opt_ans$par[2]
  phi.max <- exp(opt_ans$par[1])
  # fit glm using the optimized index parameter
  fit <- cpglm.fit(fr, p.max, link.power)                  
    
  # return results
  out <- new("cpglm", 
             coefficients = fit$coefficients, residuals = fit$residuals,
             fitted.values = fit$fitted.values, weights = fit$weights,
             linear.predictors =  fit$linear.predictors,
             df.residual = as.integer(fit$df.residual),
             deviance = fit$deviance, call = call("foo"),
             aic = 2 * (opt_ans$value + fit$rank),           
             formula = ~ 1, control = control, contrasts = NULL,
             p = p.max, phi = phi.max, iter = fit$iter, converged = fit$converged,
             link.power= link.power, model.frame = fr$mf, na.action = NULL,
             offset = fr$offset, prior.weights = fit$prior.weights, y = fr$Y,
             inits = NULL, vcov = summary.glm(fit)$cov.scaled)
  return(out)  
}               


# function to implement the  profile likelihood approach for big data sets
cpglm.profile.bigglm <- function(call, data){
  # reconstruct the call
  if (is.null(control <- call$control)) 
    control <- list()  
  if (is.null(link <- call$link)) 
    link <- "log"
  if (is.null(optimizer <- call$optimizer)) 
    optimizer <- "nlminb" 
  control <- do.call("cplm.control", eval(control))
  link.power <- make.link.power(eval(link))
  mc <- match(c("link", "control", "family"), names(call), 0L)
  call <- call[-c(1, mc[mc>0])]
  pstart <- 1.5
  call$family <- tweedie(var.power = pstart, link.power = link.power)
  # default maxit in bigglm: the default seems not to be enough 
  call$maxit <- ifelse(is.null(call$maxit), 50, call$maxit)
  call$formula <- eval(call$formula)
  Y <- eval(call$formula[[2]], data, parent.frame())

  # profile loglikelihood 
  llik_profile <- function(parm){
    call$family <- tweedie(var.power = parm[2], link.power = link.power)
    fit2 <- do.call("bigglm", as.list(call))
    fs <- fitted(fit2, data)
    0.5 * dtweedie.nlogl(Y, fs$fitted.values, exp(parm[1]) / fs$prior.weights, parm[2])    
  }
  
  # generate starting values for phi  
  fit <- do.call("bigglm", as.list(call))  
  phistart <- fit$qr$ss / fit$df.resid 
  parm <- c(log(phistart), pstart)
  
  # optimize the profiled loglikelihood
  opt_ans <- cplm_optim(parm, llik_profile, gr = NULL,
                        lower = c(-Inf, control$bound.p[1]),
                        upper = c(Inf, control$bound.p[2]),
                        control = control, optimizer = optimizer)
  if (opt_ans$convergence) warning(opt_ans$message)
  p.max <- opt_ans$par[2]
  phi.max <- exp(opt_ans$par[1])
  
  # fit glm using the optimized index parameter
  call$family <- tweedie(var.power = p.max, link.power = link.power)
  fit <- do.call("bigglm", as.list(call))
  fs <- fitted(fit, data)  
    
  # return results
  out <- new("cpglm",  coefficients = coef(fit), 
             residuals = fs$residuals, fitted.values = fs$fitted.values,
             linear.predictors = fs$linear.predictors,
             weights = fs$weights, df.residual = as.integer(fit$df.resid),
             aic =  2 * (opt_ans$value + fit$n - fit$df.resid),           
             deviance = fit$deviance, call = call("foo"),
             formula = ~ 1, control = control,
             contrasts = NULL, p = p.max, phi = phi.max, 
             iter = fit$iterations, converged = fit$converged,
             link.power= link.power, model.frame = as.data.frame(NULL),
             na.action = NULL, offset = fs$offset,
             prior.weights = fs$prior.weights, y = Y,
             inits = NULL, vcov = vcov(fit))
  return(out)
}

