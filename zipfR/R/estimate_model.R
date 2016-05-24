## internal: estimate parameters of LNRE model of specified class
estimate.model <- function (model, spc, param.names,
                            method, cost.function, m.max=15,
                            debug=FALSE, ...)
{
  UseMethod("estimate.model")
}

## generic estimation procedure for class 'lnre'
estimate.model.lnre <- function (model, spc, param.names,
                                 method, cost.function, m.max=15,
                                 debug=FALSE, ...)
{
  # vector of init values for parameters in transformed scale (same order as in param.names)
  param.values <- rep(0, length(param.names)) 

  # "Custom" not implemented for this LNRE model -> fall back to standard optimization method
  if (method == "Custom") {
    method <- if (length(param.names) > 1) "Nelder-Mead" else "NLM" # probably the best choices
  }
  
  compute.cost <- function (P.vector, param.names, model, spc, m.max=15, debug=FALSE)
    {
      P.trans <- as.list(P.vector)                     # translate parameter vector into list
      names(P.trans) <- param.names
      P <- model$util$transform(P.trans, inverse=TRUE) # convert parameters to normal scale
      model <- model$util$update(model, P)             # update model parameters (checks ranges)
      cost <- cost.function(model, spc, m.max)

      if (debug) {
        report <- as.data.frame(model$param)
        report$cost <- round(cost, digits=2)
        rownames(report) <- ""
        print(report)
      }
      cost
    }

  if (method == "NLM") {                # NLM = standard nonlinear minimization
    result <- nlm(compute.cost, param.values, print.level=debug, stepmax=10, steptol=1e-12,
                  param.names=param.names, model=model, spc=spc, m.max=m.max, debug=debug)

    res.code <- result$code
    if (res.code > 3) stop("parameter estimation failed (code ", res.code,")")
    if (res.code == 3) warning("estimated parameter values may be incorrect (code 3)")
    P.estimate <- as.list(result$estimate)
    names(P.estimate) <- param.names
  }
  else {                                # Nelder-Mead, SANN, BFGS = selected optim() algorithm
    result <- optim(param.values, compute.cost, method=method,
                    control=list(trace=debug, reltol=1e-12),
                    param.names=param.names, model=model, spc=spc, m.max=m.max, debug=debug)

    res.conv <- result$convergence
    if (res.conv > 1) stop("parameter estimation failed (code ", res.conv, ")")
    if (res.conv > 0)
      warning("iteration limit exceeded, estimated parameter values may be incorrect (code 1)")
    
    P.estimate <- as.list(result$par)
    names(P.estimate) <- param.names
  }
    
  model <- model$util$update(model, P.estimate, transformed=TRUE)
  model$gof <- lnre.goodness.of.fit(model, spc, n.estimated=length(param.names))
    
  model
}

