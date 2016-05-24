.calibrar = function (par, fn, gr = NULL, ..., lower = -Inf, upper = Inf, active=NULL, 
                    control = list(), hessian = FALSE, method = "default", skeleton=NULL) {
  
  skeleton = skeleton
  if(is.null(skeleton)) skeleton = as.relistable(par)
  
  npar = length(par)
 
  # check active parameters
  active = .checkActive(active=active, npar=npar)
  isActive = which(active)
  activeFlag = isTRUE(all(active))
  
  # update to active parameters only
  guess  = par
  par    = guess[isActive]
  lower  = lower[isActive]
  upper  = upper[isActive]
  
  npar = length(par)
  
  # closure for function evaluation
  fn   = match.fun(fn)
  control = .checkControl(control=control, method=method, par=guess, fn=fn, active=active, skeleton=skeleton, ...)
  
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    parx = relist(flesh = parx, skeleton = skeleton)
    fn(parx, ...)/control$fnscale
  }

  
  imethod = "default"
  optimMethods  = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                    "Brent")
  optimxMethods = c("nlm", "nlminb", "spg", "ucminf", "newuoa", "bobyqa", 
                    "nmkb", "hjkb", "Rcgmin", "Rvmmin")
    
  if(method %in% optimMethods) imethod = "optim"
  if(method %in% optimxMethods) imethod = "optimx"
  if(method == "cmaes") imethod = "cmaes"  
    
  output = 
    switch(imethod, 
           default = .optimES(par=par, fn=fn1, lower=lower, upper=upper, control=control, isActive=isActive),
           optim   = .optim(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                            hessian=hessian, method=method),
           optimx  = .optimx(par=par, fn=fn1, gr=gr, lower=lower, upper=upper, control=control, 
                             hessian=hessian, method=method),
           cmaes   = .cmaes(par=par, fn=fn1, lower=lower, upper=upper, control=control)
           )

  # reshaping full parameters
  paropt = guess
  paropt[isActive] = output$ppar 
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  
  # final outputs
  output = c(list(par=paropt), output, list(active=list(par=isActive, flag=activeFlag)))
  
  class(output) = c("optimES.result", class(output)) # change name
  
  return(output)
  
}


# wrapper for optim -------------------------------------------------------

.optim = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  if(!(method %in% c("L-BFGS-B", "Brent"))) {
    lower = -Inf
    upper = Inf
  }
  
  output = suppressWarnings(optim(par=par, fn=fn, gr=gr, method=method, lower=lower, 
              upper=upper, control=control, hessian=hessian))
  
  names(output)[names(output)=="par"] = "ppar"
  
  return(output)
  
}


# wrapper for optimx ------------------------------------------------------

.optimx = function(par, fn, gr, lower, upper, control, hessian, method) {
  
  parNames = names(par)
  
  out = suppressWarnings(optimx::optimx(par=par, fn=fn, gr=gr, method=method, lower=lower, 
                               upper=upper, control=control, hessian=hessian))
  
  par = as.numeric(out[, parNames])
  value = out$value
  counts = as.numeric(out[, c("fevals", "gevals")])
  output = list(ppar=par, value=value, counts=counts, convergence=out$convcode)
  
  return(output)
  
}


# wrapper for cma_es ------------------------------------------------------

.cmaes = function(par, fn, lower, upper, control) {
  
  npar = length(unlist(par))
  output = suppressWarnings(cmaes::cma_es(par=par, fn=fn, lower=lower, upper=upper, control=control))
  
  if(is.null(output$par)) {
    output$par = relist(rep(NA, npar), skeleton=par)
    if(!is.finite(output$value)) warning("Infinite value reached for fn.")
  }
  
  names(output)[names(output)=="par"] = "ppar"
  
  return(output)
  
}

# optimES internal --------------------------------------------------------

.optimES = function(par, fn, lower, upper, control, isActive, hessian=FALSE) {
  
  # get restart for the current phase
  restart = .restartCalibration(control) # flag: TRUE or FALSE
  if(isTRUE(restart)) {
    
    res = .getRestart(control=control)
    opt   = res$opt
    trace = res$trace
    
  } else {
    
    opt = .newOpt(par=par, lower=lower, upper=upper, control=control)
    
    trace = NULL
    
    if(control$REPORT>0 & control$trace>0) {
      
      trace = list()
      trace$control = control
      trace$par = matrix(NA, nrow=control$maxgen, ncol=length(isActive))
      trace$value = rep(NA, control$maxgen)
      trace$best  = rep(NA, control$maxgen)
      
      if(control$trace>1) {
        trace$sd = matrix(NA, nrow=control$maxgen, ncol=length(isActive))   
        trace$step = rep(NA, control$maxgen)     
      }
      
      if(control$trace>2) trace$opt = vector("list", control$maxgen)
      
    } 
    
  } 
  
  # start new optimization
  while(isTRUE(.continueEvolution(opt, control))) {
    
    opt$gen  = opt$gen + 1
    opt$ages = opt$ages + 1
    
    # create a new population
    
    if(all(opt$SIGMA==0)) break
    opt$pop = .createPopulation(opt)
    
    # evaluate the function in the population: evaluate fn, aggregate fitness
    
    opt$fitness = .calculateFitness(opt, fn=fn)
    
    # select best 'individuals'
    
    opt$selected = .selection(opt)
    
    # create the new parents: MU and SD
    
    opt = .calculateOptimalWeights(opt)
    opt = .updatePopulation(opt)
    
    # save detailed outputs
    if(control$REPORT>0 & control$trace>0) {
      
      trace$par[opt$gen, ] = opt$MU
      trace$best[opt$gen]  = opt$selected$best$fit.global
      
      if(control$trace>1) {
        trace$sd[opt$gen, ]  = opt$SIGMA
        trace$step[opt$gen]  = opt$step       
      }
      
      if(opt$gen%%control$REPORT==0) {
        trace$value[opt$gen] = control$aggFn(fn(opt$MU), control$weights)
        if(control$trace>2) trace$opt[[opt$gen]] = opt
      }
      
    }
    
    # save restart
    .createRestartFile(opt=opt, trace=trace, control=control)
    
    if(control$verbose & opt$gen%%control$REPORT==0) 
      .messageByGen(opt, trace)
    
  } # end generations loop
  
  partial = fn(opt$MU)
  value = control$aggFn(x=partial, w=control$weights) # check if necessary
  names(opt$MU) = names(par)
  opt$counts = c('function'=opt$gen*control$popsize, generations=opt$gen)
  
  output = list(ppar=opt$MU, value=value, counts=opt$counts, 
                trace=trace, partial=partial, convergence=1)
  
  return(output)
  
}




