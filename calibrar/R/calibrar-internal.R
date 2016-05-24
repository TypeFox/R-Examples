.setNewCalibration = function() {
  return(invisible())
}
.calculateObjetiveValue = function(obs, sim, info) {
  fit = NULL
  for(j in seq_len(nrow(info))) {
    if(!info$calibrate[j]) next
    var = info$variable[j]
    fit = c(fit, .fitness(obs=obs[[var]], sim=sim[[var]], FUN=info$type[j]))
  }
  names(fit) = info$variable[which(info$calibrate)]
  return(fit)
}

.fitness = function(obs, sim, FUN, ...) {
  FUN = match.fun(FUN)
  output = FUN(obs=obs, sim=sim, ...)
  return(output)
}

.messageByGen = function(opt, trace) {
  nx  = ceiling(log10(opt$control$maxgen))
  cat(sprintf(paste0("Generation %", nx, "d - Value: %.10g\n"), 
              opt$gen, trace$value[opt$gen]))
  #   print(trace$par[opt$gen, ])
}


.printSeq = function(n, preffix=NULL, suffix=NULL, sep="") {
  nx  = ceiling(log10(n))
  fmt = paste0(sprintf(paste0("%0", nx, "d"), seq_len(n)))
  out = paste(preffix, fmt, suffix, sep=sep)
  return(out)
}

.read.csv3 = function(path,  ...) {
  
  x = tryCatch(read.csv(file=path, row.names=1, ...),
               error = function(e) .errorReadingCsv(e, path))
  if(is.null(x)) stop()
  x = as.matrix(x)
  mode(x)="numeric"
  
  return(x)
}

.errorReadingCsv = function(e, path) {
  on.exit(stop())
  message(e, ".")
  message("Error reading file ", sQuote(path), ".", sep="")
  return(invisible())
}


.pasteList = function(x) {
  if(length(x)<=1) return(x)
  out = paste(paste(x[-length(x)], sep="", collapse=", "), "and", x[length(x)])
  return(out)
  
}

.messageBbyGen = function(opt, trace) {
  cat("--- Generation", opt$gen, "-----")
}

.checkActive = function(active, npar) {
  
  if(is.null(active)) return(rep(TRUE, npar))
  active = as.logical(active)
  if(all(!active)) stop("No parameter is active, at least one parameter must be optimized.")
  if(length(active)!=npar) stop("Length of 'phases' argument must match number of parameters.")
  return(active)
  
}

.checkPhases = function(phases, npar) {
  
  if(is.null(phases)) return(rep(1L, npar))

  phases = as.integer(phases)
  active = phases>0 & !is.na(phases)
  if(all(!active)) stop("No parameter is active, at least one parameter must be optimized.")
  if(length(active)!=npar) stop("Length of 'phases' argument must match number of parameters.")
  phases[phases<=0] = NA
  if(min(phases, na.rm=TRUE)!=1L) stop("No parameters has been specified for phase 1.")
  allPhases = seq(from=1, to=max(phases, na.rm=TRUE))
  ind = allPhases %in% phases
  missingPhases = .pasteList(allPhases[!ind])
  if(!all(ind)) stop(paste("No parameters for phases", missingPhases, "have been indicated."))
  
  return(phases)
  
}

.checkReplicates = function(replicates, nphases) {
  
  if(any(is.na(replicates))) stop("NAs are not allowed as 'replicates' number.")
  if(is.null(replicates)) return(rep(1, nphases))
  if(length(replicates)==1) return(rep(replicates, nphases)) 
  if(length(replicates)!=nphases) stop("Length of 'replicates' argument must match number of phases.")
  replicates = as.integer(replicates)
  
  if(any(replicates<=0)) stop("All 'replicates' must be positive integers.")
  
  return(replicates)
  
}

.checkBounds =  function(lower, upper, npar) {

  lf = uf = FALSE
  
  if(is.null(lower)) {
    lower = -Inf
    lf = TRUE
  }
  
  if(is.null(upper)) {
    upper = Inf
    uf = TRUE
  }
  
  nl = length(lower)
  nu = length(upper)
  
  if(nl==1) {
    if(npar!=1 & lower!=-Inf) 
      if(!lf) warning("Only one lower bound has been provided, used for all parameters.")
    lower = rep(lower, npar)
    nl = npar
  }
  
  if(nu==1) {
    if(npar!=1 & upper!=Inf) 
      if(!uf) warning("Only one upper bound has been provided, used for all parameters.")    
    upper = rep(upper, npar)
    nu = npar
  }
  
  if(nl!=npar) stop("Lower bounds must match parameter number.")
  if(nu!=npar) stop("Upper bounds must match parameter number.")
  
  if(any(is.na(lower))) {
    lower[is.na(lower)] = -Inf
    warning("NAs supplied in lower thresholds, replaced by -Inf")
  }
  
  if(any(is.na(upper))) {
    upper[is.na(upper)] = Inf
    warning("NAs supplied in upper thresholds, replaced by Inf")
  }
  
  if(any(lower >= upper)) stop("Lower bounds must be lower than upper bounds.")
  
  output = list(lower=lower, upper=upper)
  
  return(output)
  
}

.checkOpt = function(par, lower, upper) {
  
  if(any(par < lower, na.rm=TRUE)) 
    stop("Initial guess for parameters must be greater than lower threshold.")
  
  if(any(par > upper, na.rm=TRUE)) 
    stop("Initial guess for parameters must be lower than upper threshold.")
  
  ind = is.na(par) & !is.finite(lower) & !is.finite(upper)
  par[ind] = 0
  
  ind = is.na(par) & is.finite(lower) & is.finite(upper)
  par[ind] = 0.5*(lower + upper)[ind]
  
  ind = is.na(par) & is.finite(lower) & !is.finite(upper)
  par[ind] = lower[ind] + 0.1*abs(lower[ind])
  
  ind = is.na(par) & !is.finite(lower) & is.finite(upper)
  par[ind] = upper[ind] - 0.1*abs(upper[ind]) 
  
  return(par)
  
}

.optPopSize = function(n, selection) floor(0.5*(4 + floor(3*log(n)))/selection)

.checkControl = function(control, method, par, fn, active, skeleton, ...) {
  
  fn = match.fun(fn)
  
  con = list(trace = 0, fnscale = 1, parscale = rep.int(1L, length(which(active))), maxit = NULL, maxgen=NULL,
             abstol = -Inf, reltol = sqrt(.Machine$double.eps), REPORT = 10L, nCores=parallel::detectCores(), 
             alpha=0.05, age.max=1, selection=0.5, step=0.5, nvar=NULL, weights=1, sigma=NULL,
             method=method, aggFn=.weighted.sum, parallel=FALSE, run=NULL, master=NULL, useCV=TRUE,
             convergence=1e-6, stochastic=FALSE, verbose=FALSE, restart.file=NULL)
  
  
  popOpt = .optPopSize(n=length(which(active)), selection=con$selection)
  con$popsize = popOpt
  
  controlDef  = names(con)       # default options
  controlUser = names(control)   # user provided options
  
  con[controlUser] = control
  
  # get unknown options
  unknown = controlUser[!(controlUser %in% controlDef)]
  
  # update population size and selection rate
  if(con$popsize < popOpt) warning("'popsize' is too small, using default value.")
  con$popsize = max(con$popsize, popOpt)
  if(isTRUE(con$parallel)) {
    popOptP = ceiling(con$popsize/con$nCores)*con$nCores
    if(con$popsize != popOptP) message(paste("Optimizing 'popsize' to work with", con$nCores, "cores."))
    con$popsize = popOptP
    con$selection = round(con$selection*popOpt/popOptP, 1)
  }
  
  # check for user-provided variances (sigma)
  if(!is.null(con$sigma) & length(con$sigma)!=length(active)) 
    stop("Vector of variances (sigma) must match parameter length.")
  if(!is.null(con$sigma)) con$sigma = con$sigma[which(active)]

  if(inherits(fn, "objFn")) {
    con$nvar = attr(fn, "nvar")
    con$weights = attr(fn, "weights")
    names(con$weights) = attr(fn, "variables")
    if(con$verbose) {
      print(con$weights) 
    }
  }
  # check number of variables
  xpar = if(missing(skeleton)) par else relist(par, skeleton)
  if(is.null(con$nvar)) con$nvar = length(fn(xpar, ...))

  
  # update maximum number of function evaluations and generations
  if(!is.null(con$maxit) & !is.null(con$maxgen)) 
    warning("'maxit' and 'maxgen' provided, ignoring 'maxit'.")
  if(!is.null(con$maxit) & is.null(con$maxgen)) con$maxgen = floor(con$maxit/con$popsize)
  if(is.null(con$maxit) & is.null(con$maxgen)) con$maxgen = 2000L
  
  con$maxit = con$popsize*con$maxgen

  if(method=="default") {
    # update and check weights
    if(length(con$weights)==1) con$weights = rep(con$weights, con$nvar)
    if(length(con$weights)!=con$nvar) stop("Vector of weights should match the length of the output of fn.")
    if(any(con$weights<0)) stop("Weights should be positive numbers.")
    if(any(is.na(con$weights))) stop("Weights cannot be NA.")    
  }
  if(method=="cmaes") con$weights = NULL 
  # aggregation function for global fitness
  con$aggFn = match.fun(con$aggFn)
  
  if(length(unknown)!=0) 
    warning("Unknown control parameters: ", paste0(paste(unknown, collapse = ", ")),".")
  
  return(con)
  
}

.checkConvergence = function(control, nphases) {
  
  maxgen      = control$maxgen
  maxiter     = control$maxiter
  convergence = control$convergence
  
  if(length(maxgen)==1) maxgen = rep(maxgen, nphases)
  if(length(maxiter)==1) maxiter = rep(maxiter, nphases)
  if(length(convergence)==1) convergence = rep(convergence, nphases)
  
  if(!is.null(maxgen) & length(maxgen)!=nphases) 
    stop("'maxgen' length must match number of calibration phases.")
  
  if(!is.null(maxiter) & length(maxiter)!=nphases) 
    stop("'maxiter' length must match number of calibration phases.")
  
  if(!is.null(convergence) & length(convergence)!=nphases) 
    stop("'convergence' length must match number of calibration phases.")
  
  return(list(maxgen=maxgen, maxiter=maxiter, convergence=convergence))
  
}

.setWorkDir = function(run, i) {
  
  if(is.null(run)) return(invisible())
  work.dir = file.path(run, paste0("i", i))
  setwd(work.dir)
  return(work.dir)
  
}

.updateControl = function(control, opt, method) {
  # update sigma...
  return(control)
}
