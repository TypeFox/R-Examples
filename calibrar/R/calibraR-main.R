# calibrar package: Automated Calibration for Complex (Ecological) Models --------

#' Automated Calibration for Complex (Ecological) Models
#' 
#' Automated Calibration for Complex (Ecological) Models
#' 
#' @name calibrar-package
#' @aliases calibrar-package calibrar
#' @docType package
#' @author Ricardo Oliveros-Ramos Maintainer: Ricardo Oliveros-Ramos
#' <ricardo.oliveros@@gmail.com>
#' @references calibrar: an R package for the calibration of ecological models (Oliveros-Ramos and Shin 2014)
#' @keywords calibration
#' @examples
#' \dontrun{
#' require(calibrar)
#' set.seed(880820)
#' path = NULL # NULL to use the current directory
#' # create the demonstration files
#' demo = calibrarDemo(model="PoissonMixedModel", L=5, T=100) 
#' # get calibration information
#' calibrationInfo = getCalibrationInfo(path=demo$path)
#' # get observed data
#' observed = getObservedData(info=calibrationInfo, path=demo$path)
#' # read forcings for the model
#' forcing = read.csv(file.path(demo$path, "master", "environment.csv"), row.names=1)
#' # Defining 'runModel' function
#' runModel = function(par, forcing) {
#' output = calibrar:::.PoissonMixedModel(par=par, forcing=forcing)
#' # adding gamma parameters for penalties
#' output = c(output, list(gammas=par$gamma)) 
#' return(output)
#' }
#' # real parameters
#' cat("Real parameters used to simulate data\n")
#' print(demo$par)
#' # objective functions
#' obj  = createObjectiveFunction(runModel=runModel, info=calibrationInfo, 
#'                                observed=observed, forcing=forcing)
#' cat("Starting calibration...\n")
#' control = list(weights=calibrationInfo$weights, maxit=3.6e5) # control parameters
#' cat("Running optimization algorithms\n", "\t", date(), "\n")
#' cat("Running optim AHR-ES\n")
#' ahr = calibrate(par=demo$guess, fn=obj, lower=demo$lower, upper=demo$upper, control=control)
#' summary(ahr)
#' } 
NULL

# calibrate ---------------------------------------------------------------

#' @title Sequential parameter estimation for the calibration of models
#' @description This function performs the optimization of a function, possibly 
#' in sequential phases of increasing complexity, and it is designed for the 
#' calibration of a model, by minimizing the error function \code{fn} associated to it.  
#' @param par A numeric vector. The length of the par argument defines the 
#' number of parameters to be estimated (i.e. the dimension of the problem).
#' @param fn The function to be minimized.
#' @param gr the gradient of \code{fn}. Ignored, added for portability with
#' other optimization functions.
#' @param \dots Additional parameters to be passed to \code{fn}.
#' @param method The optimization method to be used. The 'default' method
#' is the AHR-ES (Oliveros & Shin, 2016). All the methods from stats::optim,
#' optimx::optimx and cmaes::cma_es are available.
#' @param lower Lower threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{-Inf}. By default \code{-Inf} is used (unconstrained).
#' @param upper Upper threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{Inf}. By default \code{Inf} is used (unconstrained). 
#' @param control Parameter for the control of the algorithm itself, see details.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Currently not implemented. 
#' @param phases An optional vector of the same length as \code{par}, 
#' indicating the phase at which each parameter becomes active. If omitted, 
#' default value is 1 for all parameters, performing a single optimization.
#' @param replicates The number of replicates for the evaluation of \code{fn}.
#' The default value is 1. A value greater than 1 is only useful for stochastic
#' functions.
#' @details In the control list, \code{aggFn} is a function to aggregate \code{fn} to 
#' a scalar value if the returned value is a vector. Some optimization algorithm can 
#' exploite the additional information provided by a vectorial output from \code{fn}.
#' @author Ricardo Oliveros-Ramos
#' @examples
#' calibrate(par=rep(NA, 5), fn=SphereN)
#' \dontrun{
#' calibrate(par=rep(NA, 5), fn=SphereN, replicates=3)
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=3, lower=-5, upper=5)
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=3, lower=-5, upper=5, phases=c(1,1,1,2,3))
#' calibrate(par=rep(0.5, 5), fn=SphereN, replicates=c(1,1,4), lower=-5, upper=5, phases=c(1,1,1,2,3))
#' }
#' @export
calibrate = function(par, fn, gr = NULL, ..., method = "default",
                     lower = NULL, upper = NULL, control = list(), 
                     hessian = FALSE, phases = NULL, replicates=1) {

  # check for a restart file
  restart = .restartCalibration(control, type="results")
  
  skeleton = as.relistable(par)
  par    = unlist(par)
  lower  = unlist(lower)
  upper  = unlist(upper)
  phases = unlist(phases)
  
  npar = length(par)
  
  # checking conformance of all arguments
  if(is.null(names(par))) names(par) = .printSeq(npar, preffix="par")
  
  fn = match.fun(fn)
  
  phases     = .checkPhases(phases=phases, npar=npar)
  bounds     = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess      = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)
  
  par     = guess
  lower   = bounds$lower
  upper   = bounds$upper
  nphases = max(phases, na.rm=TRUE)
  
  replicates = .checkReplicates(replicates, nphases) 
  
  output = if(isTRUE(restart)) .getResults(control=control) else list(phase=1)
  
  conv = .checkConvergence(control, nphases)
  
  # start the sequential parameter estimation
  for(phase in seq(from=output$phase, to=nphases)) {
    
    if(output$phase > nphases) break
  
    control$maxgen      = conv$maxgen[phase]
    control$maxiter     = conv$maxiter[phase]
    control$convergence = conv$convergence[phase]

    active = (phases <= phase) # NAs are corrected in .calibrar 
    # call optimizers handler .calibrar
    temp = .calibrar(par=par, fn=fn, gr = NULL, ..., method = method, 
                   lower = lower, upper = upper, control = control, 
                   hessian = hessian, active=active, skeleton=skeleton)
    
    output$phases[[phase]] = temp # trim?
    output$phase = phase + 1
    
    .createOutputFile(output, control) 
    
    par = temp$par #update parameter guess
    control = .updateControl(control=control, opt=temp, method=method)  # update CVs? 
    
    cat(sprintf("\nPhase %d finished (%d of %d parameters active)\n",
                phase, sum(active, na.rm=TRUE), npar))
    cat(sprintf("Function value: %g \n", temp$value))
    print(par[which(active)])
    cat("\n")
  }
  
   isActive = (phases>0) & !is.na(phases)
   paropt = output$phases[[nphases]]$par # save parameters of last phase
  
#   newNames = rep("*", npar)
#   newNames[isActive] = ""
#   
#   names(paropt) = paste0(names(paropt), newNames)

  paropt = relist(paropt, skeleton)
  
  final = list(par=paropt, value=output$phases[[nphases]]$value, 
               counts=output$phases[[nphases]]$counts, 
               partial=output$phases[[nphases]]$partial, 
               active=isActive, fn=fn)
  
  output = c(final, output)
  class(output) = c("calibrar.results")
  .createOutputFile(output, control) 
  
  return(output)
  
}

# optimES -----------------------------------------------------------------

#' @title Optimization using Evolutionary Strategies
#' @description This function performs the optimization of a function using 
#' evolutionary strategies, by default the AHR-ES (Oliveros & Shin, 2015). 
#' @param par A numeric vector. The length of the par argument defines the 
#' number of parameters to be estimated (i.e. the dimension of the problem).
#' @param fn The function to be minimized.
#' @param gr the gradient of \code{fn}. Ignored, added for portability with
#' other optimization functions.
#' @param \dots Additional parameters to be passed to \code{fn}.
#' @param lower Lower threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{-Inf}. By default \code{-Inf} is used (unconstrained).
#' @param upper Upper threshold value(s) for parameters. One value or a vector 
#' of the same length as par. If one value is provided, it is used for all 
#' parameters. \code{NA} means \code{Inf}. By default \code{Inf} is used (unconstrained). 
#' @param active A boolean vector of the same length of \code{par}. If \code{TRUE}, the
#' parameter is optimized, if \code{FALSE} the parameter is fixed to the value specified
#' in \code{par}.
#' @param control Parameter for the control of the algorithm itself, see details.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' Currently not implemented. 
#' @param method The optimization method to be used. Currently, the only implemented
#' is the 'default' method, corresponding to the AHR-ES (Oliveros & Shin, 2015).
#' @author Ricardo Oliveros-Ramos
#' @examples
#' optimES(par=rep(1, 5), fn=SphereN)
#' @export
optimES = function (par, fn, gr = NULL, ..., lower = -Inf, upper = Inf, active=NULL, 
                    control = list(), hessian = FALSE, method = "default") {
  
  
  npar = length(par)

  if(is.null(names(par))) names(par) = .printSeq(npar, preffix="par")
  
  active = .checkActive(active=active, npar=npar)
  bounds = .checkBounds(lower=lower, upper=upper, npar=npar)
  guess  = .checkOpt(par=par, lower=bounds$lower, upper=bounds$upper)
  
  isActive = which(active)
  activeFlag = isTRUE(all(active))
  
  par    = guess[isActive]
  lower  = bounds$lower[isActive]
  upper  = bounds$upper[isActive]
  
  npar = length(par)
  
  # closure for function evaluation
  fn   = match.fun(fn)
  
  control = .checkControl(control=control, method=method, par=par, fn=fn, active=active, ...)
  
  fn1  = function(par) {
    parx = guess
    parx[isActive] = par
    fn(parx, ...)/control$fnscale
  }
  
  output = .optimES(par=par, fn=fn1, lower=lower, upper=upper, control=control, isActive = isActive)
  
  paropt = guess
  paropt[isActive] = output$ppar 
  
  if(is.null(names(paropt))) names(paropt) = .printSeq(npar, preffix="par")
  
  output = list(par=paropt, output, active=list(par=isActive, flag=activeFlag))
  
  class(output) = c("optimES.result", class(output))
  
  return(output)
  
}

# getObservedData ---------------------------------------------------------

#' @title Get observed data for the calibration of a model 
#' 
#' @description Create a list with the observed data with the 
#' information provided by its main argument. 
#' 
#' @param info A data.frame with the information about the calibration, 
#' normally created with the \code{\link{getCalibrationInfo}} function. 
#' See details.
#' @param path Path to the directory to look up for the data.
#' @param data.folder folder in the path containing the data.
#' @param \dots Additional arguments to \code{read.csv} function 
#' to read the data files.
#' @return A list with the observed data needed for a calibration, to be used 
#' in combination with the \code{\link{createObjectiveFunction}}.
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{createObjectiveFunction}}, \code{\link{getCalibrationInfo}}.
#' @export
getObservedData = function(info, path, data.folder="data", ...) {
  
  observed  = list()
  variables = info$variable
  
  useData       = as.logical(info$useData)
  
  cat("Creating observed data list for calibration...","\n")
  
  for(var in 1:nrow(info)) {
    
    cat(paste0("Variable: ", variables[var], "\n"))
    var.path        = file.path(path, data.folder, paste0(variables[var],".csv"))
    datos           = if(useData[var]) .read.csv3(var.path, ...) else NA
    observed[[var]] = datos
    
  }
  
  names(observed) = variables
  
  return(observed)
  
}

# getCalibrationInfo ------------------------------------------------------

#' Get information to run a calibration using the \code{calibrar} package.
#' 
#' A wrapper for \code{read.csv} checking column names and data types 
#' for the table with the calibration information.
#' 
#' @param path The path to look for the file.
#' @param file The file with the calibration information, see details.
#' @param stringsAsFactors To be passed to \code{read.csv}.
#' @param \dots Additional arguments to \code{read.csv} function.
#' @return A data.frame with the information for the calibration of a 
#' model, to be used with the \code{\link{createObjectiveFunction}} 
#' and \code{\link{getObservedData}}.
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{createObjectiveFunction}}, \code{\link{getObservedData}}.
#' @export
getCalibrationInfo = function(path, file="calibrationInfo.csv", stringsAsFactors=FALSE, ...) {
  
  caliPath = file.path(path, file)
  calibrationInfo = read.csv(caliPath, stringsAsFactors=FALSE, ...)
  
  fullNames = c("variable", "type", "calibrate", "weights", "useData")  
  doesNotMatch = !(names(calibrationInfo) %in% fullNames)
  dnm = names(calibrationInfo)[doesNotMatch]
  
  isMissing = !(fullNames %in% names(calibrationInfo))
  im = fullNames[isMissing]
  
  sdnm = if(length(dnm)>1) " columns do " else " column does "
  sim  = if(length(im)>1) " variables are " else " variable is "
  msg1 = paste0("Error in ", caliPath, " file (", paste(sapply(dnm, sQuote), collapse=", "), 
                sdnm, "not match).")
  msg2 = paste0("Error in ", caliPath, " file (", paste(sapply(im, sQuote), collapse=", "), 
                sim, "missing).")
  
  if(any(doesNotMatch)) stop(msg1)
  if(any(isMissing)) stop(msg2)
  
  # cating correct data types
  calibrationInfo$variable  = as.character(calibrationInfo$variable)
  calibrationInfo$type      = as.character(calibrationInfo$type)
  calibrationInfo$calibrate = as.logical(calibrationInfo$calibrate)
  calibrationInfo$weights   = as.numeric(calibrationInfo$weights)
  calibrationInfo$useData   = as.logical(calibrationInfo$useData)
  
  return(calibrationInfo)
}

# createObjectiveFunction -------------------------------------------------

#' Create an objective function to be used with optimization routines
#' 
#' Create a new function, to be used as the objective function in the 
#' calibration, given a function to run the model within R, observed data 
#' and information about the comparison with data.
#' 
#' @param runModel Function to run the model and produce a list of outputs.
#' @param info A data.frame with the information about the calibration, 
#' normally created with the \code{\link{getCalibrationInfo}} function. 
#' See details.
#' @param observed A list of the observed variables created with the 
#' function \code{\link{getObservedData}}
#' @param aggFn A function to aggregate \code{fn} to a scalar value if the
#' returned value is a vector. Some optimization algorithm can explote the
#' additional information provided by a vectorial output from \code{fn}
#' @param aggregate boolean, if TRUE, a scalar value is returned using the 
#' \code{aggFn}.
#' @param \dots More arguments passed to the \code{runModel} function.
#' @return A function, integrating the simulation of the model and the 
#' comparison with observed data. 
#' @author Ricardo Oliveros-Ramos
#' @seealso \code{\link{getObservedData}}, \code{\link{getCalibrationInfo}}.
#' @export
createObjectiveFunction = function(runModel, info, observed, aggFn=.weighted.sum, 
                                   aggregate=FALSE, ...) {

  fn   = match.fun(runModel)
  aggFn = match.fun(aggFn)
  
  force(observed)
  force(info)
  force(aggregate)
  
  weights = info$weights[info$calibrate]

  # check for names in observed and simulated
  fn1  = function(par) {
    aggFn = match.fun(aggFn)
    simulated = fn(par, ...)
    # apply fitness to all outputs
    output = .calculateObjetiveValue(obs=observed, sim=simulated, info=info)
    if(isTRUE(aggregate)) output = aggFn(x=output, w=weights)
    return(output)
  }
  
  class(fn1) = c(class(fn1), "objFn")

  fnx = function(par) {
    simulated = fn(par, ...)
    return(simulated)
  }
  
  attr(fn1, "nvar") = sum(info$calibrate)
  attr(fn1, "weights") = weights
  attr(fn1, "variables") = info$variables[info$calibrate]
  attr(fn1, "fn") = fnx
  return(fn1) 
  
}
