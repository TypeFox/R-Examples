# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

#' @title Create Combination Matrix
#' @description Generates a matrix containing all the combination listings
#' @param n The number of variables.
#' @details Port expand.grid to C++ at a later time...
#' @return Returns a binary matrix (e.g. 1 or 0) entries
#' @keywords internal
comb.mat = function(n){
  c = rep(list(1:0), n)
  expand.grid(c)
}

#' @title TS Model Checks
#' @description Stops the process in R if there is an issue with the model desc
#' @param desc A \code{character} vector containing \code{ts.model} terms from desc. 
#' @details Checks if there are two or more objects found of type: DR, QN, RW, or WN. In addition, it currently forbids ARMA Models.
#' @return Returns nothing if there is no issues. Otherwise, it will kill the process. 
#' @keywords internal
select.desc.check = function(desc){
  models.active = count_models(desc)
  
  # Identifiability issues
  if(any( models.active[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, one of the supplied models will have identifiability issues. Please submit a new model.")
  }
  
  if(models.active["ARMA"] > 0){
    stop("Model selection with ARMA terms is NOT supported currently.")
  }
  
}

#' @title Formats the model score matrix
#' @description The model score matrix receives the appropriate model numbering, col descriptors, and ordering of models.
#' @param out         A \code{list} containing the model matrix
#' @param model.names A \code{character vector} that contains a list of the model names.
#' @return An updated matrix in the ith position of the list. 
#' @keywords internal
cust.model.score = function(out, model.names){

  colnames(out) = c("Obj Fun", "Optimism", "Criterion", "GoF P-Value")
  rownames(out) = sapply(model.names, FUN = paste0, collapse=" ")
  
  #out = out[order(out[,3]), ]

  out
}

#' @title Formats the rank.models (auto.imu) object
#' @description Creates the correct GMWM object and rank.models (auto.imu) summary object
#' @param out         A \code{list} containing the model matrix
#' @param model.names A \code{character vector} that contains the names of the models fit.
#' @param scales      A \code{vector} containing the 2^(1:J) scales.
#' @param N           A \code{int} indicating the length of the time series.
#' @param alpha       A \code{double} indicating the CI confidence.
#' @param robust      A \code{bool} indicating whether to use robust (T) or to use classical (F)
#' @param eff         A \code{double} indicating the efficiency for robust estimation.
#' @param B           A \code{int} to indicate how many bootstraps should occur when generating the V matrix.
#' @param G           A \code{int} to indicate how many guesses should be performed during the grid search
#' @param seed        A \code{seed} to recreate the same GMWM estimator results.
#' @param freq        A \code{double} that represents the frequency between observations.
#' @return An updated matrix in the ith position of the list. 
#' @keywords internal
output.format = function(out, model.names, scales, N, alpha, robust, eff, B, G, seed, freq = 1){
  desc = model.names[out[[1]][[2]]]
  
  model.ts = desc.to.ts.model(desc[[1]])
  
  out[[1]] = cust.model.score(out[[1]][[1]], desc)
  
  gmwm.obj = out[[2]]
  estimate = gmwm.obj[[1]]
  rownames(estimate) = model.ts$process.desc
  init.guess = gmwm.obj[[2]]
  rownames(init.guess) = model.ts$process.desc
  
  model.hat = model.ts
  
  model.hat$starting = F  
  
  
  if(any(model.hat$desc == "GM")){
    estimate[,1] = conv.ar1.to.gm(estimate[,1], model.hat$process.desc, freq)
    init.guess[,1] = conv.ar1.to.gm(init.guess[,1], model.hat$process.desc, freq) 
  }
  
  model.hat$theta = as.numeric(estimate)
  model.ts$theta = as.numeric(init.guess)
  
  
  # Release model
  out[[2]] = structure(list(estimate = estimate,
                                 init.guess = init.guess,
                                 wv.empir = gmwm.obj[[3]], 
                                 ci.low = gmwm.obj[[4]], 
                                 ci.high = gmwm.obj[[5]],
                                 orgV = gmwm.obj[[7]],
                                 V = gmwm.obj[[6]],
                                 omega = gmwm.obj[[12]],
                                 obj.fun = gmwm.obj[[11]],
                                 theo = gmwm.obj[[9]],
                                 decomp.theo = gmwm.obj[[10]],
                                 scales = scales, 
                                 robust = robust,
                                 eff = eff,
                                 model.type = "imu",
                                 compute.v = "fast",
                                 augmented = F,
                                 alpha = alpha,
                                 expect.diff = gmwm.obj[[8]],
                                 N = N,
                                 G = G,
                                 H = B,
                                 K = 1,
                                 model = model.ts, # ADD THIS IN AT A LATER TIME!
                                 model.hat = model.hat,
                                 starting = TRUE,
                                 seed = seed,
                                 freq = freq), class = "gmwm")
  
  out
}


#' @title Automatically select appropriate model for a set of models
#' @description 
#' Runs through a model selection algorithm to determine the best model in a given set
#' @param data       A \code{vector}, \code{data.frame}, \code{matrix}, or \code{gts} object with 1 column.
#' @param ...        Different \code{ts.model}s to be compared.
#' @param nested     A \code{bool} that indicates whether the ts.model objects are nested within a large object given within the list. If not, the a full model will be created.
#' @param bootstrap  A \code{bool} that is either true or false to indicate whether we use bootstrap or asymptotic By default, we use asymptotic.
#' @param model.type A \code{string} indicating whether the model should be a \code{"ssm"} or \code{"imu"}.
#' @param alpha      A \code{double} that indicates the level of confidence for the WV CI.
#' @param robust     A \code{boolean} that indicates whether to use robust estimation.
#' @param eff        A \code{double} between 0 and 1 that indicates the efficiency for the robust estimation.
#' @param B          A \code{integer} that contains the amount of bootstrap replications
#' @param G          A \code{integer} that indicates the amount of guesses for caliberating the startup.
#' @param seed       A \code{integer} that is used to set a seed for reproducibility.
#' @param freq       A \code{double} that represents the frequency between observations.
#' @details 
#' The models MUST be nested within each other. 
#' If the models are not nested, the algorithm creates the "common denominator" model.
#' 
#' To supply the models, enter them as:
#' AR1()+WN(), AR1(), 3*AR1()
#' 
#' Any parameter that you wish to use must then be specified.
#' e.g. to specify nested, you must use nested = T. Otherwise, it the function will stop.
#' 
#' Due to the structure of \code{rank.models}, you cannot mix and match \code{AR1()} and \code{GM()} objects.
#' So you must enter either AR1() or GM() objects. 
#' @return A \code{rank.models} object.
rank.models = function(data, ..., nested = F, bootstrap = F, 
                       model.type="ssm", alpha = 0.05, robust = F, eff = 0.6, B = 50, G = 100000, freq = 1, seed = 1337){
  
  set.seed(seed)
  
  models = list(...)
  numObj = length(models)
  desc = vector("list", numObj) 
  
  gmterm = FALSE
  
  for(i in 1:numObj){
    mod = models[[i]]
    
    if(!is.ts.model(mod)){
      stop("Ill-formed ... request. Detected non-ts.model object in ... Please specify parameters with names")
    }
    
    # Prevent mixing
    t1 = any(mod$desc == "GM")
    if(t1 == TRUE){ gmterm = TRUE }
    if(t1 && any(mod$desc == "AR1") ){
      stop("Please use either `GM()` or `AR1()` terms. Not both.")
    }
  
    select.desc.check(mod$desc)
    
    desc[[i]] = mod$desc
  }
  
  if(is.gts(data)){
    freq = attr(data, 'freq')
    
  }else if(is.data.frame(data) || is.matrix(data) || is.imu(data)){
    if(ncol(data) > 1){
      stop("`rank.models()` is not supported for multiple columns.")
    }
    if(is.imu(data)){
      freq = attr(data, 'freq')
    }
  }else if(is.lts(data)){
    data = data[,ncol(data)]
  }
  
  if(nested == F){
    full.str = .Call('gmwm_find_full_model', PACKAGE = 'gmwm', x = desc)
    
    if(!any(sapply(desc, function(x, want) isTRUE(all.equal(x, want)),  full.str)) ){
      print("Creating a Common Denominator Model!")
      desc[[length(desc)+1]] = full.str
    }
    desc = vector_to_set(desc)
  }else{
    
    full.str = models[[1]]$desc
    
    m = as.matrix(comb.mat(length(full.str)))
    m = m[-nrow(m),]
    
    desc = build_model_set(m, full.str)
    
  }
  
  out = .Call('gmwm_rank_models', PACKAGE = 'gmwm', data, model_str=desc, full_model=full.str, alpha, compute_v = "fast", model_type = model.type, K=1, H=B, G, robust, eff, bootstrap)
  
  N = length(data)
  nlevels =  floor(log2(N))
  scales = .Call('gmwm_scales_cpp', PACKAGE = 'gmwm', nlevels)
  
  out[[1]] = output.format(out[[1]], desc, scales, N, alpha, robust, eff, B, G, seed, freq)  
  
  class(out) = "rank.models"
  
  out
}

#' @title Automatically select appropriate model for IMU
#' @description Runs through a model selection algorithm to determine the best model
#' @param model     A \code{ts.model} object that is the largest model to be tested.
#' @param data      A \code{vector}, \code{matrix}, \code{data.frame}, or \code{imu} object with either 1, 3, or 6 columns. 
#' @param bootstrap A \code{bool} that is either true or false to indicate whether we use bootstrap or asymptotic By default, we use asymptotic.
#' @param alpha     A \code{double} that indicates the level of confidence for the WV CI.
#' @param robust    A \code{boolean} that indicates whether to use robust estimation.
#' @param eff       A \code{double} between 0 and 1 that indicates the efficiency for the robust estimation.
#' @param B         A \code{integer} that contains the amount of bootstrap replications
#' @param G         A \code{integer} that indicates the amount of guesses for caliberating the startup.
#' @param seed      A \code{integer} that controls the reproducibility of the auto model selection phase.
#' @return A \code{auto.imu} object.
#' @details 
#' The \code{auto.imu} object stores two important features for each signal:
#' \itemize{
#' \item{[[1]]}{A matrix containing model output}
#' \item{[[2]]}{The best \code{gmwm} object.}
#' }
#' To access it for each signal use:
#' \code{object[[i]][[1]]} or \code{object[[i]][[2]]}, where \eqn{i} denotes the signal.  
#' @author JJB
#' @examples 
#' \dontrun{
#' if(!require("imudata")){
#' install_imudata()
#' library("imudata")
#' }
#' 
#' data(imu6)
#' 
#' # Example 1
#' test1 = imu(imu6, gyros = 1:3, accels = NULL, axis = c('X', 'Y', 'Z'), freq = 100)
#' 
#' m = auto.imu(test1)
#' 
#' # Process 1's model table
#' m[[1]][[1]]
#' 
#' # Process 1's best fitting gmwm object
#' m[[1]][[2]]
#' 
#' }
auto.imu = function(data, model = 3*AR1()+WN()+RW()+QN()+DR(), bootstrap = F, alpha = 0.05, robust = F, eff = 0.6, B = 50, G = 100000, seed = 1337){
  
  # Check object
  if(!is.imu(data) ) {
    stop('Object must an imu object via imu()')
  }
  
  # Prevent mixing
  if(any(model$desc == "AR1") && any(model$desc == "GM") ){
    stop("Please use either GM() or AR1() terms. Not both.")
  }
  
  # Extract data for IMU Injection
  sensors = attr(data, 'sensor')
  num.sensor = attr(data, 'num.sensor')
  axis = attr(data, 'axis')

  # Set seed for reproducible results
  # Need to figure out a way to do set same seed on each model generation.
  set.seed(seed)
  
  # Set up data and models for automatic processing
  full.str = model$desc
  m = as.matrix(comb.mat(length(full.str)))
  m = m[-nrow(m),]
  
  out = .Call('gmwm_auto_imu', PACKAGE = 'gmwm', data, combs=m, full_model=full.str, alpha, compute_v = "fast", model_type = "imu", K=1, H=B, G, robust, eff, bootstrap)
  
  # Handle post processing
  
  # Get model names
  model.names = build_model_set(m, full.str) 
  
  # Get basic data info
  N = nrow(data)
  nlevels =  floor(log2(N))
  scales = .Call('gmwm_scales_cpp', PACKAGE = 'gmwm', nlevels)
  
  # Get set statements for Wenchao
  n.gyro = num.sensor[1]
  n.acc = num.sensor[2]
  
  a.gyro = 0
  a.acc = 0
  
  freq = attr(data, 'freq')
  for(i in 1:ncol(data)){
    obj = out[[i]]
    obj = output.format(obj, model.names, scales, N, alpha, robust, eff, B, G, seed, freq)
    
    obj.gmwm = obj[[2]]
    if(a.acc != n.acc){
      a.acc = a.acc + 1 
      obj.gmwm$sensor = "Gyroscope"
      obj.gmwm$axis = axis[a.acc]
    } else if(a.gyro != n.gyro){
      a.gyro = a.gyro + 1 
      obj.gmwm$sensor = "Accelerometer"
      obj.gmwm$axis = axis[a.gyro]
    }

    obj.gmwm$num.sensor = num.sensor
    
    obj[[2]] = obj.gmwm
    out[[i]] = obj
  }
  
  class(out) = c("auto.imu","rank.models")
  out  
}

#' Print function for rank.models object
#' 
#' Prints the rank.models function nicely.
#' 
#' @param x   A \code{rank.models} object.
#' @param ... Additional parameters
#' @method print rank.models
#' @export
#' @keywords internal
print.rank.models = function(x, ...){
  summary.rank.models(x)
}

#' Print function for auto.imu object
#' 
#' Prints the auto.imu function nicely.
#' 
#' @param x   A \code{auto.imu} object.
#' @param ... Additional parameters
#' @method print auto.imu
#' @export
#' @keywords internal
print.auto.imu = function(x, ...){
  summary.auto.imu(x)
}

#' Summary function for auto.imu object
#' 
#' Provides a summary of the auto.imu results.
#' 
#' @param object A \code{auto.imu} object.
#' @param digits A \code{int} indicating how the numbers should be rounded.
#' @param ...    Additional parameters
#' @method summary auto.imu
#' @export
summary.auto.imu = function(object, digits = 4, ...){
  
  n.process = length(object)
  
  cat("There were ", n.process,"processes observed.\n\n")

  for(i in 1:n.process){
    out = object[[i]][[1]]
    cat(paste0("The model ranking for data column ", i, ": \n"))
    rownames(out) = paste0(1:nrow(out), ". ", rownames(out) )
    
    print(round(out,digits))
    
    cat("\n")
  }
}

#' Summary function for rank.models object
#' 
#' Provides a summary of the rank.models results.
#' 
#' @param object A \code{rank.models} object.
#' @param digits A \code{int} indicating how the numbers should be rounded.
#' @param ...    Additional parameters
#' @method summary rank.models
#' @export
summary.rank.models = function(object, digits = 4, ...){
  cat("The model ranking is given as: \n")
  out = object[[1]][[1]]
  
  rownames(out) = paste0(1:nrow(out), ". ", rownames(out) )
  
  print(round(out,digits))
}