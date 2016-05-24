
###############################################################################
# DivE - Diversity Estimator
# Authors: Daniel Laydon, Aaron Sim, Charles Bangham, Becca Asquith
# Date: 01 May 2014
# Version: 1.00
# 
# Contents:
#   A. SUBSAMPLE AND RAREFACTION DATA GENERATION
#     1. Other functions
#     2a-2c. Subsample sizes generation functions
#     3a-3i. Rarefaction and diversity data generation functions (incl. summary)
#     4a.    Curvature function
#   B. CURVE FITTING
#     1a-1g. Supporting fitting functions
#     2a-2b. Fitting functions
#     3a-3e. Fitting output functions (incl. summary)
#   C. SCORING AND MODEL COMPARISON 
#     1a-1f. Rounding functions
#     2a.    Scoring single model 
#     3a-3e. Scoring single model output
#     4a.    Comparison of models
#     5a-5e. Comparison of models output
#     6a.    Combining separate runs
#     7a.    Diversity estimate for different population size
###############################################################################

######################## A. SUBSAMPLE DATA GENERATION #########################
######################## A1. Function to format input #########################
# Two formats accepted:
# i. dataframe with two variables - clone id, number of clones
# ii. dataframe, vector or matrix with one dimension
format.input <- function(x, ...) {
  if (is.data.frame(x) && (dim(x)[2] == 2)) {
    as.character(rep(x[,1], x[,2]))
  } else if (is.vector(x)) {
    as.character(x)
  } else if ((is.data.frame(x) && dim(x)[2] == 1) || (is.matrix(x) && dim(x)[2] == 1)) {
    as.character(x[,1])
  } else {
    print("Incorrect sample input format - reformat and try again")
  }
}

############## A2a. Function to generate vector of subsample lengths #############
# Inputs: Main sample, desired number of subsamples
# By default produces a vector of length six with equally spaced subsample sizes
gen.subsamp.lengths <- function(main.samp, num.subsamp) {
  if ((length(num.subsamp))!=1) {
    sort(num.subsamp, decreasing=TRUE)
  } else {
    samp.int <- length(main.samp)/num.subsamp
    sort(round(seq(from=samp.int, to=length(main.samp), by=samp.int)), decreasing=TRUE)
  }
}

############### A2b. Generic generate subsample lengths function #################
divsamplenum <- function(ms, n) {UseMethod("divsamplenum")}

############### A2c. Default generate subsample lengths function #################
divsamplenum.default <- function(ms, n=6) { 
  if ((length(n)>1) && (n<2)) { 
    stop('Number of nested subsamples (subsizes) must be 2 or greater')
  }
  main.samp <- format.input(ms)
  gen.subsamp.lengths(main.samp=main.samp, num.subsamp=n)
}

############ A3a. Function to generate vector of rarefaction datapoints ############
# Inputs: Any sample, desired interval between rarefaction points, maximum data size of subsample, minimum data size of subsample
gen.rarefac.lengths <- function(samp, num.rarefac, max.rarefac=length(format.input(samp)), min.rarefac=1) {
  rarefac.int <- num.rarefac
  out <- round(seq(from=min.rarefac, to=max.rarefac, by=rarefac.int))
  if (out[length(out)]==max.rarefac) {
    return (out)
  } else {
    out <- c(out, max.rarefac)
    return (out)
  }
}

############ A3b. Function to determine diversity of a (sub)sample  ###########
# Inputs: Any vector
div.count <- function(samp) {length(unique(samp))}

##### A3c. Function to generate the sample diversity values at rarefaction datapoints for fitting ####
# Inputs: Sample, interval between rarefaction datapoints, minimum data size of subsample, Iterations, maximum data size of subsample 
gen.rarefac.div <- function(samp, num.rarefac, min.rarefac=1, B, max.rarefac=length(format.input(samp))) {
  l <- length(samp)
  rarefac.lengths <- gen.rarefac.lengths(samp, num.rarefac, max.rarefac, min.rarefac)
  s <- length(rarefac.lengths)
  order.vec <- as.vector(apply(matrix(sample(1:(l*B)), nrow=B, ncol=l), 1, order))
  shuf.mat <- matrix(samp[order.vec], nrow=B, ncol=l, byrow=TRUE)
  rarefac.div.mean <- rep(1, s)
  rarefac.div.sd <- rep(0, s)
  
  for (i in (2:s)) {
    temp <- apply(shuf.mat[,1:rarefac.lengths[i]], 1, div.count)
    rarefac.div.mean[i] <- mean(temp)
    rarefac.div.sd[i] <- sd(temp)
  }
  #rarefac.div.mean[s] <- div.count(samp)
  #rarefac.div.sd[s] <- 0
  list(RarefacXAxis=rarefac.lengths, 
       RarefacYAxis=rarefac.div.mean,
       div_sd=rarefac.div.sd,
       NResamples=B)
}

############## A3d. Generic generate subsamples/diversity function ############
divsubsamples <- function(mainsamp, nrf, minrarefac=1, maxrarefac=length(format.input(mainsamp)), NResamples=1000) {UseMethod("divsubsamples")}

############# A3e. Default generate subsamples/diversity function #############
divsubsamples.default <- function(mainsamp, nrf, minrarefac=1, maxrarefac=length(format.input(mainsamp)), NResamples=1000) { 
  samp <- format.input(mainsamp)

  xyvals <- gen.rarefac.div(samp=samp, num.rarefac=nrf, min.rarefac=minrarefac, B=NResamples, max.rarefac=maxrarefac)

  xyvals$call <- match.call()
  class(xyvals) <- "divsubsamples"
  xyvals
}

################ A3f. Print method -  subsample/diversity data ################
print.divsubsamples <- function(x, ...) {
  cat("RarefacXAxis:\n")
  print(x$RarefacXAxis)
  cat("\nRarefacYAxis (Mean diversity values):\n")
  print(x$RarefacYAxis)
}

################ A3g. Summary method: subsample/diversity data ################
summary.divsubsamples <- function(object, ...) {
  tot.div <- object$RarefacYAxis[length(object$RarefacYAxis)]
  tot.rarefac <- length(object$RarefacXAxis)
  samp.size <- object$RarefacXAxis[tot.rarefac]
  num.iter <- object$NResamples
  ave.sdpc <- mean(object$div_sd/object$RarefacYAxis)
  TAB <- cbind(Subsample.size=samp.size,
               Subsample.diversity=tot.div,
               No.of.rarefac.points=tot.rarefac,
               Iterations=num.iter,
               Ave.StdErr=ave.sdpc)
  dss.sum <- list(call=object$call,
              rsum=TAB)
  class(dss.sum) <- "summary.divsubsamples"
  dss.sum
}

############ A3h. Print summary method -  subsample/diversity data ############
print.summary.divsubsamples <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nSubsample data summary:\n")
  print(x$rsum)
}

############ A3i. Method to extract the relevant nested subset of a divsubsamples object ############
# Inputs: divsubsample object, size of nested subsample desired
divsubsamples_nested <- function(dss, maxsamp) {
  dss.temp <- dss
  end.index <- which.min(abs(maxsamp-dss$RarefacXAxis))
  dss.temp$RarefacXAxis <- dss$RarefacXAxis[1:end.index]
  dss.temp$RarefacYAxis <- dss$RarefacYAxis[1:end.index]
  dss.temp$div_sd <- dss$div_sd[1:end.index]
  return (dss.temp)
}

############ A4a. Method to calculate the curvature of a rarefaction curve ############
# Inputs: divsubsample object
Curvature <- function(dss) {
  if ((length(dss) > 1) && (class(dss[[1]])=="divsubsamples")) {
    SubSizes = c()
    for (Sub in 1:length(dss)) {
      SubSizes[Sub] = tail(dss[[Sub]]$RarefacXAxis,1)
    }
    LargestSub = order(SubSizes)[1]
    dss = dss[[LargestSub]]
  }
  x_val <- dss$RarefacXAxis
  y_val <- dss$RarefacYAxis

  x_end <- tail(x_val, 1)
  y_end <- tail(y_val, 1)
  
  AnB <- 0.5*x_end*y_end
  
  x_vec <- c(x_val[1], diff(x_val))
  
  y_vec_left <- c(0,y_val)[1:(length(y_val))]
  y_vec_right <- y_val[1:length(y_val)]
  y_vec <- (y_vec_right + y_vec_left)/2
  
  A <- sum(x_vec*y_vec) - AnB
  
  output <- A/AnB
  return (output)
}

############################# B. CURVE FITTING ################################

#require(deSolve) # To be packaged correctly
#require(FME) # To be packaged correctly

####### B1a. Function to Load model list, initial parameter guesses, parameter ranges #######
## Input: None
#load.Mod <- function() { # Give the index of the models to fit
#  load("modelset")
#}

####### B1b. Function to predict the diversity values from subsample x-values according to a given model #######
# Input: Model id, rarefaction x-axis values, model parameters
# Output: dataframe of x and y axis plotting values
pred.div <- function(model, rarefac, params) {
  model.param <- as.list(params)
  div.temp <- model(rarefac, model.param)
  if(any(div.temp=="NaN")) {  				# Ensures an infinite ssr
    div.temp <- rep(Inf, length(rarefac)) 	
  }
  data.frame(RarefacXAxis=rarefac, RarefacYAxis=div.temp)
}

############# B1c. Function to calculate the residuals (modCost) ##############
# Input: Model id, model parameters, divsubsamples object
model.cost <- function(model, params, dss) {
  xypred <- pred.div(model, dss$RarefacXAxis, params)
  xyactual <- data.frame(RarefacXAxis=dss$RarefacXAxis, RarefacYAxis=dss$RarefacYAxis)
  modCost(xypred, xyactual, x="RarefacXAxis")
}

####### B1d. Function to calculate the discrepancy as the mean of pointwise percentage error ########
# Input: Model id, model parameters, divsubsamples object
model.cost.abs <- function(model, params, dss) {
  ypred <- pred.div(model, dss$RarefacXAxis, params)$RarefacYAxis
  yactual <- dss$RarefacYAxis
  mean(abs((ypred-yactual)/yactual))
}

############ B1e. Function to convert parameter ranges to vectors #############
# Input: model.set, Param.ranges (possibly truncated)
convert.ranges <- function(model.set, param.ranges) {
  bounds.lower <- list()
  bounds.upper <- list()
  for (mod in 1:length(model.set)) {
    bounds.lower[[names(model.set)[mod]]] = param.ranges[[mod]]["Lower",]
    bounds.upper[[names(model.set)[mod]]] = param.ranges[[mod]]["Upper",]
  }
  list(lower=bounds.lower, upper=bounds.upper)
}

################### B1f. Function to perform a global Fit #####################
# Input: Function, initial parameters, lower parameter bounds, upper parameter bounds, Max no. of iterations, minimum variance
#global.fit <- function(func, init.param, lower.bd, upper.bd, numit, varleft) {
global.fit <- function(func, init.list, lower.bd, upper.bd, numit, varleft) {
  cat("Choosing optimal initial parameters for global fit", "\n")
  init.param <- (lower.bd + upper.bd)/2
  cost.lowest <- func(init.param)$model
  
  for (i in 1:dim(init.list)[1]) {
    if (all((init.list[i,] >= lower.bd)) & all((init.list[i,] <= upper.bd))) {
      cost.temp <- func(init.list[i,])$model
      if (cost.temp < cost.lowest) {
        cost.lowest <- cost.temp
        init.param <- init.list[i,]
      }
    }
  }
  cat("Performing global fit", "\n")
  modFit(f = func, p = init.param, lower=lower.bd, upper=upper.bd,       # "Pseudo" Global Fit
         method = "Pseudo", control = c(numiter = numit, varleft=varleft, verbose=TRUE))
}

#################### B1g. Function to perform a local Fit #####################
# Input: Function, initial parameters, lower parameter bounds, upper parameter bounds, Max no. of iterations, minimum variance
local.fit <- function(func, init.param, lower.bd, upper.bd) {
  cat("Performing local fit", "\n")
  localfit <- modFit(f = func, p = init.param, lower=lower.bd, upper=upper.bd,       # Local Fit
                     method = "Marq")
  return(localfit)
}

############# B1h. Function to perform a combined global-local Fit ############
# Input: Function, initial parameters, lower parameter bounds, upper parameter bounds, Max no. of iterations, minimum variance
combined.fit <- function(func, init.param, lower.bd, upper.bd, numit, varleft) {
  g.fit <- global.fit(func, init.param, lower.bd, upper.bd, numit, varleft)
  local.fit(func, g.fit$par, lower.bd, upper.bd)
}

####### B1i. Function to convert function "funk" to function of a single variable #######
# Input: function, parameters
single.var <- function(funk, params) {
  fn = function(x)
  {
    funk(x, params)
  }
  return(fn)
}

############# B1j. Function to find intersections between curves ##############
# Input: function, parameters1, parameter2, lower x-limit, upper x-limit
crossings <- function(model, param1, param2, xvals) {
  yvals1 <- model(xvals, param1)
  yvals2 <- model(xvals, param2)
  s1 <- try(SpatialLines(list(Lines(list(Line(cbind(xvals , yvals1 ))) , ID=1))), silent = TRUE)
  s2 <- try(SpatialLines(list(Lines(list(Line(cbind(xvals , yvals2 ))) , ID=1))), silent = TRUE)
  
  temp.int <- gIntersection(s1, s2)
  int.len <- length(temp.int)

  if (int.len != 0) {
    row.names(temp.int) <- seq(1,int.len)
  }

  intersections <- try(as.data.frame(temp.int), silent=TRUE)
  #intersections <- try(as.data.frame(gIntersection(s1, s2)), silent=TRUE)
  intersections
}


############### B1k. Function to evaluate distance between curves #############
# Input: function, parameters1, parameter2, lower x-limit, upper x-limit
simm <- function(model, param1, param2, lowerlimit, upperlimit, tot.length) {
  intersections <- crossings(model, param1, param2, seq(lowerlimit, upperlimit, length.out=tot.length))
  
  # No intersections
  if (length(intersections) == 0  |  (class(intersections) == "try-error")) {
    dist <- abs(as.numeric(as.character(integrate(single.var(model, param1), lower = lowerlimit, upper = upperlimit))[1])  -
      as.numeric(as.character(integrate(single.var(model, param2), lower = lowerlimit, upper = upperlimit))[1]))
    return(dist)
  } else if (length(intersections) != 0 & (class(intersections) != "try-error")) {
    # with intersections
    limits <- c(lowerlimit)
    for (i in (1:(dim(intersections)[1]))) {
      limits <- c(limits, intersections[i,1])
    }
    limits <- c(limits, upperlimit)
    
    dist <- 0
    for(interval in 1:(dim(intersections)[1]+1)) {
      dist <- dist + abs(as.numeric(as.character(integrate(single.var(model, param1), lower=limits[interval],
                                                           upper=limits[interval+1]))[1]) -
                                                             as.numeric(as.character(integrate(single.var(model, param2), lower=limits[interval], 
                                                                                               upper=limits[interval+1]))[1]))
    }
    return(dist)
  }
}

######## B1l. Function to evaluate the area under the reference curve #########
# Input: function, parameters1, lower x-limit, upper x-limit
ref.area <- function(model, param1, lowerlimit, upperlimit) {
  total.area <- abs(as.numeric(as.character(integrate(single.var(model, param1), lower = lowerlimit, upper = upperlimit))[1]))
  return(total.area)
}

####### B1l. Function to evaluate the mean squared distances between curves #######
# Input: function, parameters1, parameter2, lower x-limit, upper x-limit
msr.curves <- function(dist) {
  return(sum(dist*dist)/(2*nrow(dist)))
}

####### B1m. Function to evaluate the monotoniticy of the fitted curve ########
# Input: Model, parameters, minimum datapoint, maximum datapoint, number of points to test
is.mono <- function(model, params, min.samp, max.samp, tot.length) {
  xvals <- seq(min.samp, max.samp, length.out=tot.length)
  all(diff(model(xvals, params))>=0)
}

####### B1n. Function to evaluate the negativity of the 2nd derivative of the fitted curve ########
# Input: Model, parameters, minimum data size, maximum data size, number of points to test
is.slowing <- function(model, params, min.samp, max.samp, tot.length) {
  xvals <- seq(min.samp, max.samp, length.out=tot.length)
  all(round(diff(diff(model(xvals, params))), 5)<=0)
}

################# B2a. Fit a single sample to a single model ##################
fitsample <- function(model, init.param, lower.bd, upper.bd, dss, numit, varleft) {
  func <- function(p) {
    model.cost(model, p, dss)
  }
  combined.fit(func, init.param, lower.bd, upper.bd, numit, varleft)
}

############# B2b. Fit all samples to a single model ##############
# Input: SS=list of subsample sizes, model, initial parameters, parameter boundaries,
#       list of diversity subsample objects, numit, varleft, population total estimate, vector of original samples,
#       minimum sample size (for criteria 3 & 4)
fitallSubs <- function(SS, model.list, init.param, param.range, dSS, numit, varleft, tot.pop, v1, minsampsize) {
  
  ###### Extract model and name ######
  model <- model.list[[1]]
  model.name <- names(model.list)
  
  ###### Create lower and upper boundaries ######
  lower.bd <- param.range[1,]
  upper.bd <- param.range[2,]
  
  ###### Create storage matrices ######
  # Parameter names
  param.mat <- matrix(nrow=length(SS), ncol=ncol(init.param))
  rownames(param.mat) <- as.character(SS)
  colnames(param.mat) <- names(init.param[1,])
  
  # Criteria 1: Sum of squares
  ssr.mat <- matrix(nrow=length(SS), ncol=1)
  rownames(ssr.mat) <- as.character(SS)
  colnames(ssr.mat) <- c("Sum_of_squares")  
  
  # Criteria 1: Mean squared residual
  ms.mat <- matrix(nrow=length(SS), ncol=1)
  rownames(ms.mat) <- as.character(SS)
  colnames(ms.mat) <- c("Mean_squared_residual")
  
  # Criteria 1: Discrepancy - mean abs error
  discrep.mat <- matrix(nrow=length(SS), ncol=1)
  rownames(discrep.mat) <- as.character(SS)
  colnames(discrep.mat) <- c("Mean_absolute_error")
  
  # Criteria 2: Local diversity prediction
  local.mat <- matrix(nrow=length(SS), ncol=1)
  rownames(local.mat) <- as.character(SS)
  colnames(local.mat) <- c("Predicted_local_diversity")
  
  # Criteria 2: Global diversity prediction
  global.mat <- matrix(nrow=length(SS), ncol=1)
  rownames(global.mat) <- as.character(SS)
  colnames(global.mat) <- c("Predicted_global_diversity")  
  
  # Criteria 2: Accuracy to observed dataset diversity
  AccuracyToObserved.mat <- matrix(nrow=length(SS), ncol=1)
  rownames(AccuracyToObserved.mat) <- as.character(SS)
  colnames(AccuracyToObserved.mat) <- c("Accuracy_to_observed_dataset_diversity")
  
  # Criteria 3: Similarity of curves (local and global)
  sim.local.mat <- matrix(rep(0, length(SS)*length(SS)), nrow=length(SS), ncol=length(SS))
  rownames(sim.local.mat) <- as.character(SS)
  colnames(sim.local.mat) <- as.character(SS)
  
  sim.global.mat <- matrix(rep(0, length(SS)*length(SS)), nrow=length(SS), ncol=length(SS))
  rownames(sim.global.mat) <- as.character(SS)
  colnames(sim.global.mat) <- as.character(SS)
  
  # Criteria 3: Mean squared distance between curves (local and global)
  msr.curve.mat <- matrix(nrow=2, ncol=1)
  rownames(msr.curve.mat) <- c("local", "global")
  colnames(msr.curve.mat) <- c("MSR_between_subsample_curves")
  
  # Criteria 3: Distance from reference curve (local and global)
  sim.local.ref <- matrix(nrow=length(SS)-1, ncol=1)
  rownames(sim.local.ref) <- as.character(SS[2:length(SS)])
  colnames(sim.local.ref) <- c("Distance_from_local_reference_curve")
  
  sim.global.ref <- matrix(nrow=length(SS)-1, ncol=1)
  rownames(sim.global.ref) <- as.character(SS[2:length(SS)])
  colnames(sim.global.ref) <- c("Distance_from_global_reference_curve")
  
  # Criteria 4: Plausibility (Monotonicity and decreasing 2nd derivate)
  monotonic.local.mat <- matrix(nrow=length(SS), ncol=1)
  monotonic.global.mat <- matrix(nrow=length(SS), ncol=1)
  slowing.local.mat <- matrix(nrow=length(SS), ncol=1)
  slowing.global.mat <- matrix(nrow=length(SS), ncol=1)
  plausibility.mat <- matrix(nrow=length(SS), ncol=3)
  rownames(monotonic.local.mat) <- as.character(SS)
  rownames(monotonic.global.mat) <- as.character(SS)
  rownames(slowing.local.mat) <- as.character(SS)
  rownames(slowing.global.mat) <- as.character(SS)
  rownames(plausibility.mat) <- as.character(SS)
  colnames(monotonic.local.mat) <- c("Is monotonic (local)")
  colnames(monotonic.global.mat) <- c("Is monotonic (global)")
  colnames(slowing.local.mat) <- c("Has_decreasing_2nd_derivative (local)")
  colnames(slowing.global.mat) <- c("Has_decreasing_2nd_derivative (global)")
  colnames(plausibility.mat) <- c("Is plausible", "Monotonically increasing", "Slowing")
  
  
  ###### Loop through different samples in SS #######
  for (b in 1:length(SS)) {
    cat("Performing fitting routine for sample ", b, "\n")
    if (b>1) { # whether to use meta.list or previous fitting values
      rbind(init.param, fit$par)
    }


    fit <- try(fitsample(model, init.param, lower.bd, upper.bd, dSS[[b]], numit, varleft), silent=TRUE)
    
    if(class(fit) == "try-error") {next} 

    param.mat[b,] <- fit$par # Assign parameter values
    ssr.mat[b,] <- fit$ssr # Assign sum of squared residuals
    ms.mat[b,] <- fit$ms # Assign mean squared residuals
    discrep.mat[b,] <- model.cost.abs(model, fit$par, dSS[[b]]) # Assign discrepancy mean absolute residuals (in percentages)
    local.mat[b,] <- model(length(v1), fit$par) # Assign local diversity prediction
    global.mat[b,] <- model(tot.pop, fit$par) # Assign global diversity prediction
    AccuracyToObserved.mat[b,] <- (model(length(v1), fit$par) - length(unique(v1)))/length(unique(v1)) # Assign % accuracy of local diversity predition
    
    # Assign plausibility assessments
    monotonic.local.mat[b,] <- is.mono(model, fit$par, min.samp=minsampsize, max.samp=length(v1), tot.length=length(v1)-minsampsize+1) 
    monotonic.global.mat[b,] <- is.mono(model, fit$par, min.samp=minsampsize, max.samp=tot.pop, tot.length=min(tot.pop-minsampsize+1,2000))
    slowing.local.mat[b,] <- is.slowing(model, fit$par, min.samp=minsampsize, max.samp=length(v1), tot.length=length(v1)-minsampsize+1) 
    slowing.global.mat[b,] <- is.slowing(model, fit$par, min.samp=minsampsize, max.samp=tot.pop, tot.length=min(tot.pop-minsampsize+1,2000))
    plausibility.mat[b,2] <- ((monotonic.local.mat[b,])*(monotonic.global.mat[b,]))==1
    plausibility.mat[b,3] <- (slowing.global.mat[b,])==1
    plausibility.mat[b,1] <- (plausibility.mat[b,2]*plausibility.mat[b,3])==1
    
  }
  
  # Assign similarity measures (local and global)
  norm.fac.local <- try(ref.area(model, param.mat[1,], minsampsize, length(v1)), silent=TRUE)
  norm.fac.global <- try(ref.area(model, param.mat[1,], minsampsize, tot.pop), silent=TRUE)

  
  for (b in 1:length(SS)) {
    for (c in (b:length(SS))) {
      if (b==c) {
        sim.local.mat[b,c] <- 0
        sim.global.mat[b,c] <- 0
      } else {
        local.tmp.mat <- try(simm(model, param.mat[b,], param.mat[c,], lowerlimit=minsampsize,
                                   upperlimit=length(v1), tot.length=length(v1)-minsampsize+1), silent=TRUE)
        global.tmp.mat <- try(simm(model, param.mat[b,], param.mat[c,], lowerlimit=minsampsize,
                                    upperlimit=tot.pop, tot.length=tot.pop-minsampsize+1), silent=TRUE)
        if((class(local.tmp.mat) != "try-error") & (class(global.tmp.mat) != "try-error")) {
          sim.local.mat[b,c] <- local.tmp.mat
          sim.global.mat[b,c] <- global.tmp.mat
        }
      }
    }
  }

  sim.local.mat <- sim.local.mat + t(sim.local.mat)
  sim.global.mat <- sim.global.mat + t(sim.global.mat)
  sim.local.mat.tmp <- try(sim.local.mat/norm.fac.local, silent=TRUE)
  sim.global.mat.tmp <- try(sim.global.mat/norm.fac.global, silent=TRUE)
  if ((class(sim.local.mat.tmp) != "try-error") & (class(sim.global.mat.tmp) != "try-error")) {
    sim.local.mat <- sim.local.mat.tmp
    sim.global.mat <- sim.global.mat.tmp
  }
  
  for (b in 2:length(SS)) {
    sim.local.ref[b-1,] <- sim.local.mat[b,1]
    sim.global.ref[b-1,] <- sim.global.mat[b,1]
  }
  
  # Assign msr between sample curves
  msr.curve.mat[1,1] <- msr.curves(sim.local.mat)
  msr.curve.mat[2,1] <- msr.curves(sim.global.mat)
  
  list(param=param.mat, ssr=ssr.mat, ms=ms.mat, discrep=discrep.mat, local=local.mat, global=global.mat, AccuracyToObserved=AccuracyToObserved.mat,
       subsamplesizes=SS, datapoints=dSS, modelname=names(model.list), numparam=ncol(param.mat), sampvar=msr.curve.mat,
       mono.local=monotonic.local.mat, mono.global=monotonic.global.mat, slowing.local=slowing.local.mat,
       slowing.global=slowing.global.mat, plausibility=plausibility.mat, dist.local=sim.local.mat, dist.global=sim.global.mat,
       local.ref.dist=sim.local.ref, global.ref.dist=sim.global.ref, popsize=tot.pop, themodel=model)
}


#################### B3a. Generic fit single model function ###################
fitsinglemod <- function(model.list, init.param, param.range, main.samp, tot.pop=(100*(divsamplenum(main.samp,2)[1])),
                         numit=10^5, varleft=1e-8, data.default=TRUE, subsizes=6, dssamps=list(),
                         nrf=1, minrarefac=1, NResamples=1000, minplaus=10, fitloops=2) {
  UseMethod("fitsinglemod")
}


################### B3b. Default fit single model function ####################
fitsinglemod.default <- function(model.list, init.param, param.range, main.samp, tot.pop=(100*(divsamplenum(main.samp,2)[1])),
                                 numit=10^5, varleft=1e-8, data.default=TRUE, subsizes=6, dssamps=list(),
                                 nrf=1, minrarefac=1, NResamples=1000, minplaus=10, fitloops=2) { 
  
  v1 <- format.input(main.samp) # Convert samp to a vector of length=main sample length and with clone ids
  if (data.default) {
    cat("Creating subsamples and rarefaction data", "\n")
    SS <- divsamplenum(main.samp, subsizes)
    
    # Create divsubsample object for each sample
    samp <- v1
    dSS <- list()
    
    main.dss <- divsubsamples(mainsamp=samp, nrf=nrf, minrarefac=minrarefac, maxrarefac=SS[1])

    dSS[[1]] <- main.dss
    for (b in (2:length(SS))) {
      #samp <- sample(samp, size=SS[b], replace=FALSE)
      max.ss <- SS[b]
      dss_temp <- divsubsamples_nested(main.dss, max.ss)
      dSS[[b]] <- dss_temp
    }
  } else {
    SS <- subsizes # largest element of samp.vec must be < length(v1)
    # Create divsubsample object for each sample
    dSS <- dssamps
  }
  
  # Fit the samples
  cat("Fitting loop 1", "\n")
  fitsingle <- fitallSubs(SS, model.list, init.param, param.range, dSS, numit, varleft, tot.pop, v1, minplaus)
  
  if (fitloops > 1) {
    for (c in (2:fitloops)) {
      cat("Fitting loop", c, "\n")
      rbind(init.param, fitsingle$param)
      fitsingle <- fitallSubs(SS, model.list, init.param, param.range, dSS, numit, varleft, tot.pop, v1, minplaus)
    }
  }
  fitsingle$call <- match.call()
  class(fitsingle) <- "fitsingleMod"
  fitsingle
}


################## B3c. Print method -  fit of single model ###################
print.fitsingleMod <- function(x, ...) {
  cat("Fitting results for model:", x$modelname, "\n")
  cat("Model parameters ($param):\n")
  print(x$param)
  
  cat("\nSubsample sizes ($subsamplesizes):\n")
  print(x$subsamplesizes)
  
  #cat("\nSum-of-squares residuals ($ssr):\n")
  #print(x$ssr)
  #cat("\nMean-squared residuals ($ms):\n")
  #print(x$ms)
  #cat("\nDiscrepancy mean absolute residual percentages ($discrep):\n")
  #print(x$discrep)
  
  cat("\nMeasures for discrepancies between rarefaction data and model fits ($discrep, $ssr, $ms):\n")
  print(cbind(x$discrep, x$ssr, x$ms))
  
  #cat("\nLocal diversity predictions ($local):\n")
  #print(x$local)
  #cat("\nGlobal diversity predictions ($global):\n")
  #print(x$global)
  #cat("\nMain sample observation accuracy ($AccuracyToObserved):\n")
  #print(x$AccuracyToObserved)
  
  cat("\nDiversity predictions and the local prediction accuracy ($AccuracyToObserved, $local, $global):\n")
  print(cbind(x$AccuracyToObserved, x$local, x$global))
  
  #cat("\nDistances between sample local fits ($dist.local):\n")
  #print(x$dist.local)
  #cat("\nDistances between sample global fits ($dist.global):\n")
  #print(x$dist.global)
  
  cat("\nSimilarity measure: normalised area between curve and largest-sample curve ($local.ref.dist):\n")
  print(x$local.ref.dist)
  
  #cat("\nMean-squared residuals of sample curves ($sampvar):\n")
  #print(x$sampvar)
  
  cat("\nPlausibility measures: monotonicity and decreasing 2nd derivative ($plausibility):\n")
  print(x$plausibility)

  #cat("\nSubsample and rarefraction datapoints ($datapoints):\n")
  #print(x$datapoints)
}

################## B3d. Summary method: fit of single model ###################
summary.fitsingleMod <- function(object, ...) {
  tot.subsamp <- length(object$datapoints)
  num.rf <- length((object$datapoints[[1]])$RarefacXAxis)
  fitted.model <- object$modelname
  n.param <- object$numparam
  
  TAB <- cbind(No.of.subsamples=tot.subsamp,
               No.of.rarefaction.points=num.rf,
               Fitted.model.name=fitted.model,
               No.of.model.parameters=n.param)
  fit.sum <- list(call=object$call,
                  rsum=TAB)
  class(fit.sum) <- "summary.fitsingleMod"
  fit.sum
}

############## B3e. Print summary method -  fit of single model ###############
print.summary.fitsingleMod <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nSingle model fit setup summary:\n")
  print(x$rsum)
}


################## B3f. Plot method -  fit of single model ####################
plot.fitsingleMod <- function(x, range="local", ...) {
  points.colours <- c("red", "green", "brown", "orange", "blue", "yellow")
  
  # Define plot limits
  ylim.list <- c()
  for (i in 1:length(x$subsamplesizes)) {
    # Define fitted curve
    func1 <- function(z) {
      x$themodel(z, x$param[i,])
    }
    
    if (range=="local") {
      xvalue <- x$datapoints[[i]]$RarefacXAxis[length(x$datapoints[[i]]$RarefacXAxis)] 
    } else {
      xvalue <- x$popsize
    }
    yvalue <- x$datapoints[[i]]$RarefacYAxis[length(x$datapoints[[i]]$RarefacYAxis)]
    ypred <- func1(xvalue)
    max.y <- max(yvalue, ypred)
    ylim.list[i] <- max.y
  }
  ylim.choose <- max(ylim.list)
  
  # Plot main sample diversity
  Sample_number <- x$datapoints[[1]]$RarefacXAxis[length(x$datapoints[[1]]$RarefacXAxis)]
  Species_count <- x$datapoints[[1]]$RarefacYAxis[length(x$datapoints[[1]]$RarefacXAxis)]
  
  if (range=="local") {
    xupper <- Sample_number*1.1
    yupper <- ylim.choose*1.1
  } else if (range=="global") {
    xupper <- x$popsize*1.1
    yupper <- ylim.choose*1.1
  }
  plot(x=Sample_number, y=Species_count, pch=20, col="purple", cex=3, xlim=c(0,xupper),
       ylim=c(0, yupper), ...)
  
  # Plot curves and points
  for (j in 1:length(x$subsamplesizes)) {
    # Plot the points
    points(x=x$datapoints[[j]]$RarefacXAxis,
           y=x$datapoints[[j]]$RarefacYAxis,
           pch=4, col=points.colours[j])
    
    # Define the fitted function
    func1 <- function(z) {
      x$themodel(z, x$param[j,])
    }
    
    # Plot the fitted function
    if (range=="local") {
      xval <- seq(0, Sample_number, 0.1)
      yval <- func1(xval)
    } else if (range=="global") {
      xval <- seq(0, x$popsize, 1)
      yval <- func1(xval)
    }
    lines(x=xval, y=yval, lty=1, col=points.colours[j], lwd=2)
  }
  
  # Purpledot
  points(x=Sample_number, y=Species_count, pch=20, col="purple", cex=3)
}

####################### C. SCORING AND MODEL COMPARISON #######################
############################ C1a. Whole number check ##########################
# Input: Single number
# Output: TRUE/FALSE
IsWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  return(abs(x - round(x)) < tol)  
}

########## C1b. Binning function - Single number, absolute divisor ############
# Input: Number, divisor
# Output: Number rounded UP to multiple of divisor
BinFunkSingleNumber <- function(x,d) {
  if (is.na(x))  {
    return(x)
  } else if ( x == 0 )	{	## need this line as otherwise 0 value gets 0 score, want minimum score to be 1, regardless of divisor.
    return(1)
  }	else {
    X = abs(x)
    res = X%%d
    if(IsWholeNumber(X/d))	{	
      return((	((X-res)/d) ) )				
    }	else {
      return( (	((X-res)/d) + 1	) 	)
    }
  }
}

######### C1c. Binning function - Vector of numbers, absolute divisor #########
# Input: Vector of numbers, divisor
# Output: Vector of number rounded UP to multiple of divisor
BinFunk <- function(xx,dd) {  	# this rounds each entry of xx DOWN to the last multiple of dd. 
  return(sapply(xx , BinFunkSingleNumber, d = dd))
}

######### C1d. Binning function - Single number, Percentage precision #########
# Input: Number, precison
# Output: Number rounded to specified precision
Bfunprop = function(a, percent) {
  if (is.na(a)) {
    return(a)
  } else {
    if (a >= 0) {
      return( floor((1/percent) * a)/(1/percent)	)
    }
    if (a < 0) {
      return( ceiling((1/percent) * a)/(1/percent)	)
    }  
  }
}

####### C1e. Binning function - Vector of numbers, Percentage precision #######
# Input: Vector of numbers, precison
# Output: Vector of numbers rounded to specified precision
BBprop=function(aa,percentwecareabout) {  
  if(percentwecareabout == 0)	{
    return(aa)
  } else	{
    return(sapply(aa , Bfunprop, percent = percentwecareabout))
  }
}

####### C1f. Geometric mean function #######
# Input: Top estimates
# Output: Geometric mean
geo.mean <- function(x) {
  if (any(is.na(x))) {warning('global diversity values are unreliable due to one or more poor model fits amongst top models')}
  return(exp(mean(log(x), na.rm=TRUE)))
}

##################### C2a. Scoring function - Single model ####################
# Input: FitsingleMod, Sample aggregation method ("Equal", "Proportional", "InvProportional"), 
#         Criteria rounding values c(fit, accuracy, local similarity), plausibility penalty score
# Output: Single model scores - Five criteria - 1. Fit, 2. Accuracy, 3. Local similarity, 
#         4. Monotonicity (local-global combined), 5. Plausible-global
single.mod.score <- function(fsm, precision.lv = c(0.0001, 0.005, 0.005), plaus.pen=500) {
  nS <- length(fsm$subsamplesizes) # No. of samples
  samp.wts <- rep(1/nS, nS)

  
  ## Un-aggregated scores
  fit.scores <- BinFunk(fsm$discrep, precision.lv[1])
  accuracy.scores <- BinFunk(fsm$AccuracyToObserved, precision.lv[2])
  similarity.scores <- BinFunk(fsm$local.ref.dist, precision.lv[3])
  plausible.scores <- rep(500, nS)
  
  # Monotonic and plausibility scores penalties
  for (i in 1:nS) {
#     if (is.na(fsm$mono.local[i,]) || is.na(fsm$mono.global[i,])) {
#       monotonic.scores[i] <- NA
#     } else if (fsm$mono.local[i,]==FALSE || fsm$mono.global[i,]==FALSE) {
#       monotonic.scores[i] <- plaus.pen
#     }
#     if (is.na(fsm$slowing.global[i,])) {
#       plausible.scores[i] <- NA
#     } else if (fsm$slowing.global[i,] == FALSE) {
#       plausible.scores[i] <- plaus.pen
#     }
#     if (plausible.scores[i] == plaus.pen && monotonic.scores[i] == plaus.pen) {
#       monotonic.scores[i] <- plaus.pen/2
#       plausible.scores[i] <- plaus.pen/2
#     }
    
    if (is.na(fsm$plausibility[i,1]) || fsm$plausibility[i,1]==FALSE) {
      plausible.scores <- plaus.pen
    }
  }
  
  ## Aggregated scores
  fit.agg <- sum(fit.scores*samp.wts)
  accuracy.agg <- sum(accuracy.scores*samp.wts)
  similarity.agg <- sum(similarity.scores*samp.wts[2:length(samp.wts)])/sum(samp.wts[2:length(samp.wts)])
  #monotonic.agg <- sum(plausible.scores*samp.wts)
  plausible.agg <- sum(plausible.scores*samp.wts)
  list(fit=fit.agg, accuracy=accuracy.agg, similarity=similarity.agg, plausibility=plausible.agg,
       binsize=precision.lv, plausibility.penalty=plaus.pen,
       modname=fsm$modelname)
}

################### C3a. Generic score single model function ##################
scoresinglemod <- function(fsm, precision.lv = c(0.0001, 0.005, 0.005), plaus.pen=500) {
  UseMethod("scoresinglemod")
}

################### C3b. Default score single model function ##################
scoresinglemod.default <- function(fsm, precision.lv = c(0.0001, 0.005, 0.005), plaus.pen=500) {
  scoresingle <- single.mod.score(fsm, precision.lv, plaus.pen)
  scoresingle$call <- match.call()
  class(scoresingle) <- "scoresingleMod"
  scoresingle
}

################# C3c. Print method -  score of single model ##################
print.scoresingleMod <- function(x, ...) {
  cat("Scoring results ($[colnames]):\n")
  score.output <- matrix(c(x$fit, x$accuracy, x$similarity, x$plausibility), nrow=1, ncol=4)
  rownames(score.output) <- as.character(x$modname)
  colnames(score.output) <- c("discrepancy", "accuracy", "similarity", "plausibility")
  print(score.output)
}

################# C3d. Summary method: score of single model ##################
summary.scoresingleMod <- function(object, ...) {
  fitted.model <- object$modname
  binsize <- object$binsize
  plauspen <- object$plausibility.penalty
  
  TAB <- cbind(Model.name=fitted.model,
               precision.fit=binsize[1],
               precision.accuracy=binsize[2],
               precision.similarity=binsize[3],
               Plausibility.penalty=plauspen)
  score <- list(call=object$call,
              rsum=TAB)
  class(score) <- "summary.scoresingleMod"
  score
}

############# C3e. Print summary method -  score of single model ##############
print.summary.scoresingleMod <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nSingle model scoring setup summary:\n")
  print(x$rsum)
}


################# C4a. Combine attributes of a single model ###################
# Input: scoresingleMod object, Criteria weights (zero=exclude from scoring)
# Output: Single model score
combine.criteria <- function(ssm, crit.wts=c(1.0, 1.0, 1.0, 1.0)) {
  vect <- c(ssm$fit, ssm$accuracy, ssm$similarity, ssm$plausibility)
  vect <- sum((vect*crit.wts)/sum(crit.wts))
  return(vect)
}

################# C4b. Fit, score and compare multiple models #################
# Input: List of models, List of initial parameters, List of Lower bounds, List of upper bounds,
#       main.samp, tot.pop, numit, varleft, subsizes, dssamps, nrf, minrarefac, NResamples, minplaus,
#       precision.lv, plaus.pen, crit.wts, fitloops, numpred
# Output: List of i) Model scores ii) List of dss iii) List of ssm iv) predictions x3 vii) original model list
multiple.scoring <- function(models, init.params, param.ranges, main.samp, tot.pop, numit, varleft,
                             subsizes, dssamps, nrf, minrarefac, NResamples, minplaus, precision.lv, plaus.pen, crit.wts,
                             fitloops, numpred) {
  num.mod <- length(models)
  fmm <- list()
  ssm <- matrix(rep(NA, num.mod*4), nrow=num.mod, ncol=4)
  colnames(ssm) <- c("fit", "accuracy", "similarity", "plausibility")
  mod.score <- matrix(rep(NA, num.mod), nrow=num.mod, ncol=1)
  mod.rownames <- as.character(rep(NA, num.mod))
  
  # Create the samples and divsubsample object
  if (length(dssamps)==0) {
    cat("Creating subsamples and rarefaction data", "\n")
    mas.dss <- list() # Create 3 master samples
    mas.SS <- divsamplenum(main.samp, subsizes) # Create master vector of 3 samples sizes
    samp <- format.input(main.samp) # Format main sample
    
    dss.main <- divsubsamples(mainsamp=samp, nrf=nrf, minrarefac=minrarefac, maxrarefac=mas.SS[1], NResamples=NResamples)
    mas.dss[[1]] <- dss.main
    for (b in (2:length(mas.SS))) { # Create master list of divsubsamples
      #samp <- sample(samp, size=mas.SS[b], replace=FALSE)
      max.ss <- mas.SS[b]
      dss <- divsubsamples_nested(dss.main, max.ss)
      mas.dss[[b]] <- dss  
    }
  } else {
    cat("Loading predefined subsamples", "\n")
    mas.SS <- divsamplenum(main.samp, subsizes)
    samp <- format.input(main.samp) # Format main sample

    mas.dss <- dssamps
  }

  timeremain = "tbd"
  ptm <- proc.time()[3]
  # Loop throught the different models
  for (i in (1:num.mod)) {
    # Fit model
    cat("Fitting model", i, "of", length(models),"(Est. time remaining:", timeremain, "mins)", "\n")
    fsm.temp <- fitsinglemod(model.list=models[i], init.param=init.params[[i]], param.range=param.ranges[[i]],
                             main.samp, data.default=FALSE, tot.pop= tot.pop, numit=numit, varleft=varleft,
                             subsizes=mas.SS, dssamps=mas.dss, nrf=nrf, minrarefac=minrarefac, NResamples=NResamples,
                             minplaus=minplaus, fitloops=fitloops)
    
    # Score model criteria
    cat("Scoring model", i, "\n")
    ssm.temp <- try(scoresinglemod(fsm=fsm.temp, precision.lv = precision.lv, plaus.pen=plaus.pen), silent=TRUE)
    if(class(ssm.temp) == "try-error") {next}  
    
    # Aggregate criteria score
    cat("Aggregating scoring for model", i, "\n")
    mod.score.temp <- combine.criteria(ssm=ssm.temp, crit.wts=crit.wts)
    
    fmm[[fsm.temp$modelname]] <- fsm.temp
    ssm[i,1] <- ssm.temp$fit
    ssm[i,2] <- ssm.temp$accuracy
    ssm[i,3] <- ssm.temp$similarity
    ssm[i,4] <- ssm.temp$plausibility
    mod.score[i,] <- mod.score.temp
    mod.rownames[i] <- fsm.temp$modelname
    
    # Calculate time
    timeremain <- round(((proc.time()[3] - ptm)*(num.mod - i)/i)/60)
  }
  rownames(mod.score) <- mod.rownames
  colnames(mod.score) <- "Combined score"
  rownames(ssm) <- mod.rownames
  
  # calculate estimate
  num.models <- length(fmm)
  m <- min(numpred, num.models)
  topx_scores <- sort(unique(mod.score))[1:m]
  lenx <- length(topx_scores)
  topx_index <- which(mod.score %in% topx_scores)
  predx.vector <- c()
  for (i in 1:lenx) {
    tmp <- ((fmm[[topx_index[i]]])$global)[1,]
    predx.vector <- c(predx.vector, tmp)
  }
  predx <- geo.mean(predx.vector)
  predx_high <- max(predx.vector)
  predx_low <- min(predx.vector)
  
  list(model.scores=mod.score, fmm=fmm, ssm=ssm, estimate=predx, upper_estimate=predx_high, lower_estimate=predx_low, models=models, m=m)
}


############### C5a. Generic fit/score multiple model function ################
DiveMaster <- function(models, init.params, param.ranges, main.samp, tot.pop=(100*(divsamplenum(main.samp,2)[1])), numit=10^5, varleft=1e-8,
                       subsizes=6, dssamps=list(), nrf=1, minrarefac=1, NResamples=1000, minplaus=10, 
                       precision.lv=c(0.0001, 0.005, 0.005), plaus.pen=500,
                       crit.wts=c(1.0, 1.0, 1.0, 1.0), fitloops=2, numpred=5) {
  UseMethod("DiveMaster")
}

################ C5b. Default fit/score single model function #################
DiveMaster.default <- function(models, init.params, param.ranges, main.samp, tot.pop=(100*(divsamplenum(main.samp,2)[1])), numit=10^5, varleft=1e-8,
                              subsizes=6, dssamps=list(), nrf=1, minrarefac=1, NResamples=1000, minplaus=10, 
                              precision.lv=c(0.0001, 0.005, 0.005), plaus.pen=500,
                              crit.wts=c(1.0, 1.0, 1.0, 1.0), fitloops=2, numpred=5) {
  mainfit <- multiple.scoring(models, init.params, param.ranges, main.samp, tot.pop, numit, varleft,
                              subsizes, dssamps, nrf, minrarefac, NResamples, minplaus, precision.lv, plaus.pen, crit.wts, fitloops, numpred)
  mainfit$call <- match.call()
  class(mainfit) <- "DiversityMaster"
  mainfit
}

################ C5c. Print method - score of multiple models #################
print.DiversityMaster <- function(x, ...) {
  cat("Model score ($model.scores):\n")
  print(x$model.scores)
  cat("\nPopulation diversity estimate ($estimate):\n")
  print(x$estimate)
}

############# C5d. Summary method: fit/scores of multiple models ##############
summary.DiversityMaster <- function(object, ...) {
  num.models <- length(object$fmm)
  
  best.index <- which.min(object$model.scores)
  best.model <- (object$fmm[[best.index]])$modelname
  div.pred.from.top.model <- ((object$fmm[[best.index]])$global)[1,]
  
#   m <- min(5, num.models)
#   top5_scores <- sort(unique(object$model.scores))[1:m]
#   len5 <- length(top5_scores)
#   top5_index <- which(object$model.scores %in% top5_scores)
#   pred5.vector <- c()
#   for (i in 1:len5) {
#     tmp <- ((object$fmm[[top5_index[i]]])$global)[1,]
#     pred5.vector <- c(pred5.vector, tmp)
#   }
#   pred5 <- geo.mean(pred5.vector)
  
  TAB <- data.frame(Number_of_models_fitted=as.numeric(num.models),
               Best_scoring_model=best.model,
               Estimated_global_diversity=as.numeric(object$estimate),
               Diversity_upperbound=as.numeric(object$upper_estimate),
               Diversity_lowerbound=as.numeric(object$lower_estimate))

  res <- list(call=object$call,
              rsum=TAB)
  class(res) <- "summary.DiversityMaster"
  res
}

######### C5e. Print summary method -  fits/scores of multiple models #########
print.summary.DiversityMaster <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nMultiple Model fitting and scoring summary:\n")
  print(x$rsum)
}

######### C6. Combining several DiveMaster objects into a single DiveMaster object #########
# Input: List of DiveMaster objects
# Output: Single combined DiveMaster object
comb.dm <- function(dmlist) {
  l <- length(dmlist) # number of DiveMaster objects
  out <- dmlist[[1]]
  if (l>1) {
    for (i in 2:l){
      out$model.scores <- rbind(out$model.scores, dmlist[[i]]$model.scores)
      out$ssm <- rbind(out$ssm, dmlist[[i]]$ssm)
      out$fmm <- c(out$fmm, dmlist[[i]]$fmm)
      out$models <- c(out$models, dmlist[[i]]$models)
      out$m <- min(out$m, dmlist[[i]]$m)
    }
    
    # Recompute estimate
    num.models <- length(out$fmm)
    m <- min(out$m, num.models)
    topx_scores <- sort(unique(out$model.scores))[1:m]
    lenx <- length(topx_scores)
    topx_index <- which(out$model.scores %in% topx_scores)
    predx.vector <- c()
    for (j in 1:lenx) {
      tmp <- ((out$fmm[[topx_index[j]]])$global)[1,]
      predx.vector <- c(predx.vector, tmp)
    }
    predx <- geo.mean(predx.vector)
    predx_high <- max(predx.vector)
    predx_low <- min(predx.vector)
    
    out$estimate <- predx
    out$upper_estimate <- predx_high
    out$lower_estimate <- predx_low
  }
  class(out) <- "DiversityMaster"
  
  return(out)
}


######### C7. Give a diversity estimate for a different population size based on top-x score models #########
# Input: DiveMaster object, population size, Number of top scores models to include in estimate
# Output: Single population estimate
popdiversity = function (dm, popsize, TopX = NULL) {
  if (is.null(TopX))	{	
    m <- dm$m	
  } else {
    if (TopX <= length(dm$fmm)) {	
      m <- TopX
    } else 	{ 	
      m <- min(5, length(dm$fmm))		
      warning("Number of models specified is greater then number of models fitted")
    }	 
  }
  
  topx_scores <- sort(unique(dm$model.scores))[1:m]
  lenx <- length(topx_scores)
  topx_index <- which(dm$model.scores %in% topx_scores)
  predx.vector <- c()
  for (i in 1:lenx) {
    tmp <- dm$models[[topx_index[i]]](popsize, dm$fmm[[topx_index[i]]]$param[1, ])
    predx.vector <- c(predx.vector, tmp)
  }
  predx <- geo.mean(predx.vector)
  predx_high <- max(predx.vector)
  predx_low <- min(predx.vector)
  return(list(estimate = predx, upper_estimate = predx_high, 
              lower_estimate = predx_low))
}






