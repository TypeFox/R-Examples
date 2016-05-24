### This file contains two type of methods: MLE and posterior mean for
### estimating an appropriate initial value of expected gene expression phi
### for one gene given sequence information/counts.

### - fitlist is a list beta estimation from fitMultinom(). amino acid -> beta.
### - reu13.list.g is a gene list. amino acid -> codon -> position.
### - y.g is a gene list. amino acid -> codon count.
### - n.g is a gene list. amino acid -> total codon count.

### Get the specific function according to the options.
get.my.estimatePhiOne <- function(init.Phi, model){
  if(!any(init.Phi[1] %in% .CF.CT$init.Phi)){
    stop("initial method is not found.")
  }
  if(!any(model[1] %in% .CF.CT$model)){
    stop("model is not found.")
  }
  ret <- eval(parse(text = paste("my.estimatePhiOne.", init.Phi[1],
                                 ".", model[1], sep = "")))
  assign("my.estimatePhiOne", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.estimatePhiOne().


### For mle methods.
### For ROC + NSEf model.
my.estimatePhiOne.MLE.rocnsef <- function(fitlist, reu13.list.g, y.g, n.g,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  tmp <- optim(E.Phi, my.objectivePhiOne.nlogL.rocnsef, gr = NULL,
               fitlist, reu13.list.g, y.g, n.g, method = .CF.OP$optim.method[1],
               lower = lower.optim, upper = upper.optim, control = control)
  ret <- tmp$par
  ret
} # End of my.estimatePhiOne.MLE.rocnsef().

### For ROC model.
my.estimatePhiOne.MLE.roc <- function(fitlist, reu13.list.g, y.g, n.g,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  tmp <- optim(E.Phi, my.objectivePhiOne.nlogL.roc, gr = NULL,
               fitlist, reu13.list.g, y.g, n.g, method = .CF.OP$optim.method[1],
               lower = lower.optim, upper = upper.optim, control = control)
  ret <- tmp$par
  ret
} # End of my.estimatePhiOne.MLE.roc().

### For NSEf model.
my.estimatePhiOne.MLE.nsef <- function(fitlist, reu13.list.g, y.g, n.g,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  tmp <- optim(E.Phi, my.objectivePhiOne.nlogL.nsef, gr = NULL,
	       fitlist, reu13.list.g, y.g, n.g, method = .CF.OP$optim.method[1],
	       lower = lower.optim, upper = upper.optim, control = control)
  ret <- tmp$par
  ret
} # End of my.estimatePhiOne.MLE.nsef().


### For posterior mean methods.
### For ROC + NSEf model.
my.estimatePhiOne.PM.rocnsef <- function(fitlist, reu13.list.g, y.g, n.g,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  ### Check lower first.
  lphilv.g <- -my.objectivePhiOne.nlogphiL.rocnsef(lower.optim, fitlist,
                 reu13.list.g, y.g, n.g)

  ### set default values.
  add.lphilv.g <- 0.0
  mlphilv.g <- 0.0
  add.llv.g <- 0.0
  mllv.g <- 0.0

  ### Find the maximum value of log phi * L and adjust accordingly.
  if(lphilv.g <= .CF.OP$stable.min.exp || lphilv.g >= .CF.OP$stable.max.exp){
    tmp <- optim(E.Phi, my.objectivePhiOne.nlogphiL.rocnsef, gr = NULL,
                 fitlist, reu13.list.g, y.g, n.g,
                 method = .CF.OP$optim.method[1],
                 lower = lower.optim, upper = upper.optim, control = control)
    mlphilv.g <- -tmp$value
    if(mlphilv.g <= .CF.OP$stable.min.exp ||
       mlphilv.g >= .CF.OP$stable.max.exp){
      add.lphilv.g <- .CF.OP$stable.max.exp - mlphilv.g
    }

    ### Find the maximum value of log L and adjust accordingly.
    if(mlphilv.g - log(tmp$par) <= .CF.OP$stable.min.exp ||
       mlphilv.g - log(tmp$par) >= .CF.OP$stable.max.exp){
      tmp <- optim(tmp$par, my.objectivePhiOne.nlogL.rocnsef, gr = NULL,
                   fitlist, reu13.list.g, y.g, n.g,
                   method = .CF.OP$optim.method[1],
                   lower = lower.optim, upper = upper.optim, control = control)
      mllv.g <- -tmp$value
      if(mllv.g <= .CF.OP$stable.min.exp || mllv.g >= .CF.OP$stable.max.exp){
        add.llv.g <- .CF.OP$stable.max.exp - mllv.g
      }
    }
  }

  ### Integrate scaled phi * L over phi.
  numerator <- my.integrate(Vectorize(my.objectivePhiOne.xLfp.rocnsef, "phi"),
                            lower.integrate, upper.integrate,
                            fitlist, reu13.list.g, y.g, n.g, add.lphilv.g,
                            control = control)

  ### Integrate scaled L over phi.
  denominator <- my.integrate(Vectorize(my.objectivePhiOne.Lfp.rocnsef, "phi"),
                              lower.integrate, upper.integrate,
                              fitlist, reu13.list.g, y.g, n.g, add.llv.g,
                              control = control)

  ### Compute the posterior mean.
  ret <- (numerator$value / denominator$value) * exp(add.llv.g - add.lphilv.g)
  ret
} # End of my.estimatePhiOne.PM.rocnsef().

### For ROC model.
my.estimatePhiOne.PM.roc <- function(fitlist, reu13.list.g, y.g, n.g,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  ### Check lower first.
  lphilv.g <- -my.objectivePhiOne.nlogphiL.roc(lower.optim, fitlist,
                 reu13.list.g, y.g, n.g)

  ### set default values.
  add.lphilv.g <- 0.0
  mlphilv.g <- 0.0
  add.llv.g <- 0.0
  mllv.g <- 0.0

  ### Find the maximum value of log phi * L and adjust accordingly.
  if(lphilv.g <= .CF.OP$stable.min.exp || lphilv.g >= .CF.OP$stable.max.exp){
    tmp <- optim(E.Phi, my.objectivePhiOne.nlogphiL.roc, gr = NULL,
                 fitlist, reu13.list.g, y.g, n.g,
                 method = .CF.OP$optim.method[1],
                 lower = lower.optim, upper = upper.optim, control = control)
    mlphilv.g <- -tmp$value
    if(mlphilv.g <= .CF.OP$stable.min.exp ||
       mlphilv.g >= .CF.OP$stable.max.exp){
      add.lphilv.g <- .CF.OP$stable.max.exp - mlphilv.g
    }

    ### Find the maximum value of log L and adjust accordingly.
    if(mlphilv.g - log(tmp$par) <= .CF.OP$stable.min.exp ||
       mlphilv.g - log(tmp$par) >= .CF.OP$stable.max.exp){
      tmp <- optim(tmp$par, my.objectivePhiOne.nlogL.roc, gr = NULL,
                   fitlist, reu13.list.g, y.g, n.g,
                   method = .CF.OP$optim.method[1],
                   lower = lower.optim, upper = upper.optim, control = control)
      mllv.g <- -tmp$value
      if(mllv.g <= .CF.OP$stable.min.exp || mllv.g >= .CF.OP$stable.max.exp){
        add.llv.g <- .CF.OP$stable.max.exp - mllv.g
      }
    }
  }

  ### integrate phi * L over phi.
  numerator <- my.integrate(Vectorize(my.objectivePhiOne.xLfp.roc, "phi"),
                            lower.integrate, upper.integrate,
                            fitlist, reu13.list.g, y.g, n.g, add.lphilv.g,
                            control = control)

  ### integrate L over phi.
  denominator <- my.integrate(Vectorize(my.objectivePhiOne.Lfp.roc, "phi"),
                              lower.integrate, upper.integrate,
                              fitlist, reu13.list.g, y.g, n.g, add.llv.g,
                              control = control)

  ### Compute the posterior mean.
  ret <- (numerator$value / denominator$value) * exp(add.llv.g - add.lphilv.g)
  ret
} # End of my.estimatePhiOne.PM.roc().

### For NSEf model.
my.estimatePhiOne.PM.nsef <- function(fitlist, reu13.list.g, y.g, n.g,
    E.Phi = .CF.OP$E.Phi, lower.optim = .CF.OP$lower.optim,
    upper.optim = .CF.OP$upper.optim, lower.integrate = .CF.OP$lower.integrate,
    upper.integrate = .CF.OP$upper.integrate, control = list()){
  ### Check lower first.
  lphilv.g <- -my.objectivePhiOne.nlogphiL.nsef(lower.optim, fitlist,
                 reu13.list.g, y.g, n.g)

  ### set default values.
  add.lphilv.g <- 0.0
  mlphilv.g <- 0.0
  add.llv.g <- 0.0
  mllv.g <- 0.0

  ### Find the maximum value of log phi * L and adjust accordingly.
  if(lphilv.g <= .CF.OP$stable.min.exp || lphilv.g >= .CF.OP$stable.max.exp){
    tmp <- optim(E.Phi, my.objectivePhiOne.nlogphiL.nsef, gr = NULL,
                 fitlist, reu13.list.g, y.g, n.g,
                 method = .CF.OP$optim.method[1],
                 lower = lower.optim, upper = upper.optim, control = control)
    mlphilv.g <- -tmp$value
    if(mlphilv.g <= .CF.OP$stable.min.exp ||
       mlphilv.g >= .CF.OP$stable.max.exp){
      add.lphilv.g <- .CF.OP$stable.max.exp - mlphilv.g
    }

    ### Find the maximum value of log L and adjust accordingly.
    if(mlphilv.g - log(tmp$par) <= .CF.OP$stable.min.exp ||
       mlphilv.g - log(tmp$par) >= .CF.OP$stable.max.exp){
      tmp <- optim(tmp$par, my.objectivePhiOne.nlogL.nsef, gr = NULL,
                   fitlist, reu13.list.g, y.g, n.g,
                   method = .CF.OP$optim.method[1],
                   lower = lower.optim, upper = upper.optim, control = control)
      mllv.g <- -tmp$value
      if(mllv.g <= .CF.OP$stable.min.exp || mllv.g >= .CF.OP$stable.max.exp){
        add.llv.g <- .CF.OP$stable.max.exp - mllv.g
      }
    }
  }

  ### integrate phi * L over phi.
  numerator <- my.integrate(Vectorize(my.objectivePhiOne.xLfp.nsef, "phi"),
                            lower.integrate, upper.integrate,
                            fitlist, reu13.list.g, y.g, n.g, add.lphilv.g,
                            control = control)

  ### integrate L over phi.
  denominator <- my.integrate(Vectorize(my.objectivePhiOne.Lfp.nsef, "phi"),
                              lower.integrate, upper.integrate,
                              fitlist, reu13.list.g, y.g, n.g, add.llv.g,
                              control = control)

  ### Compute the posterior mean.
  ret <- (numerator$value / denominator$value) * exp(add.llv.g - add.lphilv.g)
  ret
} # End of my.estimatePhiOne.PM.nsef().


### Adopted from integrate() for more flexible controls.
my.integrate <- function(f, lower, upper, ..., control = list()){
  if(is.null(control$subdivisions)){
    subdivisions <- 100L
  } else{
    subdivisions <- control$subdivisions
  }

  if(is.null(control$rel.tol)){
    rel.tol <- .Machine$double.eps^0.25
  } else{
    rel.tol <- control$rel.tol
  }

  if(is.null(control$abs.tol)){
    abs.tol <- rel.tol
  } else{
    abs.tol <- control$abs.tol
  }

  if(is.null(control$stop.on.error)){
    stop.on.error <- TRUE
  } else {
    stop.on.error <- control$stop.on.error
  }

  keep.xy <-
  if(is.null(control$keep.xy)){
    keep.xy <- FALSE
  } else{
    keep.xy <- control$keep.xy
  }

  if(is.null(control$aux)){
    aux <- NULL
  } else{
    aux <- control$aux
  }

  ret <- integrate(f, lower, upper, ...,
                   subdivisions = subdivisions,
                   rel.tol = rel.tol, abs.tol = abs.tol,
                   stop.on.error = stop.on.error, keep.xy = keep.xy,
                   aux = aux)
  ret
} # End of my.integrate().
