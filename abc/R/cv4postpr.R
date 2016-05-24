######################################################################
#
# cv4postpr.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/abc package
# Contains: cv4postpr, is.cv4postpr, summary.cv4postpr, plot.cv4postpr 
#
######################################################################
## cross-validation for model selection
cv4postpr <- function(index, sumstat, postpr.out = NULL, nval, tols,
                      method, subset = NULL, kernel = "epanechnikov",
                      numnet = 10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit = 500, ...){

  linout <- TRUE
  ## checks:
  ## ######
  if(missing(nval)) stop("'nval' must be supplied.", call.=F)
  if(is.null(postpr.out) && missing(method)) stop("Method must be supplied when 'postpr.out' is NULL.", call.=F)
  if(length(index)!=na.omit(length(index))) stop("'index' contains missing values. Models must be specified for each simulation.", call.=F)

  if(!prod(table(index)>nval)) stop("'nval' has to be smaller or equal to number of simulations for any of the models. Choose a smaller 'nval'.", call.=F)
  
  ## set random seeds
  ## ################
  if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  
  ## define defaults:
  ## #################
  
  if(!is.null(postpr.out)){
    subset <- postpr.out$na.action
    method <- postpr.out$method
    kernel <- "epanechnikov"
  }
  
  ## checks, numbers of stats, params, sims
  if(is.vector(sumstat)){
      numstat <- 1
      sumstat <- matrix(sumstat, ncol=1)
  }
  else numstat <- dim(sumstat)[2]
  numsim <- length(index)
  
  ## names
  if(!is.null(postpr.out)){ # indexnames & statnames from postpr.out
    if(numstat != postpr.out$numstat || numsim != length(postpr.out$na.action)){
      stop("The number of summary statistics, or simulations provided in 'sumstat' are not the same as in 'postpr.out'.", call.=F)
    }
    else if(!prod(unique(index) %in% postpr.out$names$models)){
      stop("Models in 'index' are not the same as in 'postpr.out', or different names are used.", call.=F)
    }
    else if(!prod(colnames(sumstat) %in% postpr.out$names$statistics.names)){
      stop("Summary statistics in 'sumstat' are not the same as in 'postpr.out', or different names are used.", call.=F)
    }
    else{
      mymodels <- postpr.out$names$models
      statnames <- postpr.out$names$statistics.names
    }
  }
  else{ # statnames o/w
    mymodels <- levels(factor(index))
    if(length(colnames(sumstat))){
      statnames <- colnames(sumstat)
    }
    else{
      warning("No statistics names are given, using S1, S2, ...", call.=F, immediate=T)
      statnames <- paste("S", 1:numstat, sep="")
    }
  }

  ## ## order data by models
  ## index <- factor(index)
  ## sumstat <- sumstat[order(index), ]
  ## index <- index[order(index)]
  
  ## indices for the CV sample and check that the sample is not actually an NA
  gwt <- rep(TRUE,length(sumstat[,1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if(is.null(subset)) subset <- rep(TRUE,length(sumstat[,1]))
  gwt <- as.logical(gwt*subset)
  ## CV samples
  cvsamp <- unlist(tapply(c(1:length(index))[gwt], index[gwt], sample, nval))

  ## if tols is a vector have to loop through all values
  tols <- sort(tols)
  allprobs <- list()
  mycall <- list()
  for(mytol in tols){
      res <- matrix(ncol = length(unique(index)), nrow = length(cvsamp))
      for(i in 1:length(cvsamp)){
        ## things to over-write from original call: tolerances, target, index, sumstat
        mysamp <- cvsamp[i]
        mytrue <- index[mysamp]
        mytarget <- sumstat[mysamp,]
        myindex <- index[-mysamp]
        mysumstat <- sumstat[-mysamp,]
        mysubset <- subset[-mysamp]
        subres <- postpr(target = mytarget, index = myindex, sumstat = mysumstat, tol=mytol,subset = mysubset, method = method, kernel = kernel)

        if(subres$method=="rejection") res[i,] <- summary.postpr(subres, print = F, ...)$Prob
        if(subres$method=="mnlogistic") res[i, ] <- summary.postpr(subres, print = F, ...)$mnlogistic$Prob
        if(subres$method=="neuralnet") res[i, ] <- summary.postpr(subres, print = F, ...)$neuralnet$Prob
      }
      colnames(res) <- mymodels
      rownames(res) <- index[cvsamp]
      allprobs[[paste("tol", mytol, sep="")]] <- res
      allnames <- lapply(allprobs, apply, 1, function(xx) mymodels[which(xx==max(xx))])
      
      mycall[[paste("tol", mytol, sep="")]] <-
        call("postpr", target = quote(target), index = quote(index), sumstat = quote(sumstat), tol= mytol,
             subset = quote(subset), method = subres$method, kernel = subres$kernel)
  }
  
  cv4postpr.out <-
    list(calls = mycall, cvsamples = cvsamp, tols = tols, true = index[cvsamp],
         estim = allnames,
         model.probs = allprobs,
         method = method, names = list(models = mymodels, statistics.names = statnames), seed = seed)

  class(cv4postpr.out) <- "cv4postpr"
  invisible(cv4postpr.out)
  
}

is.cv4postpr <- function(x){
    if (inherits(x, "cv4postpr")) TRUE
    else FALSE
}

summary.cv4postpr <- function(object, probs=TRUE, print = TRUE, digits = max(3, getOption("digits")-3), ...){

  if (!inherits(object, "cv4postpr")) 
    stop("Use only with objects of class \"cv4postpr\".", call.=F)
  
  cv4postpr.out <- object
  tols <- cv4postpr.out$tols
  numtols <- length(tols)
  true <- cv4postpr.out$true
  estim <- cv4postpr.out$estim
  method <- cv4postpr.out$method
  model.probs <- cv4postpr.out$model.probs
  nmodels <- length(cv4postpr.out$names$models)
  nval <- length(true)/nmodels

  if(print) cat("Confusion matrix based on ", nval, " samples for each model.\n\n", sep="")
  cm <- lapply(estim, function(x) table(true, x))
  cm <- lapply(cm, function(x) {attributes(dimnames(x))$names <- NULL; x})
  if(print) print(cm); cat("\n")

  if(probs){
      if(print) cat(paste("Mean model posterior probabilities (", method ,")\n\n", sep="")) 
      myprs <- lapply(model.probs, apply, 2, tapply, true, mean)
      if(print) print(lapply(myprs, round, digits=digits))
      out <- list(conf.matrix=cm, probs=myprs)
  }
  else out <- cm
  
  invisible(out)
}

plot.cv4postpr <- function(x, probs=FALSE, file = NULL, postscript = FALSE, onefile = TRUE, ask = !is.null(deviceIsInteractive()), caption = NULL, ...){

  if (!inherits(x, "cv4postpr")) 
    stop("Use only with objects of class \"cv4postpr\".", call.=F)
  
  cv4postpr.out <- x
  tols <- cv4postpr.out$tols
  numtols <- length(tols)
  true <- cv4postpr.out$true
  estim <- cv4postpr.out$estim
  method <- cv4postpr.out$method
  model.probs <- cv4postpr.out$model.probs
  nmodels <- length(cv4postpr.out$names$models)
  nval <- length(true)/nmodels
  
  if(is.null(caption)) caption <- "Confusion matrix"
  
  ## Devices
  save.devAskNewPage <- devAskNewPage()
  if(!is.null(file)){
      file <- substitute(file)
      if(!postscript) pdf(file = paste(file, "pdf", sep="."), onefile=onefile)
      if(postscript) postscript(file = paste(file, "ps", sep="."), onefile=onefile)
  }
  else{
      if (ask && 1 < numtols) {
          devAskNewPage(TRUE)
      }
  }
  
  par(cex = 1, cex.main = 1.2, cex.lab = 1.1)
  for(i in 1:numtols){
      if(probs){
          mym <- lapply(model.probs, apply, 2, tapply, true, mean)
          barplot(t(mym[[paste("tol", tols[i], sep="")]]), ylab="Mean model probability", ...)
      }
      else{
          mym <- lapply(estim, table, true)
          barplot(mym[[paste("tol", tols[i], sep="")]], ylab="Frequency", ...)
      }
      title(caption, sub=paste("Tolerance rate = ", tols[i], sep=""))
  }
  
  if(!is.null(file)){
      dev.off()
  }
  else devAskNewPage(save.devAskNewPage)
  invisible()
  
}



