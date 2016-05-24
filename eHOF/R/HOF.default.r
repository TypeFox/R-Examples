HOF.default <- function(
		occ, 
		grad, 
		M = max(occ), 
		y.name, 
		family=binomial, 
		lim=100, 
		bootstrap=100, 
		test = c('AICc', 'BIC', 'AIC','Dev'), 
    modeltypes=eHOF.modelnames,
		...)  {
  if(any(c('data.frame', 'matrix','list') %in% class(occ))) stop('Performance data for HOF.default must be a vector.')
  x.name <- deparse(substitute(grad))
  if (missing(y.name)) y.name <- deparse(substitute(occ))
  if(any(is.na(occ))) stop('NA in occurrence vector is not allowed!')
  if(!is.numeric(grad)) print('Gradient must be a numeric vector')
  if(!is.null(bootstrap))
    if(bootstrap == 0) stop('If you do not want to bootstrap your data use "bootstrap=NULL not 0"!') else
      options(eHOF.bootselectmessage = FALSE)
  out <- HOF.model(occ, grad, M, y.name, x.name, family=family, lim=lim,...)
  
  IC.weights <- function(x, test = 'AICc') {
    p <- sapply(x$models, function(x) length(x$par))
    k <- if(test == 'BIC') log(x$nobs) else 2
    if(test == 'AICc') ic <- -2 * logLik(x) + k * p + (2*k*(k + 1))/(x$nobs - k - 1)
    if(test %in% c('AIC', 'BIC'))   ic <- -2 * logLik(x) + k * p
    if (test == "Dev")   ic <- deviance(x)
    ic.W <- round(exp(-0.5 * ic)/ sum(exp(-0.5 * ic), na.rm=TRUE), 4)
	  return(ic.W)
  }
  
  if(!is.null(bootstrap)) {
    test <- match.arg(test)
    rejectedmodels <- sapply(out$models, function(x) is.na(x$deviance))
    modeltypes <- modeltypes[!modeltypes %in% eHOF.modelnames[rejectedmodels]]
    bootmodels <- character(length=bootstrap)
    mods <- vector('list', length=bootstrap)
	  weights <-  matrix(nrow=bootstrap, ncol=7); colnames(weights) <- eHOF.modelnames
    pb <- txtProgressBar (min = 0, max = bootstrap, char = '.',  width = 45, style = 3)
    for(i in 1:bootstrap) {
      take <- sample(length(grad), replace=TRUE)
      mods[[i]] <- HOF.model(occ[take], grad[take], M=M, y.name, x.name, bootstrap=NULL, family=family, lim=lim, ...)
      bootmodels[i] <- pick.model(mods[[i]], test=test, selectMethod = 'pick.model', modeltypes=modeltypes, ...)
	    weights[i,] <- IC.weights(mods[[i]], test=test)
      setTxtProgressBar(pb, bootstrap - (bootstrap - i))
#      for(m in 1:7) mods[[i]]$models[[m]]['fitted'] <- NULL
    }
    close (pb) ## Close progress bar
    out$call <- match.call()
#    out$bootmods <- mods
    out$bootstraptest <- test
    out$bootstrapmodels <- bootmodels
	  out$ICweights <- weights
  } else options(bootselectmessage = TRUE)
  out
}

