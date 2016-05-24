"print.HOF" <- function (
		x,
		...,
		test = 'AICc', 
		penal = 'df', 
		selectMethod = c('bootselect.lower', 'IC.weight', 'pick.model'), 
		gam = FALSE, 
		k=4) {
	  selectMethod <- match.arg(selectMethod)
    if(length(penal)==1) {
      if(penal == 'df') penal <- sapply(x$models, function(x) length(x$par))
# if(penal == 'new') penal <- c('I'=1,'II'=2,'III'=2,'IIIb'=3,'IV'=3,'IVb'=3,'V'=4,'Vb'=4)    
# penal <- c('I'=1,'II'=2,'III'=2,'IIIb'=3,'IV'=3,'IVb'=3,'V'=4,'Vb'=4) else
# possible manipulation of model penalization
    }
    cat('Response of: ', x$y.name, "\n")
    cat('Deviances and information criteria:\n')
    if(gam) {
      gamfun <- function(x, bs = 'cr', ...) gam(x$y ~ s(x$x, bs=bs, k=k, ...),family=get(x$family), scale=0, ...)
      pg <- gamfun(x, ...)
      x$models$GAM <- pg
      penal <- c(penal, GAM = sum(pg$edf) )
    }
    dev <- deviance(x)
    ll <- logLik(x)
    AICc <- -2 * ll + 2 * penal + 2 * penal *(penal + 1)/(x$nobs - penal - 1) 
    d.AICc <- AICc - min(AICc, na.rm=TRUE)
    # nAICc <- AICc/sum(AICc, na.rm=TRUE)
    AICc.W <- round(exp(-0.5*AICc)/ sum(exp(-0.5*AICc), na.rm=TRUE),4)
	  BIC  <- -2 * ll + log(x$nobs) * penal
    d.BIC <- BIC - min(BIC, na.rm=TRUE)

    out <- cbind(Deviance = dev, logLik=ll, AICc=AICc, AICc.Diff = d.AICc, AICc.W = AICc.W, BIC.Diff = d.BIC)
    printCoefmat(out, na.print="")
    if(!is.null(x$bootstrapmodels)) {
      cat('Percentage of model types after bootstrapping:')
      print(round(table(x$bootstrapmodels)/length(x$bootstrapmodels)*100))
  	  cat('Sum of bootstrapped model weights:\n')
  	  print(round(colSums(x$ICweights, na.rm=TRUE),2))

	  if(selectMethod == 'bootselect') {
		  best <- pick.model(x, k=k, selectMethod = selectMethod, silent = TRUE, ...)
		  cat("\nTest used during bootstrapping: ", x$bootstraptest, sep='')
	  }
	  if(selectMethod == 'IC.weight') best <- names(which.max(colSums(x$ICweights))) else
                            best <- pick.model(x, test=test, k=k, silent = TRUE, ...)
    
     cat("\nSuggested best model (", test, ", " , selectMethod, "): ", best, '\n', sep='')
    } else {
  	  best <- pick.model(x, test=test, k=k, silent = TRUE, ...)
  	  cat("\nSuggested best model (",test, "information criterion ): ", best, '\n\n')
  }
	    invisible(out)
}
