summary.randomLCA <-
function(object,...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
    out <- list()
	out$probit <- object$probit
	out$logLik <- logLik(object)
	out$AIC <- AIC(object)
	out$BIC <- BIC(object)
	out$AIC3 <- AIC3(object)
	out$penlogLik <- object$penlogLik
	out$nclass <- object$nclass
	out$classp <- object$classp
	names(out$classp) <- paste("Class ",1:object$nclass)
	out$outcomep <- as.data.frame(object$outcomep)
	row.names(out$outcomep) <- paste("Class ",1:object$nclass)
	names(out$outcomep) <- names(object$patterns)
	out$random <- object$random
	out$level2 <- object$level2
	  if (object$constload) blocksize <-  1
	  else {
		if (object$level2) blocksize <- object$level2size
		else blocksize <- object$blocksize
	  }  
	if (object$random) {
		out$lambdacoef <- object$lambdacoef
		if (!object$byclass) out$lambdacoef <- t(out$lambdacoef)
		if (object$byclass) names1 <- paste("Class ",1:object$nclass)
		else names1 <- ''
    #browser()
		names2 <- names(object$patterns)[1:blocksize]
		names2 <- strsplit(names2,"\\.")
		x <- NULL
    #browser()
		for (i in 1:blocksize) {
			x <- c(x,names2[[i]][1])
		}
		if (blocksize==1) names2 <- ""
		else names2 <- x
    #browser()
		dimnames(out$lambdacoef) <- list(names1,names2)
		if (object$level2) {
			out$taucoef <- object$taucoef
			if (object$byclass) {
				names(out$taucoef) <- paste("Class ",1:object$nclass)
			} else {
				names(out$taucoef) <-''
			}
		}
		margp <- calcMargProb(object)
		out$margoutcomep <- data.frame(t(matrix(margp$outcomep,ncol=object$nclass)))
		row.names(out$margoutcomep) <- paste("Class ",1:object$nclass)
		names(out$margoutcomep) <- names(object$patterns)
	}
	class(out) <- "summary.randomLCA"
	out
}


print.summary.randomLCA <- function(x, ...)
{
  if (!inherits(x, "summary.randomLCA"))
    stop("Use only with 'randomLCA summary' objects.\n")
  if (x$probit) link <- "Probit"
	else link <- "Logit"
	if (x$random) print(data.frame(Classes = x$nclass, AIC = x$AIC, BIC = x$BIC, AIC3 = x$AIC3,
		logLik = c(x$logLik),penlogLik = c(x$penlogLik),Link=link,row.names = " ") )
	else print(data.frame(Classes = x$nclass, AIC = x$AIC, BIC = x$BIC, AIC3 = x$AIC3,
		logLik = c(x$logLik),penlogLik = c(x$penlogLik),row.names = " ") )
    cat("Class probabilities","\n")
	print(round(x$classp,4))
	if (x$random)  cat("Conditional outcome probabilities","\n")
    else cat("Outcome probabilities","\n")
	print(round(x$outcomep,4))
	if (x$random) {
		cat("Marginal Outcome Probabilities","\n")
		print(round(x$margoutcomep,4))
		cat("Loadings","\n")
		if (length(x$lambdacoef)==1) cat(sprintf("%g\n",x$lambdacoef))
		else print(round(x$lambdacoef,4))
		if (x$level2) {
			cat("Tau","\n")
			if (length(x$taucoef)==1) cat(sprintf("%g\n",round(x$taucoef,4)))
			else print(round(x$taucoef,4))
		}		
	}
	invisible(x)
}

print.randomLCA <- function(x, ...)
{   
  if (!inherits(x, "randomLCA"))
    stop("Use only with 'randomLCA' objects.\n")
  print(summary(x),...)
  invisible()
}
  