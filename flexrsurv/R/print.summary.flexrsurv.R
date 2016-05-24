
print.summary.flexrsurv <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    symbolic.cor = x$symbolic.cor,
                                    signif.stars = getOption("show.signif.stars"),
                                    ...)
{

  model <- x$call[["model"]]
  if (is.null(x$call[["model"]])) {
    model <- eval(formals(flexrsurv)$model)[1]
  }

  
  cat("\nCall: \n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  
    if(length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
      cat("\nCoefficients:\n")
      coefs <- x$coefficients
      aliased <- x$aliased
      if(!is.null(aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4L,
                        dimnames=list(cn, colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
      }
      printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
    }

  if("flexrsurv.mle" %in% attr(x, "fitclass")){
    if(is.null(dim(x$cov))){
        cat("\nUnable to compute the estimated covariance matrix of the estimated coefficients.\n")
      }
    else if ( attr(x$cov, "type") == "none"){
        cat("\nThe estimated covariance matrix has not been computed; use vartype = \"oim\" or vartype = \"opg\".\n")
      } 
  } else if ("flexrsurv.glmiterative" %in% attr(x, "fitclass") ){
      cat('\nThe estimated covariance matrix of the estimated coefficients is not available with method = "glm" and multiplicative NPHNLL() effect.\n')
      cat("\nYou can obtain the estimated covariance matrix by recalling flersurv() with the same options except for method = \"MLE\", int_method = \"GLM\".\n") 
    }

#  else {
    # code from #  File src/library/stats/R/glm.R 
    correl <- x$correlation
    if(!is.null(correl)) {
	p <- NCOL(correl)
	if(p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {
		print(symnum(correl, abbr.colnames = NULL))
	    } else {
		correl <- format(round(correl, 2L), nsmall = 2L,
                                 digits = digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote = FALSE)
	    }
	}
    }
    cat("\n")
#  }

  if(!is.null(x$optim.control$reltol)){
    digitsll <- round(abs(log10(x$optim.control$reltol)))
    if ( digitsll < 0) {
      digitsll <- digits
    }
  } else {
      digitsll <- max(digits, 8)
  }
  cat("\nLog-likelihood:",  format(x$loglik, digits = digitsll), "\n")
 
  cat("\n")
  invisible(x)
  
}

