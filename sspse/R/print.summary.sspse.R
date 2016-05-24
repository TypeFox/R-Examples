#' Summarizing Population Size Estimation Model Fits
#' 
#' This is the \code{print} method for the \code{summary} class method for class \code{"sspse"} objects.
#' These objects encapsulate an estimate of the posterior distribution of
#' the population size based on data collected by Respondent Driven Sampling.
#' The approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008).  As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' 
#' \code{print.summary.sspse} tries to be smart about formatting the
#' coefficients, standard errors, etc. and additionally gives
#' \sQuote{significance stars} if \code{signif.stars} is \code{TRUE}.
#' 
#' Aliased coefficients are omitted in the returned object but restored by the
#' \code{print} method.
#' 
#' Correlations are printed to two decimal places (or symbolically): to see the
#' actual correlations print \code{summary(object)$correlation} directly.
#' 
#' @method print summary.sspse
#' @export
#' @aliases print.summary.sspse
#' @param x an object of class \code{"summary.sspse"}, usually, a result of a
#' call to \code{summary.sspse}.
#' @param digits the number of significant digits to use when printing.
#' @param correlation logical; if \code{TRUE}, the correlation matrix of the
#' estimated parameters is returned and printed.
#' @param covariance logical; if \code{TRUE}, the covariance matrix of the
#' estimated parameters is returned and printed.
#' @param signif.stars logical. If \code{TRUE}, \sQuote{significance stars} are
#' printed for each coefficient.
#' @param eps.Pvalue number; indicates the smallest p-value. 
#' \code{\link[stats]{printCoefmat}}.
#' @param \dots further arguments passed to or from other methods.
#' @return The function \code{summary.sspse} computes and returns a two row matrix of
#' summary statistics of the prior and estimated posterior distributions. 
#' The rows correspond to the \code{Prior} and the \code{Posterior}, respectively. 
#' The rows names are \code{Mean}, \code{Median}, \code{Mode}, 
#' \code{25\%}, \code{75\%}, and \code{90\%}.
#' These correspond to the distributional mean, median, mode, lower quartile,
#' upper quartile and 90\% quantile, respectively. 

#' @seealso The model fitting function \code{\link{posteriorsize}},
#' \code{\link{summary}}.
#' 
#' Function \code{\link[stats]{coef}} will extract the matrix of coefficients with
#' standard errors, t-statistics and p-values.
#' @references
#'
#' Gile, Krista J. (2008) \emph{Inference from Partially-Observed Network
#' Data}, Ph.D. Thesis, Department of Statistics, University of Washington.
#' 
#' Gile, Krista J. and Handcock, Mark S. (2010) \emph{Respondent-Driven
#' Sampling: An Assessment of Current Methodology}, Sociological Methodology
#' 40, 285-327.
#' 
#' Gile, Krista J. and Handcock, Mark S. (2014) \pkg{sspse}: Estimating Hidden 
#' Population Size using Respondent Driven Sampling Data
#' R package, Los Angeles, CA.  Version 0.5, \url{http://hpmrg.org}.
#' 
#' Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
#' \url{http://statnetproject.org}.
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2014)
#' \emph{Estimating Hidden Population Size using Respondent-Driven Sampling
#' Data}, Electronic Journal of Statistics, 8, 1, 1491-1521
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2015)
#' \emph{Estimating the Size of Populations at High Risk for HIV using Respondent-Driven 
#' Sampling Data}, Biometrics.
#' @keywords models
#' @examples
#' 
#' \dontrun{
#' N0 <- 200
#' n <- 100
#' K <- 10
#' 
#' # Create probabilities for a Waring distribution 
#' # with scaling parameter 3 and mean 5, but truncated at K=10.
#' probs <- c(0.33333333,0.19047619,0.11904762,0.07936508,0.05555556,
#'            0.04040404,0.03030303,0.02331002,0.01831502,0.01465201)
#' probs <- probs / sum(probs)
#' 
#' # Look at the degree distribution for the prior
#' # Plot these if you want
#' # plot(x=1:K,y=probs,type="l")
#' # points(x=1:K,y=probs)
#' #
#' # Create a sample
#' #
#' set.seed(1)
#' pop<-sample(1:K, size=N0, replace = TRUE, prob = probs)
#' s<-sample(pop, size=n, replace = FALSE, prob = pop)
#'  
#' out <- posteriorsize(s=s,interval=10)
#' plot(out, HPD.level=0.9,data=pop[s])
#' summary(out, HPD.level=0.9)
#' # Let's look at some MCMC diagnostics
#' plot(out, HPD.level=0.9,mcmc=TRUE)
#' }
#' 
print.summary.sspse <- function (x,
              digits = max(3, getOption("digits") - 3),
              correlation=FALSE, covariance=FALSE,
              signif.stars= getOption("show.signif.stars"),
              eps.Pvalue=0.0001, ...)
{
  cat("\n==========================\n")
  cat("Summary of model fit\n")
  cat("==========================\n\n")
  
  cat("Formula:   ")
  print(x$formula)
  cat("\n")
  
  cat ("Newton-Raphson iterations: ", x$iterations, "\n")
  if(x$pseudolikelihood){
    if (x$independence) {
      cat ("\nMaximum Likelihood Results:\n")
    } else {
      cat ("\nMaximum Pseudolikelihood Results:\n")
    }
  }else{
    cat ("MCMC sample of size", x$samplesize, "\n")
    cat ("\nMonte Carlo MLE Results:\n")
  }
    
  if(!is.null(x$randomeffects)){ 
     if(!is.matrix(x$randomeffects)){
       cat ("\n Activity random effects:\n  Variances:\n")
       print(x$randomeffects)
     }else{
      cat ("\nSender and Receiver random effects:\n  Covariances:\n")
      print(x$randomeffects)
      corr <- x$randomeffects[1,2]/sqrt(x$randomeffects[1,1]*x$randomeffects[2,2])
      corr <- max(min(1,corr),-1)
      cat (paste("\n  Correlation between sender and receiver:  ",
          round(corr,5)),"\n\n")
     }
  }

  printCoefmat(x$coefs, digits=digits, signif.stars=signif.stars,
               P.values=TRUE, has.Pvalue=TRUE, na.print="NA",
               eps.Pvalue=eps.Pvalue, ...)
  
  if(!is.null(x$message)){ 
     cat(x$message)
  }

  cat("\n")
  cat(x$devtable)

  cat(paste("AIC:", format(x$aic, digits = 5), "  ", 
            "BIC:", format(x$bic, digits = 5), "\n", sep=" "))
  

  if(any(x$drop)){
    cat("\n Warning:\n")
    for(i in names(x$coefs[x$offset,1])){
     cat(paste("  The term",i,
     "is degenerate and has an infinite coefficient estimate.\n",
      sep=" "))
    }
  }

  if(any(x$offset&!x$drop)){
    cat("\n Warning:\n")
    for(i in names(x$coefs[x$offset,1])){
    cat(paste("  The term",i,
     "has been offset and was not estimated from the data.\n",
      sep=" "))
    }
  }

  if(!is.null(x$degeneracy.value) && !is.na(x$degeneracy.value)){
   if(is.infinite(x$degeneracy.value)){
    cat("\n Warning: The diagnostics indicate that the model is very unstable.\n   They suggest that the model is near degenerate,\n   and that the numerical summaries are suspect.\n")
   }else{
    if(x$degeneracy.value > 1){
      cat("The instability of the model is: ",
        format(x$degeneracy.value, digits=2),"\n")
      cat("Instabilities greater than 1 suggest the model is near degenerate.\n")
    }
   }
  }

  if (covariance == TRUE)
    {
      cat("Asymptotic covariance matrix:\n")
      print(x$asycov)
    }
  
  if (correlation == TRUE)
    {
      cat("\nAsymptotic correlation matrix:\n")
      asycor <- x$asycov / crossprod(x$asyse)
      dimnames(asycor) <- dimnames(x$asycov)
      print(asycor)
    }
  
  invisible(x)
}
