#'@title Bootstrapping of matrixpls function
#'
#'@description
#'\code{matrixpls.boot} is a convenience method that implements bootstrapping of \code{matrixpls}
#'with \code{\link[boot]{boot}} method of the \code{boot} package.
#'
#'@section Analyzing the boostrap results:
#'
#'The output can be analyzed with any of the functions provided by the \code{\link[boot]{boot}} package.
#'For example, the \code{\link[boot]{boot.ci}} function can be used for calculating confidence 
#'intervals and the \code{\link[boot]{empinf}} function can be used to calculate influence values that may
#'be useful for identifying outliers.
#'
#'The class \code{matrixplsboot} provides only a \code{summary} function, which calculates a set
#'of statistics that are commonly of interest after boostrapping. This includes standard errors,
#'t statistics (estimate / SE) and p-values based on Student's t distribution and standard normal
#'distribution. Because the sampling distribution of the parameter estimates calculated by 
#'\code{matrixpls} are not always known, the p-values cannot be expected to be unbiased.
#'
#'This concern
#'applies particularly when using PLS weights. Because the PLS literature provides conflicting
#'advice on which probability distribution to use as a reference, the \code{summary} method of 
#'\code{matrixplsboot} produces two-tailed p-values based on four different probalility
#'distributions. The \emph{regression} p values are based on comparing the t statistic against
#'the references distribution used in regression analysis, namely Student's t distribution with 
#'\code{n - k - 1} degrees of freedom. The \emph{Hair} p values are based on 
#'Hair et al's (2014, p. 134) recommendation to ignore the number of independent variables 
#'\code{k} and set the degrees of freedom to \code{n - 1}. The \emph{Henseler} p values are based
#'on the recommendation by
#'Henseler et al (2009, p. 305) that the degrees of freedom should be set as \code{n + m - 2},
#'where \code{m}is always 1 and \code{n} is the number of bootstrap samples. The \code{z} p values
#'are based on comparing the t statistic against the standard normal distibution. This choice
#'can be motivated by asymptotic normality of the PLS estimates in certain conditions.
#'
#'
#'@param data Matrix or data frame containing the raw data.
#'
#'@param R Number of bootstrap samples to draw.
#'
#'@param ... All other arguments are passed through to \code{\link{matrixpls}}.
#'
#'@param signChangeCorrection Sign change correction function.
#'
#'@param extraFun A function that takes a \code{matrixpls} object and returns a numeric vector. The
#'vector is appended to bootstrap replication. Can be used for boostrapping additional
#'statistics calculated based on the estimation results.
#'
#'@param stopOnError A logical indicating whether boostrapping should be continued when error occurs
#' in a replication.
#'
#'@inheritParams boot::boot
#'
#'
#'@return An object of class \code{matrixplsboot} and \code{\link[boot]{boot}}.
#'
#'@seealso
#'\code{\link[boot]{boot}}
#'
#'Sign change corrections: \code{\link{signChange.individual}}; \code{\link{signChange.construct}}
#'
#'@export
#'
#'@example example/matrixpls.boot-example.R 
#'
#'@reference
#'
#'Hair, J. F., Hult, G. T. M., Ringle, C. M., & Sarstedt, M. (2014). \emph{A primer on partial least squares structural equations modeling (PLS-SEM)}. Los Angeles: SAGE.
#'
#'Henseler, J., Ringle, C. M., & Sinkovics, R. R. (2009). The use of partial least squares path modeling in international marketing. \emph{Advances in International Marketing}, 20, 277–319.
#'
#'Rönkkö, M., & Evermann, J. (2013). A critical examination of common beliefs about partial least squares path modeling. \emph{Organizational Research Methods}, 16(3), 425–448. https://doi.org/10.1177/1094428112474693
#'
#'Rönkkö, M., McIntosh, C. N., & Antonakis, J. (2015). On the adoption of partial least squares in psychological research: Caveat emptor. \emph{Personality and Individual Differences}, (87), 76–84. https://doi.org/10.1016/j.paid.2015.07.019
#'
#'
matrixpls.boot <- function(data, ..., R = 500,
                           signChangeCorrection = NULL,
                           parallel = c("no", "multicore", "snow"),
                           ncpus = getOption("boot.ncpus", 1L),
                           stopOnError = FALSE,
                           extraFun = NULL){
  
  if(! requireNamespace("boot")) stop("matrixpls.boot requires the boot package")
  
  if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
  
  data <- as.matrix(data)

  arguments <- list(...)
  
  # Prepare sign change corrections
  
  if(! is.null(signChangeCorrection)){
    Worig <- attr(matrixpls(stats::cov(data),...), "W")
    
    # Get the original weight function
    weightFunction <- arguments[["weightFunction"]]
    if(is.null(weightFunction)) weightFunction <- weight.pls
    
    # Wrap inside sign change correction
    arguments[["weightFunction"]] <- function(S, ...){
      Wrep <- weightFunction(S, ...)
      W <- signChangeCorrection(Worig, Wrep)
      W <- scaleWeights(S, W)
      attributes(W) <- attributes(Wrep)
      W
    }
    
  }
  
  # Bootstrap
  
  boot.out <- boot::boot(data,
                   function(data, indices){
                     
                     S <- stats::cov(data[indices,])
                     arguments <- c(list(S),arguments)
                     
                     if(stopOnError){
                       boot.rep <- do.call(matrixpls, arguments)
                     }
                     else{
                       tryCatch(
                         boot.rep <- do.call(matrixpls, arguments)
                       )
                     }
                     
                     # Add additional statistics
                     
                     if(!is.null(extraFun)){
                       a <- attributes(boot.rep)
                       boot.rep <-c(boot.rep,extraFun(boot.rep))
                       a$names <- names(boot.rep)
                       attributes(boot.rep) <- a
                     }
                     
                     # If the indices are not sorted, then this is not the original sample
                     # and we can safely omit all attributes to save memory
                     
                     if(is.unsorted(indices)) attributes(boot.rep) <- NULL
                     
                     boot.rep
                   },
                   R, parallel = parallel, ncpus = ncpus)
  
  class(boot.out) <- c("matrixplsboot", class(boot.out))
  boot.out
}

# These are not in use 

#'@S3method print matrixplsboot

print.matrixplsboot <- function(x, ...){
  matrixpls.out <- x$t0
  attr(matrixpls.out,"boot.out") <- x
  print(matrixpls.out, ...)
}

#'@S3method summary matrixplsboot

summary.matrixplsboot <- function(object, ...){
	matrixpls.res <- object$t0
	out <- summary(matrixpls.res)
	attr(out,"boot.out") <- object
	
	cat("\nCalculating confidence intervals.\n")
	
	# Omit CIs for the weights
	parameterIndices <- 1:(length(matrixpls.res) -
	          sum(attr(matrixpls.res,"W")!=0))
	
	cis <- lapply(parameterIndices, function(i){
	  ci <- boot::boot.ci(object,...,index = i, type = c("norm","basic", "perc"))
	  cis <- c(matrixpls.res[i],ci$normal[2:3],ci$basic[4:5], ci$percent[,4:5],NA,NA)
	  
	  # BCa intervals cannot be always calculated
	  
	  tryCatch(
	    cis[8:9] <- boot::boot.ci(object,...,index = i, type = "bca")$bca[,4:5],
	    error = function(e){
	      warning(e)
	    })
	  
	  cis
	})
	
	cis <- do.call(rbind,cis)
	rownames(cis) <- names(matrixpls.res)[parameterIndices]
	colnames(cis) <- c("Estimate","Norm low", "Norm up","Basic low", "Basic up",
	                   "Perc low", "Perc up", "BCa low", "BCa up")
	
	out$ci <- cis
	
	# P-values
	
	# There is some disagreement on the degrees of freedom of the t-statistic in the PLS literature
	
	# Hair et al. (2014, p. 134) states that "the test statistic follows a t distribution with degrees of
	# freedom […] equal to the number of observations minus 1" 
	
	dfHa <- nrow(object$data) - 1
	  
	# Henseler et al. (2009, p. 305) "the degrees of freedom for the test is  n + m – 2, where m is
	# always 1 and n is the number of bootstrap samples. 
	
	dfHe <- object$R + 1 - 2
	
	# We need to get the number of IVs in the regression to calculate the regression p value
	
	IVs <- unlist(lapply(attr(matrixpls.res, "model"),function(modelmat){
	  ivs <- rowSums(modelmat)
	  mm <- sweep(modelmat,1,ivs,"*")
	  mm[mm!=0]
	}))
	
	ps <- lapply(parameterIndices, function(i){
	  est <- object$t0[i]
	  se <- stats::sd(object$t[,i])
	  t <- est/se
	  
	 c(est,se,t,
	    (1-stats::pt(t,dfHa-IVs[i]))*2, # Regression
	    (1-stats::pt(t,dfHa))*2, # Hair
	    (1-stats::pt(t,dfHe))*2, # Henseler
	    (1-stats::pnorm(t))*2) # Standard normal
	 
	})
	
	ps <- do.call(rbind,ps)
	rownames(ps) <- names(matrixpls.res)[parameterIndices]
	
	colnames(ps) <- c("Estimate","SE","t","p (regression)","p (Hair)", "p (Henseler)","p (z)")
	out$p <- ps
	
	class(out) <- "matrixplsbootsummary"
	out
}

#'@S3method print matrixplsbootsummary

print.matrixplsbootsummary <- function(x, ...){
  class(x) <- "matrixplssummary"
  print(x)
  
  cat("\n Bootstrap SEs and significance tests\n")
  print(x$p, ...)

  cat("\n Bootstrap confidence intervals\n")
  ci <- x$ci
  ci <- data.frame(Estimate = ci[,1],
                " (",
                ci[,2],
                Normal=ci[,3],
                ") (",
                ci[,4],
                Basic=ci[,5],
                ") (",
                ci[,6],
                Percentile=ci[,7],
                ") (",
                ci[,8],
                BCa=ci[,9],
                ")")
  colnames(ci)[c(2,3,5,6,8,9,11,12,14)] <- " "
  print(ci, ...)
  
}

