#' Sensitivity indices based on simple and standard regression analyses
#' @param X a design dataframe
#' @param Y a dataframe of responses
#' @param model the regression model
#' @value a list of regression results
#' @note This function is essentially a (possibly repeated) call to the "src" function of the R sensitivity package

regressionSI <- function(X, Y, rank=FALSE, nboot=100, conf=0.95, ...){

  ## PRELIMINARIES
  nf <- ncol(Y)
  Ynames <- names(Y)
  results <- vector(mode="list", length=nf)

  ## MAIN CALCULATIONS
  for(i in seq(nf)){
    resultsA <- src(X, Y[,i], rank=rank, nboot=nboot, conf=conf)
    if(!rank)
      results[[i]] <- resultsA[c("SRC")]
    else
      results[[i]] <- resultsA[c("SRRC")]
  }
  
  ## Output
  output <- list(main=results, information=list(Ynames=Ynames,rank=rank,nboot=nboot,conf=conf))
  return(output)
}

print.regressionSI <- function(x, ...){
  ## PRELIMINARIES
  results <- x$main
  Ynames <- x$information$Ynames
  rank <- x$information$rank
  nf <- length(results)

  ## MAIN CALCULATIONS
  cat("\nStandardized Regression Coefficients (SRC)\n")
  cat(" calculated by the src function of the sensitivity library\n")
  cat(" with parameters rank = ", x$information$rank,
      ", nboot = ", x$information$nboot,
      ", conf = ", x$information$conf,
      "\n\n")
  for(i in seq(nf)){
    cat("Response variable: ", Ynames[i], "\n") 
    print(results[[i]])
  }
  invisible()
}

plot.regressionSI <- function(x,y, ...){
  ## PRELIMINARIES
  results <- x$main
  Ynames <- x$information$Ynames
  rank <- x$information$rank
  nf <- length(results)
 
  ## MAIN CALCULATIONS
  nb.col <- ceiling(sqrt(nf))
  nb.row <- ceiling(nf/nb.col)
  keep.mfrow <- par(mfrow=c(nb.row,nb.col))

  for(i in seq(nf)){
    src.obj.i <- results[[i]]
    class(src.obj.i) <- "src"
    plot(src.obj.i, ...)
    title(xlab=Ynames[i])
    abline(h=0)
  }
  par(keep.mfrow)
  
  invisible()
}

