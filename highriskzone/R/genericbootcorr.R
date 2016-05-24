#' Print Brief Details of a bootstrap correction for a high-risk zone
#'
#' Prints a very brief description of the bootstrap correction for a high-risk zone.
#'
#' A very brief description of the bootstrap correction x for a high-risk zone is printed.
#' This is a method for the generic function \code{\link[base]{print}}.
#'
#' @param x bootstrap correction for of a high-risk zone (object of class "\code{bootcorr}")
#' @param ... ignored
#' @method print bootcorr
#' @export print bootcorr
#' @seealso \code{\link[base]{print}}, \code{\link{summary.bootcorr}}

print.bootcorr <- function(x, ...){
  cat("resulting value for alpha (cutoff):", x$alphastar, " \n")
}





#' Summary of a the bootstrap correction for a high-risk zone
#'
#' Prints a useful summary of the bootstrap correction for a high-risk zone.
#'
#' A useful summary of the bootstrap correction x for a high-risk zone is printed.
#' This is a method for the generic function \code{\link[base]{summary}}.
#'
#' @param object bootstrap correction for a high-risk zone (object of class "\code{bootcorr}")
#' @param ... ignored
#' @method summary bootcorr
#' @export summary bootcorr
#' @seealso \code{\link[base]{summary}}, \code{\link{print.bootcorr}}, \code{\link{plot.bootcorr}}

summary.bootcorr <- function(object, ...){
  cat("resulting value for alpha (cutoff):", object$alphastar, " \n \n")
  
  cat("values for alpha (cutoff) which were tested: \n \n")
  
  numtest <- max(object$course$k)
  
  for(j in 1:numtest){
    dataj <- subset(object$course, object$course$k==j)
    maxi <- max(dataj$i)
    alphastarj <- max(dataj$alphastar)
    maxnumoutj <- max(dataj$numout)
    cat(alphastarj, " (", maxi, " iterations, final numout: ", maxnumoutj, ") \n", sep="")
  }
  
}




#' Visualize the bootstrap correction for a high-risk zone.
#'
#' Plot a visualization of the bootstrap correction for a high-risk zone.
#' The different values tested for alpha are plotted.
#'
#' This is the plot method for the class \code{bootcorr}.
#'
#' @param x bootstrap correction for a high-risk zone (object of class "\code{bootcorr}")
#' @param ... extra arguments passed to the generic \code{\link[graphics]{plot}} function.
#' @method plot bootcorr
#' @importFrom graphics points
#' @export plot bootcorr
#' @seealso \code{\link[graphics]{plot}}, \code{\link{print.bootcorr}}, \code{\link{summary.bootcorr}}


plot.bootcorr <- function(x, ...){

  resultdf <- x$course
  resultdf$iter <- 1
  resultdf$isum <- cumsum(resultdf$iter)

  lastid <- length(resultdf$isum)


  plot(resultdf$isum, resultdf$alphastar, type="s", xlab="iteration", ylab="cutoff", ...)
  points(resultdf$isum[lastid], resultdf$alphastar[lastid], pch=16, ...)
  
}
