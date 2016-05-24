#' Print cds Object
#' 
#' Print method for \code{cds} objects.
#' 
#' @param x A \code{cds} object.
#' @param \dots Unimplemented.
#' @keywords utility
#' @method print cds
#' @export
print.cds <- function(x, ...) {
  cat(paste("\t Constrained Dual Scaling with", x$K, "groups\n\n"))
  cat(paste("Minimum loss: ", format(x$minloss), ";\t Nr of iterations: ", 
            x$iter,"\n\n", sep = ""))
  cat(paste("\tDistribution of ", length(x$loss.G), " loss values (", length(unique(x$loss.G))," unique):\n", sep = ""))
  print(quantile(x$loss.G))
  cat("\n\tEstimated alpha's: \n\n")
  ratios <- with(x, cbind(alphamat[,3]/alphamat[,2], alphamat[,3]/alphamat[,4]))
  rstype <- cut(atan2(y = ratios[, 2] - 2, x = ratios[, 1] - 2), c(-pi, -pi/2, 0, pi/2, pi), 
                labels = c("Extreme","Acquiescence", "Midpoint", "Disacquiescence"))
  tmp <- data.frame(round(x$alphamat, 3), Style = rstype)
  print(tmp)
  if (!is.null(x$confusion)) {
    cat("\n\tConfusion matrix:\n(actual/estimated)\n")
    print(x$confusion)
    cat(paste("\n\tHitrate:", round(100*x$hitrate, 3), "%\n"))
  }
  cat("\n\tGroup sizes:\n")
  print(table(x$grp))
  cat("\nCalculation time:\n")
  cat(round(x$time.total, 1), "seconds\n")
}
