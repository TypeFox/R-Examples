#' Print method for objects of class \code{sscsample}
#' 
#' This function provides a print summary method for the output of
#' \code{sscsample}. The \code{sscsample} produces a large number of samples
#' from a fixed population using either simple random, stratified, or cluster
#' sampling. This function provides the means of each sample plus the number of
#' observations from each ethnicity stratum in the sample.
#' 
#' 
#' @param x an object of class \code{sscsamp} produced by \code{sscsample}
#' @param \dots any other arguments that are to be passed to \code{cat}
#' @author James Curran
#' @seealso \code{\link{sscsample}}
#' @export
print.sscsamp = function(x, ...){
  cat("Sample   Mean    Stratum 1  Stratum 2  Stratum 3\n", ...)
  cat("------  -------  ---------  ---------  ---------\n", ...)

  n.samples = length(x$means)
  fmt = '%6d  %7.4f  %9d  %9d  %9d\n'
  
  for (r in 1:n.samples) {
    s = sprintf(fmt, r, round(x$means[r], 4), x$s.strata[r,1],
                x$s.strata[r, 2], x$s.strata[r, 3])
    cat(s, ...)
  }
}
