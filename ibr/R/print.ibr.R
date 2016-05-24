print.ibr <- function(x, digits = max(2, getOption("digits") - 4), ...) {
  cat("\nInitial df:", format(round(x$initialdf,digits)), "; Final df:", format(round(x$finaldf,digits)), "\n")
 if (x$parcall$critmethod=="aggregation") crit4iter <- paste("aggregation of",x$parcall$criterion) else crit4iter <- x$parcall$criterion[1]
  cat("Number of iterations:", x$iter, "chosen by", crit4iter, "\n")
}
