#'@method print prof
#'@export

print.prof <- function(x, ...){
  if (!inherits(x, "prof")) stop("Use only with 'prof' objects.\n")
	cat("Subscore Reliability Estimates:\n")
	cat("\n")
	print(x$reliability)
  invisible(x)
}
