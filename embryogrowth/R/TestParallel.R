#' TestParallel estimates the likelihood of a set of parameters for nest incubation data with or without parallel computing option 
#' @title Estimate the likelihood of a set of parameters for nest incubation data with or without parallel computing option
#' @author Marc Girondot
#' @return The gain or loss of computing time using parallel version
#' @param result A object obtained after searchR or likelihoodR
#' @description Estimate the likelihood of a set of parameters for nest incubation data with or without parallel computing option.
#' It uses the user time from the print result of system.time() function.
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(resultNest_4p)
#' TestParallel(resultNest_4p)
#' }
#' @export


TestParallel <-
	function(result=stop("A ResultNest object must be provided")) {

print("Run with parallel computing")
x1 <- summary(system.time(likelihoodR(result, parallel=TRUE, hessian=FALSE)))
print(paste("End in", sprintf("%.3f", x1[[1]]),"seconds."))
print("Run without parallel computing")
x2 <- summary(system.time(likelihoodR(result, parallel=FALSE, hessian=FALSE)))
print(paste("End in", sprintf("%.3f", x2[[1]]),"seconds."))

if (x1[[1]]<x2[[1]]) {
	print("The option of parallel computing is interesting for you.")
	print(paste("You will gain", sprintf("%.3f",100*(x2[[1]]-x1[[1]])/x2[[1]]), "% of time to use it in searchR() and likelihoodR()."))
} else {
	print("The option of parallel computing is not interesting for you.")
	print(paste("You will lost", sprintf("%.3f", 100*(x1[[1]]-x2[[1]])/x2[[1]]), "% of time to use it in searchR() and likelihoodR()."))
}
return(invisible(100*(x1[[1]]-x2[[1]])/x2[[1]]))
}
