
#' Summarize Permutation Tests for fsQCA Data
#'
#' Displays observed values, confidence intervals, and raw and adjusted p-scores for both consistency and counterexamples following permutation test of fsQCA data.
#' @param object Object returned by \code{\link{fsQCApermTest}}.
#' @param ... Additional parameters to pass on.
#' @return Two matrices of values for counterexamples and consistency.
#' @keywords permutation test fsQCA
#' @export
#' @examples
#' data(social.revolutions)
#' attach(social.revolutions)
#'
#' intersect <- pmin(breakdown, pop.ins)
#' intersect2 <- pmin(breakdown, (1-pop.ins))
#' intersect3 <- pmin((1-breakdown), pop.ins)
#' intersect4 <- pmin((1-breakdown), (1-pop.ins))
#'
#' test <- fsQCApermTest(y=soc.rev, configs=list(BI=intersect, Bi=intersect2, 
#'    bI=intersect3, bi=intersect4), total.configs=4)
#' summary(test)




summary.fsQCApt <- function(object, ...){
	cat("Call:\n")
    print(object$call)
    cat("\n")
    cat("Counterexamples\n")
	printCoefmat(object$result.cex[,c(1,3,2,5,7)])
    cat("\n")
    cat("Consistency\n")
	printCoefmat(object$result.con[,c(1,2,3,5,7)])
	cat("\n")
    cat(paste("Total number of configurations:", object$total.configurations), "\n")
    cat(paste("Number of permutations:", length(object$permutations.cex[[1]])), "\n")
    cat(paste("p-value adjustment method:", object$p.adj.method))
    cat("\n")
    cat("\n")
	}
