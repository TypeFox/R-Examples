
#' Summarize Binomial Tests for csQCA Data
#'
#' Displays number of confirming cases and raw and adjusted p-scores 
#' following Binomial test of csQCA data.
#' @param object Object returned by \code{\link{csQCAbinTest}}.
#' @param ... Additional parameters to pass on.
#' @return Matrix of values for counterexamples and consistency.
#' @keywords binomial test csQCA
#' @export
#' @examples
#' test <- csQCAbinTest(freq.y=0.7, configs=list(aB=5, bCD=3, Ce=2),
#'    total.configs=20)
#' summary(test)


summary.csQCAbt <- function(object, ...){
	cat("Call:\n")
    print(object$call)
    cat("\n")
    cat("Counterexamples\n")
	printCoefmat(object$result[,1:3], P.values=TRUE, has.Pvalue=TRUE, tst.ind=1, signif.legend=FALSE)
    cat(paste("Total number of configurations:", object$total.configurations), "\n")
    cat(paste("p-value adjustment method:", object$p.adj.method))
    cat("\n")
    cat("\n")
	}
