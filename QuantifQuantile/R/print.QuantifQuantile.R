#' Print of QuantifQuantile results
#'
#' This function displays a small description of QuantifQuantile results.
#' 
#' @param x An object of class \code{QuantifQuantile}, which is the result of 
#' the \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} or 
#' \code{\link{QuantifQuantile.d}} functions.
#' @param \dots Not used.
#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimation through optimal quantization}, 
#' Journal of Statistical Planning and Inference, 2015 (156), 14-30.
#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimator based on optimal 
#' quantization: from theory to practice}, Submitted.

#' @seealso \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} and
#'  \code{\link{QuantifQuantile.d}}
#' @seealso \code{\link{plot.QuantifQuantile}}, 
#' \code{\link{summary.QuantifQuantile}}
#' @author Isabelle Charlier, Davy Paindaveine, Jerome Saracco
#' @examples
#' set.seed(644972)
#' n <- 300
#' X <- runif(300,-2,2)
#' Y <- X^2+rnorm(n)
#' res <- QuantifQuantile(X,Y,testN=seq(10,25,by=5))
#' print(res)
#' 
#' @export print.QuantifQuantile
#' @S3method print QuantifQuantile
print.QuantifQuantile <- function(x, ...) {
    stopifnot(class(x)=="QuantifQuantile")
    cat(paste("** QuantifQuantile results **", 
        "The result object is a list of length 9.", 
        "Most interesting components are the following:\n\n", 
        sep = "\n"))
    
    
    res <- array("", c(3, 2), list(1:3, c("name", "description")))
    res[1, ] <- c("$N_opt", 
                  "optimal value for N selected among testN by our criterion")
    res[2, ] <- c("$hatq_opt", 
                  "the estimated conditional quantiles obtained with N_opt")
    res[3, ] <- c("$hatISE_N", 
                  "estimated ISE as a function of N")
    
    print(res)
    
    cat(paste("\n For more information about QuantifQuantile outputs see the 
              QuantifQuantile help page\n"))
} 
