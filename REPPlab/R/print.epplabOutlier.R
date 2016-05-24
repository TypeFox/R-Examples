#' Print an epplabOutlier Object
#' 
#' Prints an \code{epplabOutlier} object.
#' 
#' The \code{print} function displays only the \code{outlier} matrix from the
#' \code{epplabOutlier} object.
#' 
#' @name print.epplabOutlier
#' @aliases print.epplabOutlier print,epplabOutlier-method
#' @docType methods
#' @param x Object of class \code{epplabOutlier}.
#' @param ... Additional parameters
#' @author Klaus Nordhausen
#' @keywords methods print
#' @examples
#' 
#' # creating data with 3 outliers
#' n <-300 
#' p <- 10
#' X <- matrix(rnorm(n*p),ncol=p)
#' X[1,1] <- 9
#' X[2,4] <- 7 
#' X[3,6] <- 8
#' # giving the data rownames, obs.1, obs.2 and obs.3 are the outliers.
#' rownames(X) <- paste("obs",1:n,sep=".")
#' 
#' PP<-EPPlab(X,PPalg="PSO",PPindex="KurtosisMax",n.simu=20, maxiter=20)
#' OUT<-EPPlabOutlier(PP, k = 3, location = median, scale = mad)
#' OUT
#' 
#' @export
print.epplabOutlier <- function(x, ...)
    {
    print(x$outlier,...)
    }
