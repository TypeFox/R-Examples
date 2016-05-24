#' Print an ezsim Object. See \code{\link{ezsim}} for details.
#' @name print.ezsim
#' @aliases print.ezsim
#' @title Print an ezsim Object.
#' @method print ezsim
#' @param x An ezsim Object
#' @param \dots unused
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \code{\link{ezsim}}
#' @export  
print.ezsim <-
function(x,...){
    cat("Number of Simulation :",x$m ,"\n")
    cat("\n")
    cat("Estimator : \n")
    print(x$estimator)
    cat("\n")
    cat("dgp : \n")
    print(x$dgp)
    cat("\n")
    if (is.function(x$true_value)){
        cat("True value of estimator : \n")
        print(x$true_value)
        cat("\n")
    }
    print(x$parameter_def)
}
