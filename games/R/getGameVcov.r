## TODO: get rid of this function
##
## Calculates the variance-covariance matrix for a fitted model, including a
## procedure for catching the error (and returning a matrix of NAs) in case the
## Hessian is non-invertible.
##
getGameVcov <- function(hessian, fixed)
{
    hes <- hessian[!fixed, !fixed, drop = FALSE]
    vv <- tryCatch(solve(-hes), error = function(e) e)
    if (inherits(vv, "error")) {
        warning("variance-covariance matrix could not be calculated: ",
                vv$message)
        vv <- matrix(NA, nrow(hes), nrow(hes))
    }
    ans <- hessian
    ans[] <- NA
    ans[!fixed, !fixed] <- vv
    return(ans)
}
