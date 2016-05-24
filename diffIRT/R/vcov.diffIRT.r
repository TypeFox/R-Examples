vcov.diffIRT = function (object, ...) {
    if (is.null(object$hessian))
        stop("to obtain the covariance matrix of the parameter estimates, you should re-fit the model using 'se = TRUE'.\n")
    covmat <- solve(object$hessian)
    covmat
}
