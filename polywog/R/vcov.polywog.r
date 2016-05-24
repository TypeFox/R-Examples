##' @S3method vcov polywog
vcov.polywog <- function(object, ...)
{
    ncf <- length(coef(object))
    if (!is.null(object$boot.matrix)) {
        ans <- var(t(as.matrix(object$boot.matrix)))
    } else {
        ans <- matrix(NA, nrow = ncf, ncol = ncf)
    }
    return(ans)
}
