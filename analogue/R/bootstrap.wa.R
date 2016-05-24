`bootstrap.wa` <- function(object, n.boot = 1000,
                           verbose = TRUE, ...) {
    retval <- predict(object, object$orig.x, CV = "bootstrap",
                      n.boot = n.boot, verbose = verbose)
    .call <- match.call()
    .call[[1]] <- as.name("bootstrap")
    retval$call = .call
    class(retval) <- "bootstrap.wa"
    return(retval)
}
