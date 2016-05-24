is.present.svisit <-
function (object, n.imp = 0, ...)
{
    if (inherits(object, "svocc")) {
        rval <- fitted(object)
    } else {
        stop("not yet implemented for abundance models")
        z1 <- 1 - object$zif.probabilities
        p1 <- 1 - exp(-fitted(object))
        obs01 <- ifelse(object$y > 0, 1, 0)
        phi <- z1 * p1
        #delta <- object$detection.probabilities
        #Py1w0 <- (phi * (1 - delta))/((1 - phi) + (phi * (1 - delta)))
        rval <- pmax(obs01, phi)
        #rval <- pmax(obs01, Py1w0)
    }
    if (n.imp) {
        rval <- rbinom(length(object$y) * n.imp, 1, rval)
        dim(rval) <- c(length(object$y), n.imp)
        rownames(rval) <- case.names(object)
        colnames(rval) <- paste("imp", 1:n.imp, sep=".")
    } else names(rval) <- case.names(object)
    rval
}
