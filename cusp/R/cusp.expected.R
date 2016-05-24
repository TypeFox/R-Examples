`cusp.expected` <-
function (alpha, beta, y, method = c("delay", "maxwell", "expected")) 
{
    method <- match.arg(method[1], c("maxwell", "delay", "expected"))
    if (method %in% c("maxwell", "delay")) 
        modi <- t(Vectorize(cusp.extrema)(alpha, beta))
    if (missing(y) && method == "delay") 
        stop("Observations y have to be provided for method=\"delay\"")
    val <- switch(method, expected = {
        .v <- t(Vectorize(function(alpha, beta) {
            val <- Vectorize(cusp.nc)(alpha, beta, 0:2)
            val <- val/c(1, val[0], val[0] * 2)
            val
        })(alpha, beta))
        colnames(.v) = c("norm.const", "E(Y)", "E(Y^2)/2")
        .v[, 2:3]
    }, delay = {
        .v <- ifelse(y < modi[, 2], modi[, 1], modi[, 3])
        .v <- as.matrix(.v)
        colnames(.v) <- "Delay"
        .v
    }, maxwell = {
        .v <- ifelse(alpha < 0, modi[, 1], modi[, 3])
        .v <- as.matrix(.v)
        colnames(.v) <- "Maxwell"
        .v
    })
    val
}

