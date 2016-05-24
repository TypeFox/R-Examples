`densityplot.residLen` <- function(x, ..., xlab = NULL,
                                   ylab = NULL) {
    ## produce the data object
    n.train <- length(x$train)
    n.pass <- length(x$passive)
    dat <- data.frame(lengths = c(x$train, x$passive),
                      type = rep(c("Training Set", "Passive"),
                      times = c(n.train, n.pass)))
    ## labels
    if(is.null(xlab)) {
        if(attr(x, "method") == "cca") {
            xlab <- expression(paste(Squared ~ chi^2 ~
                residual ~ distance))
        } else {
            xlab <- "Squared Euclidean residual distance"
        }
    }
    if(is.null(ylab)) {
        ylab <- "Density"
    }
    ## plotting
    densityplot(~ lengths | type, data = dat, layout = c(1,2),
                n = 100, xlab = xlab, from = 0, ...)
}
