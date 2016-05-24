`histogram.residLen` <- function(x, ..., xlab = NULL, ylab = NULL,
                                 type = c("percent","count","density")) {
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
    if(missing(type))
        type <- "count"
    type <- match.arg(type)
    ylab <- switch(type,
                   count = "Frequency",
                   percent = "%",
                   density = "Density")
    ## plotting
    histogram(~ lengths | type, data = dat, layout = c(1,2),
              xlab = xlab, ylab = ylab, ...)
}
