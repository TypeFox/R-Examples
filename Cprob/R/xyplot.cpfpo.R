xyplot.cpfpo <- function(x, data = NULL, conf.int = TRUE, level = 0.95,
                         odds = TRUE, intercept = TRUE, ylab, xlab,
                         lty = c(1,3,3), col = c(1,1,1), ...) {
    
    if (!inherits(x, "cpfpo")) {
        stop("'x' must be of class 'cpfpo'")
    }
    
    if (missing(ylab)) {
        if (odds==TRUE)
            ylab <- "odds-ratio"
        else
            ylab <- "log-odds-ratio"
    }
    
    if (missing(xlab))
        xlab <- "Time"
    
    ncov <- if (intercept) 1:dim(x$alpha)[2] else 2:dim(x$alpha)[2]
    dat <- lapply(ncov, function(i) {
        temp <- data.frame(coef = x$alpha[, i],
                           var = x$valpha[, i],
                           time = x$tis,
                           cov = colnames(x$alpha)[i])
        temp
    })
    dat <- do.call(rbind, dat)
    
    if (conf.int) {
        z <- qnorm(0.975)
        dat$lower <- with(dat, coef - z * sqrt(var))
        dat$upper <- with(dat, coef + z * sqrt(var))
        if (odds) {
            dat$upper <- exp(dat$upper)
            dat$lower <- exp(dat$lower)
            dat$coef <- exp(dat$coef)
        }
        aa <- xyplot(coef + lower + upper ~ time | cov, dat, type = "s",
                     col = col, lty = lty, xlab = xlab, ylab = ylab, ...)
    }
    else {
        if (odds) {
            with(dat, coef <- exp(coef))
        }
        aa <- xyplot(coef ~ time | cov, dat, type = "s",
                     col = col, lty = lty, xlab = xlab, ylab = ylab, ...)
    }
    
    aa
    
}
