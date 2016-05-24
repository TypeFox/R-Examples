kurtosis <-
function (x, na.rm = FALSE, method = "fisher", l.moment.method = "unbiased", 
    plot.pos.cons = c(a = 0.35, b = 0), excess = TRUE) 
{
    if (!(is.vector(x, mode = "numeric") && is.numeric(x))) {
        warning("'x' is not a numeric vector")
        return(NA)
    }
    wna <- which.na(x)
    if (length(wna)) {
        if (na.rm) 
            x <- x[-wna]
        else return(NA)
    }
    n <- length(x)
    method <- match.arg(method, c("fisher", "moment", "l.moment"))
    if (method == "fisher") {
        if (n < 4 || length(unique(x)) == 1) {
            warning(paste("When method='fisher',", "'x' must have at least 4", 
                "non-missing values and at least 2 must", "be distinct.\n"))
            return(NA)
        }
    }
    else if (n < 3 || length(unique(x)) == 1) {
        warning(paste("When method='", method, "',", "'x' must have at least 3", 
            "non-missing values and at least 2 must", "be distinct.\n"))
        return(NA)
    }
    kurtosis <- switch(method, moment = {
        x <- x - mean(x)
        (sum(x^4)/n)/(sum(x^2)/n)^2 - 3
    }, fisher = {
        x <- x - mean(x)
        ((n + 1) * (n - 1) * ((sum(x^4)/n)/(sum(x^2)/n)^2 - (3 * 
            (n - 1))/(n + 1)))/((n - 2) * (n - 3))
    }, l.moment = {
        l.moment.method <- match.arg(l.moment.method, c("unbiased", 
            "plotting.position"))
        lMoment(x, r = 4, method = l.moment.method, plot.pos.cons = plot.pos.cons)/lMoment(x, 
            r = 2, method = l.moment.method, plot.pos.cons = plot.pos.cons) - 
            3
    })
    if (!excess) 
        kurtosis <- kurtosis + 3
    kurtosis
}
