skewness <-
function (x, na.rm = FALSE, method = "fisher", l.moment.method = "unbiased", 
    plot.pos.cons = c(a = 0.35, b = 0)) 
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
        if (n < 3 || length(unique(x)) == 1) {
            warning(paste("When method='fisher',", "'x' must have at least 3", 
                "non-missing values and at least 2 must", "be distinct.\n"))
            return(NA)
        }
    }
    else if (n < 2 || length(unique(x)) == 1) {
        warning(paste("When method='", method, "',", " 'x' must have at least 2", 
            " non-missing distinct values.\n", sep = ""))
        return(NA)
    }
    skew <- switch(method, moment = {
        x <- x - mean(x)
        (sum(x^3)/n)/(sum(x^2)/n)^1.5
    }, fisher = {
        x <- x - mean(x)
        ((sqrt(n * (n - 1))/(n - 2)) * (sum(x^3)/n))/((sum(x^2)/n)^1.5)
    }, l.moment = {
        l.moment.method <- match.arg(l.moment.method, c("unbiased", 
            "plotting.position"))
        lMoment(x, r = 3, method = l.moment.method, plot.pos.cons = plot.pos.cons)/lMoment(x, 
            r = 2, method = l.moment.method, plot.pos.cons = plot.pos.cons)
    })
    skew
}
