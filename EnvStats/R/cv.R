cv <-
function (x, method = "moments", sd.method = "sqrt.unbiased", 
    l.moment.method = "unbiased", plot.pos.cons = c(a = 0.35, 
        b = 0), na.rm = FALSE) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    wna <- which.na(x)
    if (length(wna)) {
        if (na.rm) 
            x <- x[-wna]
        else return(NA)
    }
    method <- match.arg(method, c("moments", "l.moments"))
    if (method == "moments") {
        sd.method <- match.arg(sd.method, c("sqrt.unbiased", 
            "moments"))
        cv <- sd(x)/mean(x)
        if (sd.method != "sqrt.unbiased") {
            n <- length(x)
            cv <- sqrt((n - 1)/n) * cv
        }
    }
    else {
        l.moment.method <- match.arg(l.moment.method, c("unbiased", 
            "plotting.position"))
        cv <- lMoment(x, r = 2, method = l.moment.method, plot.pos.cons = plot.pos.cons)/mean(x)
    }
    cv
}
