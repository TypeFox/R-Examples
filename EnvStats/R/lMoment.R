lMoment <-
function (x, r = 1, method = "unbiased", plot.pos.cons = c(a = 0.35, 
    b = 0), na.rm = FALSE) 
{
    if (length(r) != 1 || !is.vector(r, mode = "numeric") || 
        r != trunc(r) || r < 1) 
        stop("'r' must be a positive integer")
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    wna <- which.na(x)
    if (length(wna)) {
        if (na.rm) 
            x <- x[-wna]
        else return(NA)
    }
    n <- length(x)
    if (n < 1) 
        stop("'x' must contain at least one non-missing, finite value")
    method <- match.arg(method, c("unbiased", "plotting.position"))
    x <- sort(x)
    index <- ifelse(r <= 4, r, 5)
    switch(index, mean(x), 2 * pwMoment(x, j = 1, method = method, 
        plot.pos.cons = plot.pos.cons) - mean(x), 6 * pwMoment(x, 
        j = 2, method = method, plot.pos.cons = plot.pos.cons) - 
        6 * pwMoment(x, j = 1, method = method, plot.pos.cons = plot.pos.cons) + 
        mean(x), 20 * pwMoment(x, j = 3, method = method, plot.pos.cons = plot.pos.cons) - 
        30 * pwMoment(x, j = 2, method = method, plot.pos.cons = plot.pos.cons) + 
        12 * pwMoment(x, j = 1, method = method, plot.pos.cons = plot.pos.cons) - 
        mean(x), {
        k <- 0:(r - 1)
        p.star <- (-1)^(r - 1 - k) * exp(lfactorial(r - 1 + k) - 
            2 * lfactorial(k) - lfactorial(r - 1 - k))
        b <- numeric(length(k))
        for (i in k) b[i + 1] <- pwMoment(x, j = i, method = method, 
            plot.pos.cons = plot.pos.cons)
        sum(p.star * b)
    })
}
