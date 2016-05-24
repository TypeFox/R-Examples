pwMoment <-
function (x, j = 0, k = 0, method = "unbiased", plot.pos.cons = c(a = 0.35, 
    b = 0), na.rm = FALSE) 
{
    if (length(j) != 1 || !is.vector(j, mode = "numeric") || 
        j != trunc(j) || j < 0 || length(k) != 1 || !is.vector(k, 
        mode = "numeric") || k != trunc(k) || k < 0) 
        stop("'j' and 'k' must be non-negative integers")
    if (j > 0 && k > 0) 
        stop("Either 'j' or 'k' (or both) must be 0")
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
    if (j == 0 && k == 0) {
        pwm <- mean(x)
    }
    else {
        x <- sort(x)
        method <- match.arg(method, c("unbiased", "plotting.position"))
        if (method == "plotting.position") {
            if (!is.vector(plot.pos.cons, mode = "numeric") || 
                length(plot.pos.cons) != 2) 
                stop("'plot.pos.cons' must be a numeric vector of length 2")
            if (any(is.na(match(c("a", "b"), names(plot.pos.cons))))) 
                names(plot.pos.cons) <- c("a", "b")
            p <- ((1:n) - plot.pos.cons["a"])/(n + plot.pos.cons["b"])
        }
        if (j > 0) {
            if (j > (n - 1)) 
                stop("'j' must be between 0 and 'n-1'")
            pwm <- switch(method, unbiased = {
                index <- (j + 1):n
                sum((x[index]/n) * exp(lchoose(index - 1, j) - 
                  lchoose(n - 1, j)))
            }, plotting.position = mean(p^j * x))
        }
        else {
            if (k > (n - 1)) 
                stop("'k' must be between 0 and 'n-1'")
            pwm <- switch(method, unbiased = {
                index <- 1:(n - k)
                sum((x[index]/n) * exp(lchoose(n - index, k) - 
                  lchoose(n - 1, k)))
            }, plotting.position = mean((1 - p)^k * x))
        }
    }
    pwm
}
