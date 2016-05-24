## Contains functions for the Pch (Piecewise constant hazards) distribution.

ppch <- function(q, cuts, levels, lower.tail = TRUE, log.p = FALSE){
    y <- Hpch(q, cuts, levels)
    if (log.p){
        if (lower.tail){
            y <- log(-expm1(-y)) ## log(-expm1(-y)) = log(1 - exp(-y))?!
        }else{
            y <- -y
        }
    }else{
        if (lower.tail){
            y <- -expm1(-y)
        }else{
            y <- exp(-y)
        }
    }
    y
}
        
dpch <- function(x, cuts, levels, log = FALSE){
    y <- hpch(x, cuts, levels) * ppch(x, cuts, levels, lower.tail = FALSE)
    if (log) y <- log(y)
    y
}

hpch <- function(x, cuts, levels, log = FALSE){
    cuts <- sort(unique(cuts))
    p <- length(levels)
    if (length(cuts) != (p - 1))stop("Must be one more level than cut.")
    if (any(cuts <= 0)) stop("all cuts must be positive")
    if (any(x < 0)) stop("x must be all positive.")
    y <- numeric(length(x))
    cuts <- c(0, cuts, Inf)
    y[(cuts[1] <= x) & (x <= cuts[1 + 1])] <- levels[1]
    if (p > 1.5){
        for (i in 2:p){
            y[(cuts[i] < x) & (x <= cuts[i + 1])] <- levels[i]
        }
    }
    if (log) y <- log(y)
    y
}
    
Hpch <- function(x, cuts, levels, log.p = FALSE){
    cuts <- sort(unique(cuts))
    p <- length(levels)
    if (length(cuts) != (p - 1))stop("Must be one more level than cut.")
    if (any(cuts <= 0)) stop("all cuts must be positive")
    if (any(x < 0)) stop("x must be all positive.")
    y <- numeric(length(x))
    cuts <- c(0, cuts, Inf)
    who <- (cuts[1] <= x) & (x <= cuts[1 + 1])
    if (sum(who)){
        z <- x[who]
        y[who] <- levels[1] * z
    }
    su <- levels[1] * cuts[2]
    if (p > 1.5){
        for (i in 2:p){
            who <- (cuts[i] < x) & (x <= cuts[i + 1])
            if (sum(who)){
                y[who] <- su + levels[i] * (x[who] - cuts[i])
            }
            su <- su + levels[i] * (cuts[i + 1] - cuts[i])
        }
    }
    if (log.p) y <- log(y)
    y
}
            
qpch <- function(p, cuts, levels, lower.tail = TRUE, log.p = FALSE){
    if (log.p) p <- exp(p)
    if (any(p >= 1)) stop("p must be < 1") 
    if (any(p <= 0)) stop("p must be > 0") 
    if (!lower.tail) p <- 1 - p
    n <- length(p)
    y <- numeric(n)
    f <- function(q, x){
        ppch(q, cuts, levels) - x
    }
    for (i in 1:n){
        y[i] <- uniroot(f, interval = c(0, 2000), x = p[i])$root
    }
    y
}


rpch <- function(n, cuts, levels){
    x <- runif(n)
    qpch(x)
}
