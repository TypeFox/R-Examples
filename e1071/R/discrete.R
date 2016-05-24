rdiscrete <- function (n, probs, values = 1:length(probs), ...)
{
    sample(values, size=n, replace=TRUE, prob=probs)
}



ddiscrete <- function (x, probs, values = 1:length(probs))
{
    
    if (length(probs) != length(values))
        stop("ddiscrete: probs and values must have the same length.")
    if (sum(probs < 0) > 0)
        stop("ddiscrete: probs must not contain negative values.")
    if (!is.array(x) && !is.vector(x) && !is.factor(x))
        stop("ddiscrete: x must be an array or a vector or a factor.")
    
    p <- probs/sum(probs)
    
    y <- as.vector(x)
    l <- length(y)
    z <- rep(0,l)
    
    for (i in 1:l)
        if (any(values == y[i]))
            z[i] <- p[values == y[i]]
    
    z <- as.numeric(z)
    if (is.array(x))
        dim(z) <- dim(x)
    
    return(z)
}


pdiscrete <- function (q, probs, values = 1:length(probs))
{
    
    if (length(probs) != length(values))
        stop("pdiscrete: probs and values must have the same length.")
    if (sum(probs < 0) > 0)
        stop("pdiscrete: probs must not contain negative values.")
    if (!is.array(q) & !is.vector(q))
        stop("pdiscrete: q must be an array or a vector")
    
    p <- probs/sum(probs)
    
    y <- as.vector(q)
    l <- length(y)
    z <- rep(0,l)
    
    for (i in 1:l)
        z[i] <- sum(p[values <= y[i]])
    
    z <- as.numeric(z)
    if (is.array(q))
        dim(z) <- dim(q)
    
    return(z)
}

qdiscrete <- function (p, probs, values = 1:length(probs))
{
    
    if (length(probs) != length(values))
        stop("qdiscrete: probs and values must have the same length.")
    if (sum(probs < 0) > 0)
        stop("qdiscrete: probs must not contain negative values.")
    if (!is.array(p) & !is.vector(p))
        stop("qdiscrete: p must be an array or a vector")
    
    probs <- cumsum(probs)/sum(probs)
    
    y <- as.vector(p)
    l <- length(y)
    z <- rep(0,l)
    
    for (i in 1:l)
        z[i] <- length(values) - sum(y[i] <= probs) + 1
    
    z <- as.numeric(z)
    z <- values[z]
    if (is.array(p))
        dim(z) <- dim(p)
    
    return(z)
  }



