## Methods to model and forecast the amount that members are spending during
## transactions.

library(hypergeo)
library(lattice)

spend.marginal.likelihood <- function(params, m.x, x) {
    
    max.length <- max(length(m.x), length(x))
    
    if (max.length%%length(m.x)) 
        warning("Maximum vector length not a multiple of the length of m.x")
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    
    dc.check.model.params(c("p", "q", "gamma"), params, "spend.marginal.likelihood")
    
    if (any(m.x < 0) || !is.numeric(m.x)) 
        stop("m.x must be numeric and may not contain negative numbers.")
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    
    if (any(x == 0) || any(m.x == 0)) {
        warning("Customers with 0 transactions or 0 average spend in spend.marginal.likelihood")
    }
    
    m.x <- rep(m.x, length.out = max.length)
    x <- rep(x, length.out = max.length)
    
    p <- params[1]
    q <- params[2]
    gamma <- params[3]
    
    result <- rep(0, max.length)
    
    ## non.zero: a vector indicating which elements have neither x == 0 or m.x == 0
    non.zero <- which(x > 0 & m.x > 0)
    
    result[non.zero] <- exp(lgamma(p * x[non.zero] + q) - lgamma(p * x[non.zero]) - 
        lgamma(q) + q * log(gamma) + (p * x[non.zero] - 1) * log(m.x[non.zero]) + 
        (p * x[non.zero]) * log(x[non.zero]) - (p * x[non.zero] + q) * log(gamma + 
        m.x[non.zero] * x[non.zero]))
    
    return(result)
}

spend.LL <- function(params, m.x, x) {
    
    max.length <- max(length(m.x), length(x))
    
    if (max.length%%length(m.x)) 
        warning("Maximum vector length not a multiple of the length of m.x")
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    
    dc.check.model.params(c("p", "q", "gamma"), params, "spend.LL")
    
    if (any(m.x < 0) || !is.numeric(m.x)) 
        stop("m.x must be numeric and may not contain negative numbers.")
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    
    if (any(x == 0) || any(m.x == 0)) {
        warning("Customers with 0 transactions or 0 average spend in spend.LL")
    }
    
    m.x <- rep(m.x, length.out = max.length)
    x <- rep(x, length.out = max.length)
    
    p <- params[1]
    q <- params[2]
    gamma <- params[3]
    
    ll <- rep(0, max.length)
    
    ## non.zero: a vector indicating which elements have neither x == 0 or m.x == 0
    non.zero <- which(x > 0 & m.x > 0)
    
    p <- params[1]
    q <- params[2]
    gamma <- params[3]
    
    ll[non.zero] <- (-lbeta(p * x[non.zero], q) + q * log(gamma) + (p * x[non.zero] - 
        1) * log(m.x[non.zero]) + (p * x[non.zero]) * log(x[non.zero]) - (p * x[non.zero] + 
        q) * log(gamma + m.x[non.zero] * x[non.zero]))
    
    return(ll)
}


spend.EstimateParameters <- function(m.x.vector, x.vector, par.start = c(1, 1, 1), 
    max.param.value = 10000) {
    
    if (any(m.x.vector < 0) || !is.numeric(m.x.vector)) 
        stop("m.x must be numeric and may not contain negative numbers.")
    if (any(x.vector < 0) || !is.numeric(x.vector)) 
        stop("x must be numeric and may not contain negative numbers.")
    
    if (any(x.vector == 0) || any(m.x.vector == 0)) {
        warning("Customers with 0 transactions or 0 average spend in spend.LL")
    }
    
    if (length(m.x.vector) != length(x.vector)) {
        stop("m.x.vector and x.vector must be the same length.")
    }
    
    spend.eLL <- function(params, m.x.vector, x.vector, max.param.value) {
        params <- exp(params)
        params[params > max.param.value] <- max.param.value
        return(-1 * sum(spend.LL(params, m.x.vector, x.vector)))
    }
    logparams <- log(par.start)
    results <- optim(logparams, spend.eLL, m.x.vector = m.x.vector, x.vector = x.vector, 
        max.param.value = max.param.value, method = "L-BFGS-B")
    estimated.params <- exp(results$par)
    estimated.params[estimated.params > max.param.value] <- max.param.value
    return(estimated.params)
}

spend.expected.value <- function(params, m.x, x) {
    
    max.length <- max(length(m.x), length(x))
    
    if (max.length%%length(m.x)) 
        warning("Maximum vector length not a multiple of the length of m.x")
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    
    dc.check.model.params(c("p", "q", "gamma"), params, "spend.expected.value")
    
    if (any(m.x < 0) || !is.numeric(m.x)) 
        stop("m.x must be numeric and may not contain negative numbers.")
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    
    m.x <- rep(m.x, length.out = max.length)
    x <- rep(x, length.out = max.length)
    
    p <- params[1]
    q <- params[2]
    gamma <- params[3]
    M <- (gamma + m.x * x) * p/(p * x + q - 1)
    return(M)
}


spend.plot.average.transaction.value <- function(params, m.x.vector, x.vector, xlab = "Average Transaction Value", 
    ylab = "Marginal Distribution of Average Transaction Value", title = "Actual vs. Expected Average Transaction Value Across Customers") {
    
    if (any(m.x.vector < 0) || !is.numeric(m.x.vector)) 
        stop("m.x must be numeric and may not contain negative numbers.")
    if (any(x.vector < 0) || !is.numeric(x.vector)) 
        stop("x must be numeric and may not contain negative numbers.")
    
    if (any(x.vector == 0) || any(m.x.vector == 0)) {
        warning("Customers with 0 transactions or 0 average spend in spend.plot.average.transaction.value have been removed before plotting.")
    }
    
    if (length(m.x.vector) != length(x.vector)) {
        stop("m.x.vector and x.vector must be the same length.")
    }
    
    # remove any customers with zero repeat transactions
    ave.spending <- m.x.vector[which(x.vector > 0)]
    tot.transactions <- x.vector[which(x.vector > 0)]
    
    f.m.x <- spend.marginal.likelihood(params, ave.spending, tot.transactions)
    plot(ave.spending, y = f.m.x, pch = 16, type = "n", xlab = xlab, ylab = ylab, 
        main = title)
    lines(density(ave.spending, bw = "nrd", adjust = 0.6), col = 1, lty = 1)
    lines(smooth.spline(ave.spending, y = f.m.x, w = tot.transactions, df = 15), 
        col = 2, lty = 2)
    legend("topright", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
    return(f.m.x)
}
 
