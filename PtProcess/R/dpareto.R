dpareto <- function(x, lambda, a, log=FALSE){
    if (any(x < a)) stop("Invalid data: x < a")
    if (any(lambda <= 0)) stop("All values of lambda must be positive")
    if (length(lambda)!=1 & length(lambda)!=length(x))
        stop("Argument lambda is wrong length")
    if (a <= 0) stop("Argument a must be positive")
    d <- lambda/a * (a/x)^(lambda+1)
    if (log==TRUE) d <- log(d)
    return(d)
}


dtappareto <- function(x, lambda, theta, a, log=FALSE){
    if (any(x < a)) stop("Invalid data: x < a")
    if (any(lambda <= 0)) stop("All values of lambda must be positive")
    if (length(lambda)!=1 & length(lambda)!=length(x))
        stop("Argument lambda is wrong length")
    if (any(theta <= 0)) stop("All values of theta must be positive")
    if (length(theta)!=1 & length(theta)!=length(x))
        stop("Argument theta is wrong length")
    if (a <= 0) stop("Argument a must be positive")
    d <- (lambda/x + 1/theta) * (a/x)^lambda * exp((a-x)/theta)
    if (log==TRUE) d <- log(d)
    return(d)
}


ppareto <- function(q, lambda, a, lower.tail=TRUE, log.p=FALSE){
    if (any(q < a)) stop("Invalid request: q < a")
    if (any(lambda <= 0)) stop("All values of lambda must be positive")
    if (length(lambda)!=1 & length(lambda)!=length(q))
        stop("Argument lambda is wrong length")
    if (a <= 0) stop("Argument a must be positive")
    p <- 1 - (a/q)^lambda
    if (lower.tail==FALSE) p <- 1-p
    if (log.p==TRUE) p <- log(p)
    return(p)
}


ptappareto <- function(q, lambda, theta, a, lower.tail=TRUE, log.p=FALSE){
    if (any(q < a)) stop("Invalid request: q < a")
    if (any(theta <= 0)) stop("All values of theta must be positive")
    if (length(theta)!=1 & length(theta)!=length(q))
        stop("Argument theta is wrong length")
    if (any(lambda <= 0)) stop("All values of lambda must be positive")
    if (length(lambda)!=1 & length(lambda)!=length(q))
        stop("Argument lambda is wrong length")
    if (a <= 0) stop("Argument a must be positive")
    p <- 1 - (a/q)^lambda *exp((a-q)/theta)
    if (lower.tail==FALSE) p <- 1-p
    if (log.p==TRUE) p <- log(p)
    return(p)
}


qpareto <- function(p, lambda, a, lower.tail=TRUE, log.p=FALSE){
    if (any(lambda <= 0)) stop("All values of lambda must be positive")
    if (length(lambda)!=1 & length(lambda)!=length(p))
        stop("Argument lambda is wrong length")
    if (a <= 0) stop("Argument a must be positive")
    if (log.p==TRUE) p <- exp(p)
    if (lower.tail==FALSE) p <- 1-p
    q <- a*(1-p)^(-1/lambda)
    return(q)
}


qtappareto <- function(p, lambda, theta, a, lower.tail=TRUE, log.p=FALSE, tol=1e-8){
    if (any(theta <= 0)) stop("All values of theta must be positive")
    if (length(theta)!=1 & length(theta)!=length(p))
        stop("Argument theta is wrong length")
    if (any(lambda <= 0)) stop("All values of lambda must be positive")
    if (length(lambda)!=1 & length(lambda)!=length(p))
        stop("Argument lambda is wrong length")
    if (a <= 0) stop("Argument a must be positive")
    if (log.p==TRUE) p <- exp(p)
    if (lower.tail==FALSE) p <- 1-p
    q <- rep(a+1, length(p))
    #   solve using Newton-Raphson method
    repeat{
        delta <- (ptappareto(q, lambda, theta, a) - p)/
                     dtappareto(q, lambda, theta, a)
        q <- pmax(q - delta, a)
        if(max(abs(delta)) < tol) break
    }
    return(q)
}


rpareto <- function(n, lambda, a){
    if (any(lambda <= 0)) stop("All values of lambda must be positive")
    if (length(lambda)!=1 & length(lambda)!=n)
        stop("Argument lambda is wrong length")
    if (a <= 0) stop("Argument a must be positive")
    return(a*exp(rexp(n, rate=lambda)))
}


rtappareto <- function(n, lambda, theta, a){
    if (any(theta <= 0)) stop("All values of theta must be positive")
    if (length(theta)!=1 & length(theta)!=n)
        stop("Argument theta is wrong length")
    if (any(lambda <= 0)) stop("All values of lambda must be positive")
    if (length(lambda)!=1 & length(lambda)!=n)
        stop("Argument lambda is wrong length")
    if (a <= 0) stop("Argument a must be positive")
    return(pmin(rpareto(n, lambda, a), (a + rexp(n, 1/theta))))
}


ltappareto <- function(data, lambda, theta, a){
    if (length(lambda)!=1 | length(theta)!=1 | length(a)!=1)
        stop("Invalid arguments")
    #   log-likelihood for tapered Pareto
    LL <- sum(log(dtappareto(data, lambda, theta, a)))
    #    calculate derivatives
    v <- lambda*theta + data
    dLL.lambda <- theta*sum(1/v) - sum(log(data/a))
    dLL.theta <- -sum(data/v)/theta - sum(a-data)/theta^2
    attr(LL, "gradient") <- c(dLL.lambda, dLL.theta)
    #    calculate hessian
    h11 <- -theta^2*sum(1/v^2)
    h12 <- sum(data/v^2)
    h22 <- 1/theta^2*sum(data*(2*lambda*theta+data)/v^2) +
           2/theta^3*sum(a-data)
    attr(LL, "hessian") <- matrix(c(h11, h12, h12, h22), nrow=2)
    return(LL)
}


