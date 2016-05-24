## *************************************************************************
## Probability density function(pdf) of Weibull Extension(WE)  distribution
dweibull.ext <- function (x, alpha, beta, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))    
        stop("Invalid arguments")
     u <- exp((x/alpha) ^ beta)
    pdf <-exp(log(beta)+(beta -1.0)*(log(x)-log(alpha))+log(u) + alpha*(1.0 - u))  
    if(log)
        pdf <- log(pdf)
    return(pdf)
}
## *************************************************************************
## Cummulative distribution function(cdf) of Weibull Extension(WE)  distribution
pweibull.ext <- function (q, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (q <= 0))    
        stop("Invalid arguments")
     u <- exp((q/alpha) ^ beta)
    cdf <- 1.0 - exp(alpha * (1.0 - u))       
    if(!lower.tail)
        cdf <- 1.0 - cdf
    if(log.p)
        cdf <- log(cdf)
    return(cdf)
}
## *************************************************************************
## Quantile function of Weibull Extension(WE)distribution
qweibull.ext <- function (p, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl <-  alpha *(log(1.0 -  (log(1.0 - p)/(alpha))) ^ (1.0/beta))    
    if (!lower.tail) 
       qtl <- alpha*((log(1.0 -  (log(p)/alpha))) ^ (1.0/beta))   
    if(log.p)
        qtl <- log(qtl)
    return(qtl)
}
## *************************************************************************
## Random variate generation from Weibull Extension (WE) distribution
rweibull.ext <- function(n, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (n <= 0))
        stop("Invalid arguments")
    return(alpha*((log(1.0 -  (log(1.0 - runif(n))/alpha))) ^ (1.0/beta)))
}
## *************************************************************************
## Reliability function of Weibull Extension(WE)  distribution
sweibull.ext <- function (x, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))
        stop("Invalid arguments")    
     u <- exp((x/alpha) ^ beta)   
    return(exp(alpha * (1.0 - u)))   
}
## *************************************************************************
## Hazard function of Weibull Extension (WE) distribution
hweibull.ext <- function (x, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))
        stop("Invalid arguments")         
    u <- exp((x/alpha) ^ beta)
    hazard <- exp(log(beta) + (beta - 1.0) * (log(x)-log(alpha))+log(u))      
    return(hazard)   
} 
## *************************************************************************
## Hazard rate average function of Weibull Extension(WE)  distribution
hra.weibull.ext <-function(x, alpha, beta)
{
  r <- sweibull.ext(x, alpha, beta)
  fra <- ((-1) * log(r))/x
  return(fra)
}
## *************************************************************************
## Conditional Hazard rate function of Weibull Extension(WE)  distribution
crf.weibull.ext <-function(x, t=0, alpha, beta)
{
    t <- t
    x <- x
    nume <- hweibull.ext(x+t, alpha, beta)
    deno <- hweibull.ext(x, alpha, beta)
    return(nume/deno)
}
## *************************************************************************
## Kolmogorov-Smirnov test (One-sample) for Weibull Extension (WE) distribution
ks.weibull.ext <- function(x, alpha.est, beta.est, 
       alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    beta <- beta.est
    res <-ks.test(x, pweibull.ext, alpha, beta, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pweibull.ext(t, alpha, beta)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## *************************************************************************
## Quantile-Quantile(QQ) plot for Weibull Extension(WE)distribution
qq.weibull.ext <- function(x, alpha.est, beta.est, main=' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    beta <- beta.est
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k - 0.5)/n   
    limx <- c(min(x),  max(x))
    Finv <- qweibull.ext(P, alpha, beta)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qweibull.ext(0.25, alpha, beta)
        y2 <- qweibull.ext(0.75, alpha, beta)
        m <- ((y2 - y1) / (x2 - x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))
}
## *************************************************************************
## Probability-Probability (PP) plot for Weibull Extension (WE) distribution
pp.weibull.ext <- function(x, alpha.est, beta.est, main=' ', line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est
    beta <- beta.est
    F <- pweibull.ext(x, alpha, beta)
    Pemp <- sort(F)
    n <- length(x)
    k <- seq(1, n, by = 1)
    Pteo <- (k - 0.5) / n
    plot(Pemp, Pteo, xlab = xlab, ylab = ylab, col = 4, 
         xlim = c(0, 1), ylim = c(0, 1), main = main, lwd = 2, ...)
    if(line)
        lines(c(0, 1), c(0, 1), col = 2, lwd = 2)
    Cor.Coeff <- cor(Pemp, Pteo)
    Determination.Coeff <- (Cor.Coeff^2) * 100
    return(list(Cor.Coeff = Cor.Coeff, Determination.Coeff = Determination.Coeff))
}
## *************************************************************************
## Akaike information criterium (AIC) and
## Bayesian information criterion (BIC) for Weibull Extension distribution
abic.weibull.ext <- function(x, alpha.est, beta.est)
{ 
    alpha <- alpha.est
    beta <- beta.est
    n <- length(x)
    p <- 2
    f <- dweibull.ext(x, alpha, beta)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p  
    BIC <- - 2 * LogLik + p * log(n)   
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## *************************************************************************
