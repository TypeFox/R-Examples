## ***************************************************************************
## Probability density function(pdf) of Generalized Power Weibull distribution
dgp.weibull <- function (x, alpha, theta, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- 1.0 + (x ^ alpha)
    pdf <- alpha * theta * (x^(alpha - 1.0))* (u^(theta - 1.0)) * exp(1.0 -(u ^ theta))  
    if(log)
        pdf <- log(pdf)
    return(pdf)  
}
## ***************************************************************************
## Cummulative distribution function(cdf) of Generalized Power Weibull distribution
pgp.weibull <- function (q, alpha, theta, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (q <= 0))    
        stop("Invalid arguments")
    u <- 1.0+(q ^ alpha)
    cdf <- 1.0 - exp(1.0 - (u ^ theta))    
    if(!lower.tail)
        cdf <- 1.0 - cdf
    if(log.p)
        cdf <- log(cdf)
    return(cdf)            
}
## ***************************************************************************
## Quantile function of Generalized Power Weibull distribution
qgp.weibull <- function (p, alpha, theta, lower.tail = TRUE, log.p = FALSE)
{
if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl <- ((1.0-log(1.0-p))^(1.0/theta)-1.0) ^ (1.0/alpha)    
    if (!lower.tail) 
      qtl <- ((1.0-log(p))^(1.0/theta)-1.0) ^ (1.0/alpha) 
    if (log.p)
        qtl <- log(qtl)    
    return(qtl)   
}
## ***************************************************************************
## Random variate generation from Generalized Power Weibull distribution
rgp.weibull <- function(n, alpha, theta)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (n <= 0))    
        stop("Invalid arguments") 
    return((((1.0-log(1.0-(runif(n))))^(1.0/theta))-1.0) ^ (1.0/alpha)) 
}
## *************************************************************************** 
## Reliability function of Generalized Power Weibull distribution
sgp.weibull <- function (x, alpha, theta)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (x <= 0))    
        stop("Invalid arguments")  
      u <- 1.0+(x ^ alpha)                
      return(exp(1.0-(u ^ theta)))   
}
## ***************************************************************************
## Hazard function of Generalized Power Weibull distribution
hgp.weibull <- function (x, alpha, theta)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (x <= 0))    
        stop("Invalid arguments")      
      u <- 1.0 + (x ^ alpha)          
      return(alpha * theta * (x^(alpha - 1.0))* (u ^(theta - 1.0)))   
} 
## ***************************************************************************
## Hazard rate average function of Generalized Power Weibull distribution
hra.gp.weibull <- function(x, alpha, theta)
{
    r <- sgp.weibull(x, alpha, theta)
    fra <-((-1)*log(r))/x
    return(fra)
}
## ***************************************************************************
## Conditional Hazard rate function of Generalized Power Weibull distribution
crf.gp.weibull <- function(x, t=0, alpha, theta)
{
    t <- t
    x <- x
    nume <- hgp.weibull(x+t, alpha, theta)
    deno <- hgp.weibull(x, alpha, theta)
    return(nume/deno)
  }
## ***************************************************************************
## Kolmogorov-Smirnov test (One-sample)for Generalized Power Weibull distribution
ks.gp.weibull <- function(x, alpha.est, theta.est, 
                      alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    theta <- theta.est
    res <- ks.test(x, pgp.weibull, alpha, theta, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pgp.weibull(t, alpha, theta)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}  
## ***************************************************************************
## Quantile-Quantile(QQ) plot for Generalized Power Weibull distribution
qq.gp.weibull <- function(x, alpha.est, theta.est, main=' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    theta <- theta.est       
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <-(k - 0.5)/n
    limx <- c(min(x), max(x))    
    Finv <- qgp.weibull(P,alpha, theta)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qgp.weibull(0.25, alpha, theta)
        y2 <- qgp.weibull(0.75, alpha, theta)
        m <- ((y2 - y1)/(x2 - x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))     
}
## ***************************************************************************
## Probability-Probability(PP) plot for Generalized Power Weibull distribution
pp.gp.weibull <- function(x, alpha.est, theta.est, main=' ',line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est
    theta <- theta.est
    F <- pgp.weibull(x,alpha, theta)
    Pemp <- sort(F)
    n <- length(x)
    k <- seq(1, n, by = 1)
    Pteo <-(k - 0.5)/n
    plot(Pemp, Pteo, xlab = xlab, ylab = ylab, col = 4, 
         xlim = c(0, 1), ylim = c(0, 1), main = main, lwd = 2, ...)
    if(line)
        lines(c(0, 1), c(0, 1), col = 2, lwd = 2)
    Cor.Coeff <- cor(Pemp, Pteo)
    Determination.Coeff <- (Cor.Coeff^2) * 100
    return(list(Cor.Coeff = Cor.Coeff, Determination.Coeff = Determination.Coeff))
}
## ***************************************************************************
## Akaike information criterium (AIC) and
## Bayesian information criterion (BIC) for Generalized Power Weibull distribution
abic.gp.weibull <- function(x, alpha.est, theta.est)
{ 
    alpha <- alpha.est
    theta <- theta.est 
    n <-length(x)
    p <-2
    f <- dgp.weibull(x,alpha, theta)
    l <-log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p 
    BIC <- - 2 * LogLik + p * log(n)                 
    return(list(LogLik=LogLik, AIC = AIC, BIC = BIC))
  } 
## ***************************************************************************
