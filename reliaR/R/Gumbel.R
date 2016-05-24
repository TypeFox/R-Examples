## **************************************************************************
## Probability density function(pdf) of Gumbel distribution
dgumbel <- function(x, mu, sigma, log = FALSE)
{
    if((!is.numeric(mu)) || (!is.numeric(sigma)) || (!is.numeric(x)))
      stop("non-numeric argument to mathematical function")
    if (min(sigma) <= 0.0)  
      stop("sigma must be positive")
    u <- - (x - mu)/sigma
    pdf <- (1.0/sigma) * exp(u) * exp( - exp(u))    
    if (log)
         pdf<- log(pdf)
    return(pdf)   
}
## **************************************************************************
## Cummulative distribution function(cdf) of Gumbel distribution
pgumbel <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(mu)) || (!is.numeric(sigma)) || (!is.numeric(q)))
      stop("non-numeric argument to mathematical function")
    if (min(sigma) <= 0.0)  
      stop("sigma must be positive")
    u <- - (q - mu)/sigma
    cdf<- exp( - exp(u))    
    if(!lower.tail)
        cdf <- 1.0 - cdf
    if(log.p)
        cdf <- log(cdf)
    return(cdf)
}
## **************************************************************************
## Quantile function of Gumbel distribution
qgumbel <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(mu)) || (!is.numeric(sigma)) || (!is.numeric(p)))
      stop("non-numeric argument to mathematical function")
    if ((min(sigma) <= 0.0) || (p <= 0) || (p > 1.0)) 
      stop("Invalid arguments")
    qtl <- mu - sigma * log( - log(p))    
    if (!lower.tail) 
        qtl <- mu - sigma * log(- log(1.0 - p))  
    if (log.p) 
        qtl <- log(qtl)    
    return(qtl)   
}
## **************************************************************************
## Random variate generation from Gumbel distribution
rgumbel <- function(n, mu, sigma)
{
    if((!is.numeric(mu)) ||(!is.numeric(sigma))||(!is.numeric(n)))
      stop("non-numeric argument to mathematical function")
    if ((min(sigma)<= 0.0) || (n <= 0))  
    stop("Invalid arguments") 
    return(mu - sigma * log( - log( runif(n)))) 
}
## ************************************************************************** 
## Reliability function of Gumbel distribution
sgumbel <- function(x, mu, sigma)
{
    if((!is.numeric(mu)) ||(!is.numeric(sigma))||(!is.numeric(x)))
      stop("non-numeric argument to mathematical function")
    if (min(sigma) <= 0.0)  
      stop("sigma must be positive")  
    u <- - (x - mu)/sigma     
    return(1.0 - exp( - exp(u)))   
}
## **************************************************************************
## Hazard function of Gumbel distribution
hgumbel <- function(x, mu, sigma)
{
    if((!is.numeric(mu)) ||(!is.numeric(sigma))||(!is.numeric(x)))
      stop("non-numeric argument to mathematical function")
    if (min(sigma) <= 0.0)  
      stop("sigma must be positive")       
    u <- - (x - mu)/sigma
    nume <- (1.0/sigma) * exp(u) * exp( - exp(u))
    deno <- 1.0 - exp( - exp(u))     
    return(nume/deno)   
} 
## **************************************************************************
## Hazard rate average function of Gumbel distribution
hra.gumbel <- function(x, mu, sigma)
{
    r <- sgumbel(x, mu, sigma)
    fra <- ((-1)* log(r))/x
    return(fra)
}
## **************************************************************************
## Conditional Hazard rate function of Gumbel distribution
crf.gumbel <- function(x, t=0, mu, sigma)
{
    t <- t
    x <- x
    nume <- hgumbel(x+t, mu, sigma)
    deno <- hgumbel(x, mu, sigma)
    return(nume/deno)
}
## **************************************************************************
## Kolmogorov-Smirnov test (One-sample)for Gumbel distribution
ks.gumbel <- function(x, mu.est, sigma.est, 
               alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    mu <- mu.est
    sigma <- sigma.est 
    res <- ks.test(x, pgumbel, mu, sigma, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pgumbel(t, mu, sigma)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## **************************************************************************
## Quantile-Quantile(QQ) plot for Gumbel distribution
qq.gumbel <- function(x, mu.est, sigma.est, main=' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    mu <- mu.est
    sigma <- sigma.est
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k - 0.5) / n
    limx <- c(min(x), max(x))
    Finv <- qgumbel(P, mu, sigma)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qgumbel(0.25, mu, sigma)
        y2 <- qgumbel(0.75, mu, sigma)
        m <- ((y2-y1) / (x2-x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))
}
## **************************************************************************
## Probability-Probability(PP) plot for Gumbel distribution
pp.gumbel <- function(x, mu.est, sigma.est, main=' ', line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    mu <- mu.est
    sigma <- sigma.est
    F <- pgumbel(x, mu, sigma)
    Pemp <- sort(F)
    n <- length(x)
    k <- seq(1, n, by = 1)
    Pteo <-(k - 0.5) / n
    plot(Pemp, Pteo, xlab = xlab, ylab = ylab, col = 4, 
         xlim = c(0, 1), ylim = c(0, 1), main = main, lwd = 2, ...)
    if(line)
        lines(c(0, 1), c(0, 1), col = 2, lwd = 2)
    Cor.Coeff <- cor(Pemp, Pteo)
    Determination.Coeff <- (Cor.Coeff^2) * 100
    return(list(Cor.Coeff = Cor.Coeff, Determination.Coeff = Determination.Coeff))
}
## **************************************************************************
## Akaike information criterium (AIC)  and
## Bayesian information criterion (BIC) for Gumbel distribution
abic.gumbel <- function(x, mu.est, sigma.est)
{ 
    mu <- mu.est
    sigma <- sigma.est
    n <- length(x)
    p <- 2
    f <- dgumbel(x, mu, sigma)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p 
    BIC <- - 2 * LogLik + p * log(n)             
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## **************************************************************************
