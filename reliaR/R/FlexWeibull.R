## ****************************************************************************
## Probability density function(pdf) of Flexible Weibull distribution
dflex.weibull <- function (x, alpha, beta, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp((alpha * x) - (beta / x))
    pdf <-(alpha + beta /(x * x)) * u * exp(-u)
    if (log) 
        pdf<- log(pdf)
    return(pdf)
}
## ****************************************************************************
## Cummulative distribution function(cdf) of flexible Weibull distribution
pflex.weibull <- function (q, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (q <= 0))    
        stop("Invalid arguments")
    u <- exp((alpha * q) - (beta / q))
    cdf <- 1.0 - exp(-u)
    if (!lower.tail) 
        cdf<- 1.0 - cdf
    if (log.p) 
        cdf<- log(cdf)
    return(cdf)
}
## ****************************************************************************
## Quantile function of flexible Weibull distribution 
qflex.weibull <- function(p, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
    if(!is.numeric(alpha) || !is.numeric(beta))
    {stop("non-numeric argument to mathematical function")}
    if (min(alpha) <=0 || min(beta) <=0 || p <= 0.0) 
    stop("Invalid arguments") 
    tmp <- log(-log(1.0 - p))
    if (!lower.tail) 
        tmp <- log(-log(p))
        qtl <- (1.0/(2 * alpha))* (tmp + (tmp^2.0 + (4.0 * alpha * beta))^0.5)
    if (log.p) 
        qtl <- log(qtl)    
    return(qtl)   
}   
## ****************************************************************************
## Random variate generation from flexible Weibull distribution
rflex.weibull <- function(n, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (n <= 0))    
        stop("Invalid arguments")
    u <- runif(n)
    tmp <- log(-log(1 - u)) 
    return((1.0/(2 * alpha)) * (tmp + (tmp^2.0 + (4.0 * alpha * beta)) ^0.5))  
}
## ****************************************************************************
## Reliability function of Burr distribution
sflex.weibull <- function (x, alpha, beta)
  {
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    return(exp(-exp((alpha * x) - (beta / x))))  
  }
## ****************************************************************************
## Hazard function of flexible Weibull distribution
hflex.weibull <- function (x, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))    
        stop("Invalid arguments") 
    u <- exp((alpha * x) - (beta / x))
    return((alpha + (beta/(x * x))) * u )  
}
## ****************************************************************************
## Hazard rate average function of flexible Weibull distribution
hra.flex.weibull <- function(x, alpha, beta)
{
    r <- sflex.weibull(x, alpha, beta)
    fra <-((-1) * log(r)) / x
    return(fra)
}
## ****************************************************************************
## Conditional Hazard rate function of flexible Weibull distribution
crf.flex.weibull <- function(x, t=0, alpha, beta)
{
    t <- t
    x <- x
    nume <- hflex.weibull(x+t, alpha, beta)
    deno <- hflex.weibull(x, alpha, beta)
    return(nume/deno)
}
## ****************************************************************************
## Kolmogorov-Smirnov test (One-sample)for flexible Weibull distribution
ks.flex.weibull <- function(x, alpha.est, beta.est, 
                   alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    beta <- beta.est
    res <- ks.test(x, pflex.weibull, alpha, beta, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pflex.weibull(t, alpha, beta)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## ****************************************************************************
## Quantile-Quantile(QQ) plot for flexible Weibull distribution
qq.flex.weibull <- function(x, alpha.est, beta.est, main=' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est 
    beta <- beta.est        
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k - 0.5)/n
    limx <- c(min(x), max(x))
    Finv <- qflex.weibull(P, alpha, beta)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qflex.weibull(0.25, alpha, beta)
        y2 <- qflex.weibull(0.75, alpha, beta)
        m <- ((y2 - y1)/(x2 - x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 4,lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))            
}
## ****************************************************************************
## Probability-Probability(PP) plot for flexible Weibull distribution
pp.flex.weibull <- function(x, alpha.est, beta.est, main=' ', line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est 
    beta <- beta.est              
    F <- pflex.weibull(x,alpha, beta)
    Pemp <- sort(F)
    n <- length(x)
    k <- seq(1, n, by = 1)
    Pteo <- (k - 0.5)/n
    plot(Pemp, Pteo, xlab = xlab, ylab = ylab, col = 4, 
         xlim = c(0, 1), ylim = c(0, 1), main = main, lwd = 2, ...)
    if(line)
        lines(c(0, 1), c(0, 1), col = 2, lwd = 2)
    Cor.Coeff <- cor(Pemp, Pteo)
    Determination.Coeff <- (Cor.Coeff^2) * 100
    return(list(Cor.Coeff = Cor.Coeff, Determination.Coeff = Determination.Coeff))
}
## ****************************************************************************
## Akaike information criterium (AIC) and
## Bayesian information criterion (BIC) for flexible Weibull distribution
abic.flex.weibull <- function(x, alpha.est, beta.est)
{ 
    alpha <- alpha.est 
    beta <- beta.est 
    n <- length(x)
    p <- 2
    f <- dflex.weibull(x, alpha, beta)
    l <- log(f)
    LogLik <- sum(l)
    AIC<- - 2 * LogLik  + 2 * p 
    BIC<- - 2 * LogLik + p * log(n)                   
    return(list(LogLik=LogLik,AIC=AIC, BIC=BIC))
} 
## **************************************************************************** 
