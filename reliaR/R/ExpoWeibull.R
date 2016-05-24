## *****************************************************************************
## Probability density function(pdf) of Exponentiated Weibull(EW) distribution
dexpo.weibull <- function (x, alpha, theta, log=FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (x <= 0))    
        stop("Invalid arguments")
     u <- exp(alpha * log(x))
    pdf <- exp(log(alpha)+log(theta)-log(x)+log(u)-u +(theta -1.0)* log(1.0 - exp(-u)))    
    if (log) 
       pdf<- log(pdf)
    return(pdf)   
}
## *****************************************************************************
## Cummulative distribution function(cdf) of Exponentiated Weibull distribution
pexpo.weibull <- function (q, alpha, theta, lower.tail=TRUE, log.p=FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (q <= 0))    
        stop("Invalid arguments")
     u <- exp(alpha * log(q))
    cdf<-exp(theta * log(1.0 - exp(-u)))    
    if (!lower.tail) 
        cdf<- 1.0 - cdf 
    if (log.p)  
        cdf<- log(cdf)    
    return(cdf)   
}
## *****************************************************************************
## Quantile function of Exponentiated Weibull distribution
qexpo.weibull <- function (p, alpha, theta, lower.tail=TRUE, log.p=FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl<- ((- log(1.0 - (p ^ (1.0/theta)))) ^ (1.0/alpha))    
    if (!lower.tail) 
    { 
      qtl<- ((- log(1.0 - ((1.0-p) ^ (1.0/theta)))) ^ (1.0/alpha)) 
    }    
    if (log.p) 
      qtl<- log(qtl)    
    return(qtl)   
}
## *****************************************************************************
## Random variate generation from Exponentiated Weibull distribution
rexpo.weibull<-function(n, alpha, theta)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (n <= 0))
        stop("Invalid arguments")
    return(((- log(1.0 - (runif(n) ^ (1.0/theta)))) ^ (1.0/alpha))) 
}
## ***************************************************************************** 
## Reliability function of Exponentiated Weibull distribution
sexpo.weibull <- function (x, alpha, theta)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (x <= 0))
        stop("Invalid arguments")    
       u <- exp(alpha * log(x))               
      return(1.0 - exp(theta * log(1.0 - exp(-u))))   
}
## *****************************************************************************
## Hazard function of Exponentiated Weibull distribution
hexpo.weibull <- function (x, alpha, theta)
{
    if((!is.numeric(alpha)) || (!is.numeric(theta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(theta) <= 0) || (x <= 0))
        stop("Invalid arguments")         
      u <- exp(alpha * log(x))  
      num<-exp(log(alpha)+log(theta)-log(x)+log(u)-u +
                     (theta -1.0)* log(1.0 - exp(-u)))      
      den <- 1.0 - exp(theta * log(1.0 - exp(-u)))         
      return(num/den)   
} 
## *****************************************************************************
## Hazard rate average function of Exponentiated Weibull distribution
hra.expo.weibull <- function(x, alpha, theta)
  {
    r <- sexpo.weibull(x, alpha, theta)
    fra <-((-1) * log(r)) / x
    return(fra)
  }
## *******************************************************************
## Conditional Hazard rate function of Exponentiated Weibull distribution
crf.expo.weibull <- function(x, t=0, alpha, theta)
{
    t <- t
    x <- x
    nume <- hexpo.weibull(x+t, alpha, theta)
    deno <- hexpo.weibull(x, alpha, theta)
    return(nume/deno)
  }
## *****************************************************************************
## Kolmogorov-Smirnov test (One-sample)for Exponentiated Weibull distribution
ks.expo.weibull <- function(x, alpha.est, theta.est, 
                  alternative = c("less","two.sided","greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    theta <- theta.est
    res <- ks.test(x, pexpo.weibull, alpha, theta, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pexpo.weibull(t, alpha, theta)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## *****************************************************************************
## Quantile-Quantile(QQ) plot for Exponentiated Weibull distribution
qq.expo.weibull <- function(x, alpha.est, theta.est, main=' ', line.qt=FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    theta <- theta.est       
    n <- length(x)
    k <- seq(1,n,by=1)
    P <- (k - 0.5) / n 
    limx <- c(min(x), max(x))   
    Finv <- qexpo.weibull(P,alpha, theta)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qexpo.weibull(0.25, alpha, theta)
        y2 <- qexpo.weibull(0.75, alpha, theta)
        m <-((y2-y1)/(x2-x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))     
  }
## *****************************************************************************
## Probability-Probability(PP) plot for Exponentiated Weibull distribution
pp.expo.weibull <- function(x, alpha.est, theta.est, main=' ', line=FALSE, ...)
{
      xlab <- 'Empirical distribution function'
      ylab <- 'Theoretical distribution function'
      alpha <- alpha.est
      theta <- theta.est
      F <- pexpo.weibull(x,alpha, theta)
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
## *********************************************************************
## Akaike information criterium (AIC) and
## Bayesian information criterion(BIC) for Exponentiated Weibull distribution
abic.expo.weibull <- function(x, alpha.est, theta.est)
{ 
    alpha <- alpha.est
    theta <- theta.est
    n <- length(x)
    p <- 2
    f <- dexpo.weibull(x,alpha, theta)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p          
    BIC <- - 2 * LogLik + p * log(n)      
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## *************************************************************************
