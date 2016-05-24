## ****************************************************************************
## Probability density function(pdf) of Exponential extension distribution
dexp.ext <- function (x, alpha, lambda, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp(log(1.0 + lambda * x))
    pdf <- exp(log(alpha) + log(lambda)+(alpha -1.0)*log(u)+ (1.0 - (u^alpha)))    
    if (log) 
        pdf<- log(pdf)
    return(pdf)   
}
## ****************************************************************************
## Cummulative distribution function(cdf) of Exponential extension distribution
pexp.ext <- function (q, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (q <= 0))    
        stop("Invalid arguments")
     u <- exp(log(1.0 + lambda * q))
    cdf <- 1.0 - exp(1.0-(u^alpha))     
    if (!lower.tail) 
          cdf<- 1.0 - cdf     
    if (log.p) 
          cdf<- log(cdf)    
    return(cdf)   
}
## ****************************************************************************
## Quantile function of Exponential extension distribution
qexp.ext <- function (p, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    tmp <-(1.0 - log(1.0 - p)) ^ (1.0/alpha)
    qtl <- (1.0/lambda) * (tmp - 1.0) 
    if (!lower.tail) 
      { 
        tmp <- (1.0 - log(p)) ^ (1.0/alpha)
        qtl <-  (1.0/lambda) * (tmp - 1.0)  
      }    
    if (log.p) 
        qtl<- log(qtl)    
    return(qtl)   
}
## ****************************************************************************
## Random variate generation from Exponential extension distribution
rexp.ext<-function(n, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (n <= 0))
        stop("Invalid arguments")
    tmp <-(1.0 - log(1.0 - runif(n))) ^ (1.0/alpha)    
    return((1.0/lambda) * (tmp - 1.0))
}
## ****************************************************************************
## Reliability function of Exponential extension distribution
sexp.ext <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))
        stop("Invalid arguments")
    u <- exp(log(1.0 + lambda * x))      
    return(exp(1.0 - (u ^ alpha)))   
}
## ****************************************************************************
## Hazard function of Exponential extension distribution
hexp.ext <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))
        stop("Invalid arguments")          
    u <- exp(log(1.0 + lambda * x))      
    return(exp(log(alpha) + log(lambda) + (alpha - 1.0) * log(u)))   
} 
## ****************************************************************************
## Hazard rate average function of Exponential extension distribution
hra.exp.ext <-function(x, alpha, lambda)
{
      r <- sexp.ext(x, alpha, lambda)
      fra <-((-1) * log(r)) / x
      return(fra)
}
## ******************************************************************
## Conditional Hazard rate function of Exponential extension distribution
crf.exp.ext <-function(x, t=0, alpha, lambda)
{
      t <- t
      x <- x
      nume <- hexp.ext(x+t, alpha, lambda)
      deno <- hexp.ext(x, alpha, lambda)
      return(nume/deno)
}
## ****************************************************************************
## Kolmogorov-Smirnov test (One-sample)for Exponential extension distribution
ks.exp.ext <- function(x, alpha.est, lambda.est, 
               alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    lambda <- lambda.est
    res <- ks.test(x, pexp.ext, alpha, lambda, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pexp.ext(t, alpha, lambda)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## ****************************************************************************
## Quantile-Quantile(QQ) plot for Exponential extension distribution
qq.exp.ext <- function(x, alpha.est, lambda.est, main=' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k-0.5)/n    
    limx <- c(min(x), max(x))
    Finv <- qexp.ext(P, alpha, lambda)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
      y1 <- qexp.ext(0.25, alpha, lambda)
      y2 <- qexp.ext(0.75, alpha, lambda)
      m <- ((y2 - y1)/(x2 - x1))
      inter <- y1 - (m * x1)
      abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))    
  }  
## ****************************************************************************
## Probability-Probability(PP) plot for Exponential extension distribution
pp.exp.ext <- function(x, alpha.est, lambda.est, main=' ', line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est
    lambda <- lambda.est
    F <- pexp.ext(x, alpha, lambda)
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
## ********************************************************************
## Akaike information criterium (AIC)
## Bayesian information criterion (BIC) 
## for Exponential extension distribution
abic.exp.ext <- function(x, alpha.est, lambda.est)
{ 
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    p <- 2
    f <- dexp.ext(x, alpha, lambda)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p           
    BIC <- - 2 * LogLik + p * log(n)    
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## **********************************************************************
