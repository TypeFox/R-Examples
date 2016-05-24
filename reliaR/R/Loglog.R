## ***************************************************************************
## Probability density function(pdf) of Loglog(Pham) Distribution  
dloglog <- function (x, alpha, lambda, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- (x ^alpha)
    pdf<-exp(log(alpha) + log(log(lambda))+(alpha-1) *log(x)+log(lambda)* u + 1- lambda^u)   
    if(log)
          pdf <- log(pdf)
    return(pdf)   
}
## ***************************************************************************
## Cummulative distribution function(cdf) of Loglog Distribution  
ploglog <- function (q, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (q <= 0))    
        stop("Invalid arguments")
     u <- (q ^alpha)
     cdf<- 1.0 - exp(1.0 - (lambda^u))   
    if(!lower.tail)
        cdf <- 1.0 - cdf
    if(log.p)
        cdf <- log(cdf)
    return(cdf)
 }
## ***************************************************************************
## Quantile function of Loglog Distribution distribution
qloglog <- function (p, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
       tmp<- log(1.0-p)
       qtl<- ((log(1.0 - tmp))/log(lambda))^ (1.0/alpha)   
    if (!lower.tail) 
        { 
          tmp<- log(p)
          qtl<- ((log(1.0 - tmp))/log(lambda))^ (1.0/alpha)
        }    
    if (log.p) 
          qtl<- log(qtl)    
    return(qtl)   
}
## ***************************************************************************
## Random variate generation from Loglog Distribution  
rloglog <- function(n, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (n <= 0))
        stop("Invalid arguments")
    return(-(1.0/lambda) * log(1.0 - (runif(n) ^ (1/alpha))))
}
## *************************************************************************** 
## Reliability function of Loglog Distribution 
sloglog <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")   
     u <- (x ^alpha)        
     return(exp(1.0 - (lambda^u)))   
}
## ***************************************************************************
## Hazard function of Loglog Distribution  
hloglog <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")     
    u <- (x ^alpha)
    hazard <- exp(log(alpha)+log(log(lambda))+(alpha-1.0)*log(x)+ u *log(lambda))          
    return(hazard)   
} 
## ***************************************************************************
## Hazard rate average function of Loglog Distribution 
hra.loglog <- function(x, alpha, lambda)
{
  r <- sloglog(x, alpha, lambda)
  return(((-1)*log(r))/x)
}
## ***************************************************************************
## Conditional Hazard rate function of Loglog Distribution 
crf.loglog <- function(x, t=0, alpha, lambda)
{
    t <- t
    x <- x
    nume <- hloglog(x+t, alpha, lambda)
    deno <- hloglog(x, alpha, lambda)
    return(nume/deno)
}
## ***************************************************************************
## Kolmogorov-Smirnov test (One-sample) for Loglog Distribution
ks.loglog <-function(x, alpha.est, lambda.est,
        alternative=c("less", "two.sided","greater"), plot = FALSE, ...)
{
    alpha<-alpha.est
    lambda<-lambda.est
    res<-ks.test(x, ploglog, alpha, lambda, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- ploglog(t, alpha, lambda)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## ***************************************************************************
## Quantile-Quantile(QQ) plot for Loglog Distribution 
qq.loglog <- function(x, alpha.est, lambda.est, main=' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k - 0.5) / n    
    limx <- c(min(x), max(x))
    Finv <- qloglog(P, alpha, lambda)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qloglog(0.25, alpha, lambda)
        y2 <- qloglog(0.75, alpha, lambda)
        m <- ((y2-y1) / (x2-x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))
}
## ***************************************************************************
## Probability-Probability(PP) plot for Loglog Distribution  
pp.loglog <- function(x, alpha.est, lambda.est, main = ' ', line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est
    lambda <- lambda.est
    F <- ploglog(x, alpha,lambda)
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
## ***************************************************************************
# Akaike information criterium (AIC) and 
# Bayesian information criterion (BIC) for Loglog Distribution
abic.loglog <- function(x, alpha.est, lambda.est)
{ 
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    p <- 2
    f <- dloglog(x, alpha, lambda)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p  
    BIC <- - 2 * LogLik + p * log(n)   
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## ***************************************************************************
