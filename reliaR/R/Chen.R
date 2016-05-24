## **************************************************************************
## Probability density function(pdf) of Chen distribution
dchen <- function (x, beta, lambda, log = FALSE)
{
    if((!is.numeric(beta)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(beta) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp(x ^ beta)
    pdf <- exp(log(beta)+ log(lambda)+ (beta -1.0)*log(x)+ log(u) + lambda *(1.0 - u))
    if(log)
        pdf <- log(pdf)
    return(pdf)
}
## **************************************************************************
## Cummulative distribution function(cdf) of Chen distribution
pchen <- function (q, beta, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(beta)) || (!is.numeric(lambda)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(beta) <= 0) || (min(lambda) <= 0) || (q <= 0))    
        stop("Invalid arguments")
     u <- exp(q ^ beta)
    cdf<- 1.0 - exp(lambda *(1.0 - u))
    if(!lower.tail)
        cdf <- 1.0 - cdf
    if(log.p)
        cdf <- log(cdf)
    return(cdf)
 }
## **************************************************************************
## Quantile function of Chen distribution
qchen <- function (p, beta, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(beta)) || (!is.numeric(lambda)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(beta) <= 0) || (min(lambda) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl <- (log((1.0 -  (log(1.0 - p)/lambda)))) ^ (1.0/beta)
    if(!lower.tail) 
        qtl <-  (log((1.0 -  (log(p)/lambda)))) ^ (1.0/beta)
    if(log.p) 
        qtl <- log(qtl)
    return(qtl)
}
## **************************************************************************
## Random variate generation from Chen distribution
rchen <- function(n, beta, lambda)
{
    if((!is.numeric(beta)) || (!is.numeric(lambda)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(beta) <= 0) || (min(lambda) <= 0) || (n <= 0))    
        stop("Invalid arguments")
    return((log((1.0 -  (log(1.0 - runif(n))/lambda)))) ^ 1.0/beta)
}
## **************************************************************************
## Reliability function of Chen distribution
schen <- function (x, beta, lambda)
{
    if((!is.numeric(beta)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(beta) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp(x ^ beta)
    return(exp(lambda *(1.0 - u)))
}
## **************************************************************************
## Hazard function of Chen distribution
hchen <- function (x, beta, lambda)
{
    if((!is.numeric(beta)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(beta) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp(x ^ beta)
    num <- exp(log(beta)+log(lambda)+(beta -1.0)*log(x)+log(u)+lambda *(1.0 - u))
    den <- exp(lambda *(1.0 - u))
    return(num/den)
}
## **************************************************************************
## Hazard rate average function of Chen distribution
hra.chen <- function(x, beta, lambda)
{
    r <- schen(x, beta, lambda)
    return(((-1)*log(r))/x)
}
## **************************************************************************
## Conditional Hazard rate function of Chen distribution
crf.chen <- function(x, t=0, beta, lambda)
{
      t <- t
      x <- x
      nume <- hchen(x+t, beta, lambda)
      deno <- hchen(x, beta, lambda)
      return(nume/deno)
}
## **************************************************************************
## Kolmogorov-Smirnov test (One-sample) for Chen distribution
ks.chen <- function(x, beta.est, lambda.est,
           alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    beta <- beta.est
    lambda <- lambda.est
    res<-ks.test(x, pchen, beta, lambda, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pchen(t, beta, lambda)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## **************************************************************************
## Quantile-Quantile(QQ) plot for Chen distribution
qq.chen <- function(x, beta.est, lambda.est, main = ' ', line.qt = FALSE, ...)
{
    xlab <-' Empirical quantiles'
    ylab <-' Theoretical quantiles'
    beta <- beta.est
    lambda <- lambda.est
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k-0.5)/n   
    limx <- c(min(x),  max(x))
    Finv <-qchen(P, beta, lambda)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qchen(0.25, beta, lambda)
        y2 <- qchen(0.75, beta, lambda)
        m <- ((y2-y1) / (x2-x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))
}
## **************************************************************************
## Probability-Probability(PP) plot for Chen distribution
pp.chen <- function(x, beta.est, lambda.est, main = ' ', line = TRUE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    beta <- beta.est
    lambda <- lambda.est
    F <- pchen(x, beta, lambda)
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
# Akaike information criterium (AIC) and
# Bayesian information criterion (BIC) for Chen distribution
abic.chen <- function(x, beta.est, lambda.est)
{
    beta <- beta.est
    lambda <-lambda.est
    n <- length(x)
    p <- 2
    f <- dchen(x, beta, lambda)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p  
    BIC <- - 2 * LogLik + p * log(n)   
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## **************************************************************************
