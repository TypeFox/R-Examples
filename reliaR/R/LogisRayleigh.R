## **************************************************************************
## Probability density function(pdf) of Logistic-Rayleigh distribution
dlogis.rayleigh <- function (x, alpha, lambda, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp((lambda * x * x)/2.0)
    num <- lambda * alpha * u * ((u - 1.0)^(alpha - 1.0)) * x
    deno <- (1.0 + (u - 1.0) ^ alpha) ^ 2.0
    pdf <- num/deno
    if(log)
        pdf <- log(pdf)
    return(pdf)
}
## **************************************************************************
## Cummulative distribution function(cdf) of Logistic-Rayleigh distribution
plogis.rayleigh <- function (q, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (q <= 0))    
        stop("Invalid arguments")
    u <- exp((lambda * q * q)/2.0)
    tmp <- 1.0 + ((u - 1.0) ^ alpha)
    cdf <- 1.0 - (1.0/tmp)
    if(!lower.tail)
        cdf <- 1.0 - cdf
    if(log.p)
        cdf <- log(cdf)
    return(cdf)
}
## **************************************************************************
## Quantile function of Logistic-Rayleigh distribution
qlogis.rayleigh <- function (p, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl <- ((2.0/lambda)*(log(1.0 +((p/(1.0-p)) ^ (1.0/alpha))))) ^ 0.5   
    if (!lower.tail) 
        qtl<- ((2.0/lambda)*log(1.0 +(((1.0-p)/p) ^ (1.0/alpha)))) ^ 0.5
    if (log.p) 
        qtl<- log(qtl)    
    return(qtl)   
}
## **************************************************************************
## Random variate generation from Logistic-Rayleigh distribution
rlogis.rayleigh<-function(n, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (n <= 0))
        stop("Invalid arguments")
    u <- runif(n)
    return(((2.0/lambda)*(log(1.0 +((u/(1.0 - u)) ^ (1.0/alpha))))) ^ 0.5)  
}
## **************************************************************************
## Reliability function of Logistic-Rayleigh distribution
slogis.rayleigh <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp((lambda * x * x)/2.0)
    tmp <- 1.0 +((u - 1.0) ^ alpha)
    relia <- 1.0/tmp      
    return(relia)   
}
## **************************************************************************
## Hazard function of Logistic-Rayleigh distribution
hlogis.rayleigh <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp((lambda * x * x)/2.0)
    nume <- lambda * alpha * u * ((u - 1.0) ^ (alpha - 1.0)) * x
    deno <- 1.0 + (u - 1.0) ^ alpha 
    return(nume/deno)   
} 
## **************************************************************************
## Hazard rate average function of Logistic-Rayleigh distribution
hra.logis.rayleigh <-function(x, alpha, lambda)
{
  r <- slogis.rayleigh(x, alpha, lambda)
  fra <-((-1) * log(r))/x
  return(fra)
}
## **************************************************************************
## Conditional Hazard rate function of Logistic-Rayleigh distribution
crf.logis.rayleigh <-function(x, t=0, alpha, lambda)
{
    t <- t
    x <- x
    nume <- hlogis.rayleigh(x + t, alpha, lambda)
    deno <- hlogis.rayleigh(x, alpha, lambda)
    return(nume/deno)
}
## **************************************************************************
## Kolmogorov-Smirnov test (One-sample)for Logistic-Rayleigh distribution
ks.logis.rayleigh <- function(x, alpha.est, lambda.est, 
                   alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    lambda <- lambda.est
    res <- ks.test(x, plogis.rayleigh, alpha, lambda, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- plogis.rayleigh(t, alpha, lambda)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## **************************************************************************
## Quantile-Quantile(QQ) plot for Logistic-Rayleigh distribution
qq.logis.rayleigh <- function(x, alpha.est, lambda.est, main = ' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k-0.5)/n   
    limx <- c(min(x),  max(x))
    Finv <- qlogis.rayleigh(P, alpha, lambda)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qlogis.rayleigh(0.25, alpha, lambda)
        y2 <- qlogis.rayleigh(0.75, alpha, lambda)
        m <- ((y2 - y1) / (x2 - x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))
}
## **************************************************************************
## Probability-Probability(PP) plot for Logistic-Rayleigh distribution
pp.logis.rayleigh <- function(x, alpha.est,lambda.est,main = ' ',line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est
    lambda <- lambda.est
    F <- plogis.rayleigh(x, alpha, lambda)
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
##**********************************************************************
## Akaike information criterium (AIC) and
## Schwartz/Bayesian information criterion (BIC) 
## for Logistic-Rayleigh distribution
abic.logis.rayleigh <- function(x, alpha.est, lambda.est)
{     
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    p <- 2
    f <- dlogis.rayleigh(x, alpha, lambda)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p  
    BIC <- - 2 * LogLik + p * log(n)   
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## **************************************************************************
