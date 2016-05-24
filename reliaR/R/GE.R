## Probability density function(pdf) of Generalized Exponential distribution
dgen.exp <- function (x, alpha, lambda, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp(log(lambda) + log(x))
    pdf <- exp(log(alpha) + log(lambda) - u + (alpha - 1.0) * log(1.0 - exp(-u)))   
    if(log)
        pdf <- log(pdf)
    return(pdf)
}
##****************************************************************************
## Cummulative distribution function(cdf) of Generalized Exponential distribution
pgen.exp <- function (q, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (q <= 0))    
        stop("Invalid arguments")
    u <- exp(log(lambda) + log(q))
    cdf <- exp(alpha * log(1.0 - exp(-u)))   
    if(!lower.tail)
        cdf <- 1.0 - cdf
    if(log.p)
        cdf <- log(cdf)
    return(cdf)
}
##****************************************************************************
## Quantile function of Generalized Exponential distribution
qgen.exp <- function (p, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl <- - (1.0 / lambda) * log(1.0 - (p ^ (1.0 / alpha)))   
    if(!lower.tail) 
        qtl <- - (1.0/lambda) * log(1.0 - ((1.0 - p) ^ (1.0 / alpha))) 
    if(log.p)
        qtl <- log(qtl)
    return(qtl)
}
##****************************************************************************
## Random variate generation from Generalized Exponential distribution
rgen.exp <- function(n, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (n <= 0))
        stop("Invalid arguments")
    return(-(1.0/lambda) * log(1.0 - (runif(n) ^ (1/alpha))))
}
##****************************************************************************
## Reliability function of Generalized Exponential distribution
sgen.exp <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))
        stop("Invalid arguments")
    u <- exp(log(lambda) + log(x))
    return(1.0 - exp(alpha * log(1.0 - exp(-u))))
}
##****************************************************************************
## Hazard function of Generalized Exponential distribution
hgen.exp <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))
        stop("Invalid arguments")     
    u <- exp(log(lambda) + log(x))
    num <- exp(log(alpha) + log(lambda) - u + (alpha -1.0) * log(1.0 - exp(-u)))   
    den <- 1.0 - exp(alpha * log(1.0 - exp(-u)))     
    return(num/den) 
}
##****************************************************************************
## Hazard rate average function of Generalized Exponential distribution
hra.gen.exp <- function(x, alpha, lambda)
{
    r <- sgen.exp(x, alpha, lambda)
    return(((-1) * log(r)) / x)
}
##****************************************************************************
## Conditional Hazard rate function of Generalized Exponential distribution
crf.gen.exp <- function(x, t = 0, alpha, lambda)
{
    t <- t
    x <- x
    nume <- hgen.exp(x + t, alpha, lambda)
    deno <- hgen.exp(x, alpha, lambda)
    return(nume/deno)
}
##****************************************************************************
## Kolmogorov-Smirnov test (One-sample)for Generalized Exponential distribution
ks.gen.exp <- function(x, alpha.est, lambda.est, 
    alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    lambda <- lambda.est
    res <- ks.test(x, pgen.exp, alpha, lambda, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pgen.exp(t, alpha, lambda)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
##****************************************************************************
## Quantile-Quantile(QQ) plot for Generalized Exponential distribution
qq.gen.exp <- function(x, alpha.est, lambda.est, main = ' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k-0.5)/n   
    limx <- c(min(x),  max(x))
    Finv <- qgen.exp(P, alpha, lambda)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qgen.exp(0.25, alpha, lambda)
        y2 <- qgen.exp(0.75, alpha, lambda)
        m <- ((y2-y1) / (x2-x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))
}
##****************************************************************************
## Probability-Probability (PP) plot for Generalized Exponential distribution
pp.gen.exp <- function(x, alpha.est,lambda.est,main = ' ',line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est
    lambda <- lambda.est
    F <- pgen.exp(x, alpha, lambda)
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
##****************************************************************************
## Akaike information criterium (AIC) and
## Bayesian information criterion (BIC)for Generalized Exponential distribution
abic.gen.exp <- function(x, alpha.est, lambda.est){     
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    p <- 2
    f <- dgen.exp(x, alpha, lambda)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p  
    BIC <- - 2 * LogLik + p * log(n)   
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
###****************************************************************************
