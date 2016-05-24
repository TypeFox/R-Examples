## **************************************************************************
## Probability density function(pdf) of Marshall-Olkin Extended Exponential(MOEE) distribution
dmoee <- function (x, alpha, lambda, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- exp(log(lambda) + log(x))
    pdf <- exp(log(alpha) + log(lambda) - u - 2 * log(1.0 - (1.0-alpha) * exp(-u)))
    if (log) 
        pdf <- log(pdf)
    return(pdf)
}     
## **************************************************************************
## Cummulative distribution function(cdf) of Marshall-Olkin Extended Exponential(MOEE) distribution
pmoee <- function(q, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (q <= 0))
        stop("Invalid arguments")
    u <- exp(log(lambda) + log(q))
    cdf <- exp(log(1.0 - exp(-u)) - log(1.0 - (1.0-alpha) * exp(-u)))
    if (!lower.tail) 
        cdf <- 1.0 - cdf
    if (log.p)
        cdf <- log(cdf)
    return(cdf)
}
## **************************************************************************
## Quantile function of Marshall-Olkin Extended Exponential(MOEE) distribution
qmoee <- function(p, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl <- (1.0/lambda) * log(1.0 + (alpha * p / (1.0 - p)))   
    if (!lower.tail) 
        qtl <- (1.0 / lambda) * log(1.0 + (alpha * (1 - p) / p)) 
    if (log.p)
        qtl <- log(qtl)
    return(qtl)   
}
## **************************************************************************
## Random variate generation from Marshall-Olkin Extended Exponential(MOEE) distribution
rmoee <- function(n, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (n <= 0))
        stop("Invalid arguments")
    return((1 / lambda) * log((alpha / runif(n))+ (1.0 - alpha)))
}
## **************************************************************************
## Reliability function of Marshall-Olkin Extended Exponential(MOEE) distribution
smoee <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))
        stop("Invalid arguments")     
    u <- exp(log(lambda) + log(x))
    cdf <- exp(log(1.0 - exp(-u)) - log(1.0 - (1.0-alpha) * exp(-u)))
    return(1.0 - cdf)
}  
## **************************************************************************
## Hazard function of Marshall-Olkin Extended Exponential(MOEE) distribution
hmoee <- function (x, alpha, lambda)
{
    if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0)) 
        stop("Invalid arguments")       
    u <- exp(log(lambda) + log(x))
    num <- exp(log(alpha) + log(lambda) -u - 2 * log(1.0 - (1.0-alpha) * exp(-u)))
    den <- 1.0 - exp(log(1.0 - exp(-u)) - log(1.0 - (1.0-alpha) * exp(-u)))     
    return(num/den) 
}  
## **************************************************************************
## Hazard rate average function of Marshall-Olkin Extended Exponential(MOEE) distribution
hra.moee <- function(x, alpha, lambda)
{
    return(-pmoee(x, alpha, lambda, lower.tail = FALSE, log.p = TRUE) / x)
}
## **************************************************************************
## Conditional Hazard rate function of Marshall-Olkin Extended Exponential(MOEE) distribution
crf.moee <- function(x, t=0, alpha, lambda)
{
    t <- t
    x <- x
    nume <- hmoee(x+t, alpha, lambda)
    deno <- hmoee(x, alpha, lambda)
    return(nume/deno)
}
## **************************************************************************
## Kolmogorov-Smirnov test (One-sample)for Marshall-Olkin Extended Exponential(MOEE) distribution
ks.moee <- function(x, alpha.est, lambda.est, 
                    alternative = c("less","two.sided", "greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    lambda <- lambda.est
    res <- ks.test(x, pmoee, alpha, lambda, alternative = alternative)
    if(plot)
    {
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs',
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pmoee(t, alpha, lambda)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## **************************************************************************
## Quantile-Quantile(QQ) plot for Marshall-Olkin Extended Exponential(MOEE) distribution
qq.moee <- function(x, alpha.est, lambda.est, main = ' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    lambda <- lambda.est       
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k - 0.5) / n
    Finv <- qmoee(P, alpha, lambda)
    quantiles <- sort(x)
    limx <- range(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, ylim = limx, 
         main = main, col = 4, lwd = 2)
    lines(c(0,limx), c(0,limx), col = 2,lwd = 2, ...)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]   # first quartile Q1
        x2 <- quant[4]   # third quartile Q3
        y1 <- qmoee(0.25, alpha, lambda)
        y2 <- qmoee(0.75, alpha, lambda)
        m <- ((y2-y1) / (x2-x1)) 
        inter <- y1 - (m*x1)
        abline(inter, m, col = 4, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))
}   
## **************************************************************************
## Probability-Probability(PP) plot for Marshall-Olkin Extended Exponential(MOEE) distribution
pp.moee <- function(x, alpha.est, lambda.est, main = ' ', line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est
    lambda <- lambda.est
    F <- pmoee(x, alpha, lambda)
    Pemp <- sort(F)
    n <- length(x)
    k <- seq(1, n, by = 1)
    Pteo <- (k - 0.5) / n
    plot(Pemp, Pteo, xlab = xlab, ylab = ylab, col = 4, 
         xlim = c(0,1), ylim = c(0,1), main = main, lwd = 2, ...)
    if(line)
        lines(c(0,1), c(0,1), col = 2, lwd = 2)
    Cor.Coeff <- cor(Pemp,Pteo)
    Determination.Coeff <- (Cor.Coeff^2) * 100
    return(list(Cor.Coeff = Cor.Coeff, Determination.Coeff = Determination.Coeff))
} 
## **************************************************************************
## Akaike information criterium (AIC) and
## Bayesian information criterion (BIC) for Marshall-Olkin Extended Exponential(MOEE) distribution
abic.moee <- function(x, alpha.est, lambda.est) { 
    alpha <- alpha.est
    lambda <- lambda.est
    n <- length(x)
    p <- 2
    f <- dmoee(x,alpha, lambda)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p           
    BIC <- - 2 * LogLik + p * log(n)       
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## **************************************************************************
