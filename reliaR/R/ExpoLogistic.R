## ***************************************************************************
## Probability density function(pdf) of exponentiated Logistic distribution
dexpo.logistic <- function (x, alpha, beta, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- -(x / beta)
    pdf <- exp(log(alpha) - log(beta) - (alpha + 1.0)* log(1.0 + exp(u)) + u) 
    if (log) 
        pdf <- log(pdf)
    return(pdf)   
}
## ***************************************************************************
## Cummulative distribution function(cdf) of exponentiated Logistic distribution
pexpo.logistic <- function (q, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (q <= 0))    
        stop("Invalid arguments")
    u <- exp(-q / beta)
    cdf <- (1 + u) ^ - alpha     
    if (!lower.tail) 
        cdf <- 1.0 - cdf 
    if (log.p) 
        cdf <- log(cdf)    
    return(cdf)             
}     
## ***************************************************************************
## Quantile function of exponentiated Logistic distribution
qexpo.logistic <- function (p, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl<- - beta * log(p^(-1.0/alpha) - 1.0) 
    if (!lower.tail) 
        qtl<- - beta * log((1.0-p)^(-1.0/alpha) - 1.0)   
    if (log.p) 
        qtl<- log(qtl)    
    return(qtl)   
}
## ***************************************************************************
## Random variate generation from exponentiated Logistic distribution
rexpo.logistic <- function(n, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (n <= 0))
        stop("Invalid arguments")
    return(- beta * log(runif(n)^(-1.0/alpha) - 1.0)) 
}
## **************************************************************************** 
## Reliability function of exponentiated Logistic distribution
sexpo.logistic <- function (x, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))
        stop("Invalid arguments")    
      u <- exp(-x / beta)               
      return(1.0 - ((1 + u) ^ - alpha))   
}
## ****************************************************************************
## Hazard function of exponentiated Logistic distribution
hexpo.logistic <- function (x, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))
        stop("Invalid arguments")       
    u <- exp(-x / beta)
    num <- (alpha / beta) * u * ((1.0 + u)^ -(alpha + 1.0))
    deno <- 1.0 - ((1 + u) ^ - alpha)        
    return(num/deno)   
} 
## ****************************************************************************
## Hazard rate average function of exponentiated Logistic distribution
hra.expo.logistic <- function(x, alpha, beta)
{
    r <- sexpo.logistic(x, alpha, beta)
    fra <- - log(r) / x
    return(fra)
}
## ***************************************************************************
## Conditional Hazard rate function of exponentiated Logistic distribution
crf.expo.logistic <- function(x, t=0, alpha, beta)
{ 
    t <- t
    x <- x
    nume <- hexpo.logistic(x+t, alpha, beta)
    deno <- hexpo.logistic(x, alpha, beta)
    return(nume/deno)
}
## ****************************************************************************
## Kolmogorov-Smirnov test (One-sample)for exponentiated Logistic distribution
ks.expo.logistic <- function(x, alpha.est, beta.est, 
                alternative = c("less","two.sided","greater"), plot = FALSE, ...)
{
    alpha <- alpha.est
    beta <- beta.est
    res<-ks.test(x, pexpo.logistic, alpha, beta, alternative = alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- pexpo.logistic(t, alpha, beta)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## ****************************************************************************
## Quantile-Quantile(QQ) plot for exponentiated Logistic distribution
qq.expo.logistic <- function(x, alpha.est, beta.est, main=' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est
    beta <- beta.est         
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <-(k - 0.5)  /n 
    limx<-c(min(x), max(x))   
    Finv <- qexpo.logistic(P,alpha, beta)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qexpo.logistic(0.25, alpha, beta)
        y2 <- qexpo.logistic(0.75, alpha, beta)
        m <- ((y2-y1)/(x2-x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))     
}
## ****************************************************************************
## Probability-Probability(PP) plot for exponentiated Logistic distribution
pp.expo.logistic <- function(x, alpha.est, beta.est, main=' ', line = FALSE, ...)
{
    xlab <- 'Empirical distribution function'
    ylab <- 'Theoretical distribution function'
    alpha <- alpha.est
    beta <- beta.est
    F <- pexpo.logistic(x, alpha, beta)
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
## Akaike information criterium (AIC) and
## Bayesian information criterion (BIC) 
## for exponentiated Logistic distribution
abic.expo.logistic <- function(x, alpha.est, beta.est)
{ 
    alpha <- alpha.est
    beta <- beta.est
    n <- length(x)
    p <- 2
    f <- dexpo.logistic(x, alpha, beta)
    l <- log(f)
    LogLik <- sum(l)
    AIC <- - 2 * LogLik  + 2 * p          
    BIC <- - 2 * LogLik + p * log(n)     
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
} 
## ************************************************************************
