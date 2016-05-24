## ***************************************************************************
## Probability density function(pdf) of LFR Distribution  
dlfr <- function (x, alpha, beta, log = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- - ((alpha * x) + (beta * x * x) / 2)
    pdf <- (alpha + (beta * x)) * exp(u)
    if (log) 
       pdf<- log(pdf)
    return(pdf)   
}
## ***************************************************************************
## Cummulative distribution function(cdf) of LFR Distribution  
plfr <- function(q, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(q)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (q <= 0))    
        stop("Invalid arguments") 
    u <- - ((alpha * q) + (beta * q * q) / 2)
    cdf <- 1.0 - exp(u)
    if(!lower.tail)
        cdf <- 1.0 - cdf
    if(log.p)
        cdf <- log(cdf)
    return(cdf)
}
## ***************************************************************************
## Quantile function of Linear Failure Rate (LFR) Distribution 
qlfr <- function(p, alpha, beta, lower.tail = TRUE, log.p = FALSE)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(p)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (p <= 0) || (p > 1))
        stop("Invalid arguments")
    qtl <- (1.0/beta)*(-alpha+((alpha^2.0)-(2.0*beta)*log(1.0-p))^0.5)  
    if (!lower.tail) 
        qtl<-(1.0/beta)*(-alpha+((alpha^2.0)-(2.0*beta)*log(1.0-(1.0-p)))^0.5)     
    if (log.p) 
        qtl<- log(qtl)    
    return(qtl)   
}
## ***************************************************************************
## Random variate generation from LFR Distribution  
rlfr <- function(n, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(n)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (n <= 0))    
        stop("Invalid arguments")
    return((1.0/beta) * (- alpha + ((alpha ^ 2.0) - (2.0 * beta)
               * log(1.0 - runif(n))) ^ 0.5))      
}
## ***************************************************************************
## Reliability function of Linear Failure Rate (LFR) Distribution 
slfr <- function (x, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))    
        stop("Invalid arguments")
    u <- - ((alpha * x) + (beta * x * x) / 2)
    relia <-  exp(u)
    return(relia)  
}
## ***************************************************************************
## Hazard function of Linear Failure Rate (LFR) Distribution 
hlfr <- function (x, alpha, beta)
{
    if((!is.numeric(alpha)) || (!is.numeric(beta)) || (!is.numeric(x)))
        stop("non-numeric argument to mathematical function")
    if((min(alpha) <= 0) || (min(beta) <= 0) || (x <= 0))    
        stop("Invalid arguments")   
    return(alpha + (beta * x))  
}
## ***************************************************************************
## Hazard rate average function of LFR Distribution 
hra.lfr <- function(x, alpha, beta)
{
     r <- slfr(x, alpha, beta)
     fra <- ((-1) * log(r)) / x
     return(fra)
}
## ***************************************************************************
## Conditional Hazard rate function of Linear Failure Rate (LFR) Distribution 
crf.lfr <- function(x, t=0, alpha, beta)
{
    t <- t
    x <- x
    nume <- hlfr(x+t, alpha, beta)
    deno <- hlfr(x, alpha, beta)
    return(nume/deno)
}

## ***************************************************************************
## Kolmogorov-Smirnov test (One-sample) for LFR Distribution 
ks.lfr <- function(x, alpha.est, beta.est, 
           alternative = c("less", "two.sided", "greater"), plot = FALSE, ...)
{
    alpha <- alpha.est 
    beta <- beta.est 
    res<-ks.test(x,plfr, alpha, beta, alternative=alternative)
    if(plot){
        plot(ecdf(x), do.points = FALSE, main = 'Empirical and Theoretical cdfs', 
            xlab = 'x', ylab = 'Fn(x)', ...)
        mini <- min(x)
        maxi <- max(x)
        t <- seq(mini, maxi, by = 0.01)
        y <- plfr(t, alpha, beta)
        lines(t, y, lwd = 2, col = 2)
    }
    return(res)
}
## ***************************************************************************
## Quantile-Quantile(QQ) plot for LFR Distribution 
qq.lfr <- function(x, alpha.est, beta.est, main=' ', line.qt = FALSE, ...)
{
    xlab <- 'Empirical quantiles'
    ylab <- 'Theoretical quantiles'
    alpha <- alpha.est 
    beta <- beta.est        
    n <- length(x)
    k <- seq(1, n, by = 1)
    P <- (k - 0.5)/n
    limx <- c(min(x), max(x))
    Finv <- qlfr(P, alpha, beta)
    quantiles <- sort(x)
    plot(quantiles, Finv, xlab = xlab, ylab = ylab, xlim = limx, 
         ylim = limx, main = main, col = 4, lwd = 2, ...)
    lines(c(0,limx), c(0,limx), col = 2, lwd = 2)
    if(line.qt){
        quant <- quantile(x)
        x1 <- quant[2]
        x2 <- quant[4]
        y1 <- qlfr(0.25, alpha, beta)
        y2 <- qlfr(0.75, alpha, beta)
        m <- ((y2-y1) / (x2-x1))
        inter <- y1 - (m * x1)
        abline(inter, m, col = 2, lwd = 2)
    }
    invisible(list(x = quantiles, y = Finv))
}
## ***************************************************************************
## Probability-Probability(PP) plot for LFR Distribution  
pp.lfr <- function(x, alpha.est, beta.est, main=' ', line = FALSE, ...)
{
    xlab<-'Empirical distribution function'
    ylab<-'Theoretical distribution function'
    alpha <- alpha.est 
    beta <- beta.est              
    F <- plfr(x, alpha, beta)
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
## Akaike information criterium (AIC)  and
## Bayesian information criterion (BIC) for LFR Distribution  
abic.lfr <- function(x, alpha.est, beta.est)
{ 
    alpha <- alpha.est 
    beta <- beta.est 
    n <- length(x)
    p <- 2
    f <- dlfr(x, alpha, beta)
    l <- log(f)
    LogLik <- sum(l)
    AIC<- - 2 * LogLik  + 2 * p 
    BIC<- - 2 * LogLik + p * log(n)                                                       
    return(list(LogLik = LogLik, AIC = AIC, BIC = BIC))
}
## **************************************************************************
  
