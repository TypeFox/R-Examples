#' Basic auto-correlation estimation of a given sequence
#' @param X sequence 
#' @return first auto-correlation coefficient
autocorr1 <- function(X) {
    #N = length(X)
    #return(cor(X[1:(N-1)],X[2:N]))
    if (!is.mcmc(X)) X=as.mcmc(X)
    autocorr(X,lags=1,relative=FALSE)[1]
}

#' Perform a stationary test to check for an initial burn-in in a sequence
#' @param x sequence
#' @param bridge bridge builder function
#' @param stat statistic of the bridge to use in the test
#' @param param sequence parameters: length 'N' and first auto-correlation coefficient 'rho', or "asymptotic" for default asymptotic parameters, or NULL for auto estimated parameters
#' @param plot boolean to ask for test plots
#' @return A list with class "htest" containing the following components:
#' statistic: the value of the test statistic.
#' p.value: the p-value of the test.
#' method: a character string indicating what type of test was performed.
#' data.name: character string giving the name(s) of the data.
transient.test <- function(x, bridge=studentbridge, stat=E.studentbridge, param=NULL, plot=FALSE ) {
    null.cdf = NULL
    if (is.null(param)) {
        rho = autocorr1(as.mcmc(x))
        N = length(x)
        parameters = paste("length =",N,", estimated auto-correlation =",rho)
        null.cdf = null.param.cdf(deparse(substitute(stat)), N,rho)
    } else {
        if (is.character(param) && param == "asymptotic") {
            null.cdf = null.lim.cdf(deparse(substitute(stat)))
            parameters = "asymptotic distribution"
        } else {
            null.cdf = null.param.cdf(deparse(substitute(stat)), param$N,param$rho)
            parameters = paste("length =",param$N,", auto-correlation =",param$rho)
        }
    }
        
    DNAME <- deparse(substitute(x))
    METHOD <- paste(deparse(substitute(stat)), "statistic for burn-in test on",deparse(substitute(bridge)), "with",parameters)
    b <- bridge(x)
    STATISTIC <- stat(b)
    names(STATISTIC) <- deparse(substitute(stat))
    PVAL <- 1-null.cdf(STATISTIC)
    
    if (isTRUE(plot)) {
        par(mfrow=c(2,2))
        col=rgb(.7,.7,.7)
        
        plot(x,type='l',col=col,xlab="i",ylab="X")
        
        d = density(x)
        plot(x=d$y,y=d$x,type='l',xlab="Density",ylab="X",col=col)
        polygon(x=c(d$y,d$y[1]),y=c(d$x,d$x[1]),col=col,border=NA)
        
        m = max(abs(b))
        plot(b,type='l',col=col,xlab="n",ylab=deparse(substitute(bridge)),ylim=c(-m,m))
        abline(h=STATISTIC,col='black')
        
        b.grid = seq(from=-m,to=m,l=100)
        null.cdf.grid = 1-null.cdf(b.grid)
        plot(x=null.cdf.grid,y=b.grid,type='l',xlab="p-value",ylab=deparse(substitute(stat)),ylim=c(-m,m),col=col)
        abline(h=STATISTIC,col='black')
        abline(v=PVAL,col='black')
    }
    
    return(structure(list(statistic = STATISTIC,
        p.value = PVAL, method = METHOD, data.name = DNAME), class = "htest"))
}
