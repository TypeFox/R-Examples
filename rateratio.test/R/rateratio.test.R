"rateratio.test"<-function (x, n, RR = 1, alternative = c("two.sided", "less", 
    "greater"), conf.level = 0.95) 
{
    ## modify checks from prop.test
    DNAME <- deparse(substitute(x))
    if (is.matrix(x)) {
            stop("'x' must be a vector")
    }
    else {
        DNAME <- paste(DNAME, "with time of", deparse(substitute(n)))
        if ((l <- length(x)) != length(n)) 
            stop("'x' and 'n' must have the same length")
    }
    OK <- complete.cases(x, n)
    x <- x[OK]
    n <- n[OK]
    if ((k <- length(x)) != 2) 
        stop("x must have a length 2")
    if (any(n <= 0)) 
        stop("elements of 'n' must be positive")
    if (any(x < 0)) 
        stop("elements of 'x' must be nonnegative")
    if (!is.null(RR)) {
        DNAME <- paste(DNAME, ", null rate ratio ", 
            deparse(substitute(RR)), sep = "")
        #if (length(RR) != l) 
        #    stop("'RR' must have length 1")
        if (RR <= 0) 
            stop("RR must be greater than 0")
    }
    alternative <- match.arg(alternative)
    if ((length(conf.level) != 1) || is.na(conf.level) || (conf.level <= 
        0) || (conf.level >= 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    Y<-x[1]
    N<-n[1]
    X<-x[2]
    M<-n[2]
    ESTIMATE <- c((Y/N)/(X/M),Y/N,X/M)
    names(ESTIMATE) <- c("Rate Ratio","Rate 1","Rate 2")
    pRR<- (N*RR)/(N*RR + M)
    pval.less<- pbinom(Y, X+Y, pRR)
    pval.greater<- pbinom(Y -1, Y+X, pRR, lower.tail = FALSE)

    PVAL <- switch(alternative, 
        less = pval.less, 
        greater = pval.greater,
        two.sided = min(1, 2*min(pval.less,pval.greater) ) )
    
    # Comment out less efficient code
    #p.conf.int<-binom.test(Y,X+Y,conf.level=conf.level,alternative=alternative)$conf.int
    #CINT <- switch(alternative, less = c(0, (p.conf.int[2]*M)/(N*(1-p.conf.int[2]))), 
    #    greater = c( (p.conf.int[1]*M)/(N*(1-p.conf.int[1])), Inf), two.sided = c(
    #    (p.conf.int[1]*M)/(N*(1-p.conf.int[1])), 
    #    (p.conf.int[2]*M)/(N*(1-p.conf.int[2])) ) ) 
    # Modify p.L and p.U from binom.test function
    p.L <- function(x,n,alpha) {
        if (x == 0) 
            0
        else qbeta(alpha, x, n - x + 1)
    }
    p.U <- function(x,n,alpha) {
        if (x == n) 
            1
        else qbeta(1 - alpha, x + 1, n - x)
    }
    CINT <- switch(alternative, less = c(0, (p.U(Y,X+Y,1-conf.level)*M)/(N*(1-p.U(Y,X+Y,1-conf.level) ))), 
        greater = c( (p.L(Y,X+Y,1-conf.level)*M)/(N*(1-p.L(Y,X+Y,1-conf.level))), Inf), two.sided = c(
        (p.L(Y,X+Y,(1-conf.level)/2)*M)/(N*(1-p.L(Y,X+Y,(1-conf.level)/2))), 
        (p.U(Y,X+Y,(1-conf.level)/2)*M)/(N*(1-p.U(Y,X+Y,(1-conf.level)/2))) ) ) 

    attr(CINT, "conf.level") <- conf.level
    names(RR)<-"rate ratio"

    RVAL <- list(
        p.value = as.numeric(PVAL), estimate = ESTIMATE, null.value = RR, 
        conf.int = CINT, alternative = alternative, method ="Exact Rate Ratio Test, assuming Poisson counts", 
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}


#rateratio.test(c(2,9),c(17877,16650))
#rateratio.test(c(2,9),c(17877,16650),alternative="greater",conf.level=.975)

#prompt(rateratio.test,filename="H:\\main\\lpd\\klion\\noninferiority\\r\\rrt\\rateratio.Rd")


