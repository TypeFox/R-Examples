#<<BEGIN>>
rtrunc <- function(distr=runif, n, linf=-Inf, lsup=Inf,...)
#TITLE Random Truncated Distributions
#DESCRIPTION
# Provides samples from classical \R distributions and \samp{mc2d} specific
# distributions truncated between \samp{linf} (excluded) and \samp{lsup} (included).
#KEYWORDS distribution
#INPUTS
#{distr}<<A function providing random data or its name as character.
#The function 'rdistr' should have a 'qdistr' form (with argument 'p') and a 'pdistr' form
#(with argument 'q'). Example : 'rnorm' (has a 'qnorm' and a 'pnorm' form), 'rbeta', 'rbinom', 'rgamma', ...>>
#{n}<<The size of the sample.>>         .
#[INPUTS]
#{linf}<<A vector of lower bounds.>>
#{lsup}<<A vector of upper bounds, with \samp{lsup < linf} (strictly).>>
#{\dots}<<All arguments to be passed to \samp{pdistr} and \samp{qdistr}.>>
#VALUE
#A vector of \samp{n} values.
#DETAILS
#The function 1) evaluates the \samp{p} values corresponding to \samp{linf} and \samp{lsup} using \samp{pdistr};
#2) samples \samp{n} values using \samp{runif(n, min=pinf, max=psup)}, and 3) takes
#the \samp{n} corresponding quantiles from the specified distribution using \samp{qdistr}.
#
#All distributions (but sample) implemented in the stats library could be used.
#The arguments in \dots should be named. Do not use 'log' or 'log.p' or 'lower.tail'.
#For discrete dictribution, rtrunc sample within \samp{(linf, lsup]}. See example.
#NOTE
#The inversion of the quantile function leads to time consuming functions for some distributions.
#WARNING: The method is flexible, but can lead to problems linked to rounding errors in some extreme situations.
#The function checks that the values are in the expected range and returns an error if not.
#It also warns some extreme situation that could lead to unexpected results.
#See Examples.
#EXAMPLE
#rtrunc("rnorm", n=10, linf=0)
#range(rtrunc(rnorm, n=1000, linf=3, lsup=5, sd=10))
### Discrete distributions
#range(rtrunc(rpois,1000,linf=2,lsup=4,lambda=1))
###Examples of rounding problems. 
###The first one will provide a warning while the results are unexpected, 
###The second will provide an error.
#\dontrun{
#table(rtrunc(rbinom, n=1000, size=10, prob=1-1E-20, lsup=9))
#table(rtrunc(rbinom, n=1000, size=10, prob=1E-14, linf=0))
#}

#CREATED 08-02-20
#REVISED 13-10-01
#--------------------------------------------
{
    linf <- as.vector(linf)
	lsup <- as.vector(lsup)
	if(!is.character(distr)) distr <- as.character(match.call()$distr)          #retrieve the name of the function
    distr <- substr(distr, 2, 1000)                                             #remove the r

    if(any(linf >= lsup)) stop("linf should be < lsup")  #recycle vectors

    pfun <- get(paste("p",distr,sep=""),mode="function")

    pinf <- as.vector(pfun(q=linf,...))
    psup <- as.vector(pfun(q=lsup,...))

    
    p <- runif(n,min=pinf,max=psup)

    qfun <- get(paste("q",distr,sep=""),mode="function")

    res <- as.vector(qfun(p,...))
    # Some possible problem (check if you think to others)
    #
    res[pinf <= 0 & res > lsup] <- NaN          #ex: rtrunc("lnorm",10,linf=-2,lsup=-1)
    res[psup>=1 & res < linf] <- NaN          #ex: rtrunc("unif",10,linf=2,lsup=4,max=1)
    res[is.na(linf) | is.na(lsup)] <- NaN   #ex: rtrunc("norm",10,sd=-2)
    
    
    #Two tests for extreme situations. None Catch all possibilities. THe error is first to avoid the warning
    if(any(res <= linf | res > lsup, na.rm=TRUE)) stop("Error in rtrunc: some values are not in the expected range (maybe due to rounding errors)")
    if(isTRUE(all.equal(pinf,1)) | isTRUE(all.equal(psup,0)) ) warning("Warning: check the results from rtrunc. It may have reached rounding errors")
    
    return(res)
}
