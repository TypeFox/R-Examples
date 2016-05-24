`gamStat` <-
function(x, return.list=TRUE)
{
    #includes final branching time.  See documentation for info(type '?gam.stat') on how
    # to omit this final branching time.
    if (!is.numeric(x)) stop("object x not of class 'numeric'")
    res <- list()
    x <- rev(sort(x))
    N <- length(x)+1
    b <- sort(x)
    z <- rev(c(b[1], diff(b)))
    T <- sum((2:N)*z)
    s1 <- (1/(N-2) * sum(cumsum((2:(N-1))*z[1:(N-2)])) - T/2)
    res$gamstat = s1/(T*sqrt(1/(12*(N-2))))
    res$pval = pnorm(res$gamstat, mean=0, sd=1)
    res$test = "one-tailed; Ho: rates have not decreased over time"
    if (return.list){
    	cat("------------------------------\n")
    	cat("Calculated gamma:", res$gamstat)
    	cat("\npvalue:", res$pval, "\n")
    	cat("test:", res$test, "\n")
    	cat("*assumes complete taxon sampling.\n");
    	return(res);
    }else{
    	return(res$gamstat);	
    }
}

