ss.aipe.sm <- function(sm, width, conf.level=.95, assurance=NULL, certainty=NULL, ...)
{ options(warn=-1)
if(!is.null(assurance)& !is.null(certainty))
    {if(assurance!=certainty) stop("'assurance' and 'certainty' must have the same value")}
if(!is.null(certainty)) assurance<- certainty 

alpha <- 1-conf.level

if(is.null(assurance))
    { # Initial starting value for n using the z distribution.
    n.0 <- (qnorm(1-alpha/2) / (width/2))^2       
    
    n <- ceiling(n.0)
    
    # To ensure that the initial n is not too big.
    n <- max(4, n-5)
    
    # Initial estimate of noncentral value. 
    # This is literally the theoretical t-value given delta and the initial estimate of sample size.
    lambda.0 <- sm*sqrt(n)
    
    # Initial confidence limits.
    lambda.limits.0 <- conf.limits.nct(ncp=lambda.0, df=n-1, conf.level=1-alpha)
    sm.limit.upper.0 <- lambda.limits.0$Upper.Limit/sqrt(n)
    sm.limit.lower.0 <- lambda.limits.0$Lower.Limit/sqrt(n)  
    
    # Initial full-width for confidence interval.
    Diff.width.Full <- abs(sm.limit.upper.0 - sm.limit.lower.0) - width
    
    while(Diff.width.Full > 0)
    { n <- n + 1
    lambda <- sm*sqrt(n)
    lambda.limits <- conf.limits.nct(ncp=lambda, df=n-1, conf.level=1-alpha)
    sm.limit.upper <- lambda.limits$Upper.Limit/sqrt(n)
    sm.limit.lower <- lambda.limits$Lower.Limit/sqrt(n)  
    Current.width <- abs(sm.limit.upper - sm.limit.lower)
    Diff.width.Full <- Current.width - width
    }
    
    #if(warn==FALSE) options(warn=1)
    
    return(n)
    }

if(!is.null(assurance))
    { if((assurance <= 0) | (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call.=FALSE)
    if(assurance <= .50) stop("The 'assurance' should be larger than 0.5 (but less than 1).", call.=FALSE)
    
    n0 <- ss.aipe.sm(sm=sm, conf.level=conf.level, width=width, assurance=NULL, ...)
    
    Limit.2.Sided <- (1/sqrt(n0) ) * conf.limits.nct(ncp=sm*sqrt(n0), df=n0-1, 
    alpha.upper=(1-assurance)/2, alpha.lower=(1-assurance)/2)$Upper.Limit

    Limit.1.Sided <- (1/sqrt(n0) ) * conf.limits.nct(ncp=sm*sqrt(n0), df=n0-1, 
    alpha.upper=1-assurance, alpha.lower=0)$Upper.Limit
      
    determine.limit <- function(current.sm.limit=current.sm.limit, samp.size=n0, sm=sm, 
    assurance=assurance)
        { Less <- pt(q=-current.sm.limit*sqrt(samp.size), df=samp.size-1, ncp=sm/sqrt(samp.size))
        Greater <- 1 - pt(q=current.sm.limit*sqrt(samp.size), df=samp.size-1, ncp=sm/sqrt(samp.size))
        Expected.Widths.Too.Large <- Less + Greater
        return((Expected.Widths.Too.Large - (1-assurance))^2)
        }
    Optimize.Result <- optimize(f=determine.limit, interval=c(Limit.1.Sided, Limit.2.Sided), 
    sm=sm, assurance=assurance)
    
    n <- ss.aipe.sm(sm=Optimize.Result$minimum, conf.level=1-alpha, width=width, assurance=NULL, ...)
    return(n)
    }

}
