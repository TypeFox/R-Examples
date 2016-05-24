ss.aipe.sc <- function(psi, c.weights, width, conf.level=.95, assurance=NULL, certainty=NULL, ...)
{ options(warn=-1)

if(!is.null(assurance)& !is.null(certainty))
    {if(assurance!=certainty) stop("'assurance' and 'certainty' must have the same value")}
if(!is.null(certainty)) assurance<- certainty 

if(sum(c.weights)!=0) stop("The sum of the coefficients must be zero")
if(sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")

alpha <- 1-conf.level
J <- length(c.weights)

if(is.null(assurance))
    { # Initial starting value for n using the z distribution.
    n.0 <- sum(c.weights^2)*(qnorm(1-alpha/2)/(width/2))^2       

    n <- ceiling(n.0)
    
    # To ensure that the initial n is not too big.
    n <- max(4, n-5)
    
    # Initial estimate of noncentral value. 
    # This is literally the theoretical t-value given psi and the initial estimate of sample size.
    lambda.0 <- psi/sqrt( sum(c.weights^2) / n)
    
    # Initial confidence limits.
    lambda.limits.0 <- conf.limits.nct(ncp=lambda.0, df=n*J-J, conf.level=1-alpha)
    psi.limit.upper.0 <- lambda.limits.0$Upper.Limit*sqrt(sum(c.weights^2) / n)
    psi.limit.lower.0 <- lambda.limits.0$Lower.Limit*sqrt(sum(c.weights^2) / n)  
    
    # Initial full-width for confidence interval.
    Diff.width.Full <- abs(psi.limit.upper.0 - psi.limit.lower.0) - width
    
    while(Diff.width.Full > 0)
    { n <- n + 1
    lambda <- psi/sqrt( sum(c.weights^2) / n )
    lambda.limits <- conf.limits.nct(ncp=lambda, df=n*J-J, conf.level=1-alpha)
    
    psi.limit.upper <- lambda.limits$Upper.Limit*sqrt(sum(c.weights^2) / n)
    psi.limit.lower <- lambda.limits$Lower.Limit*sqrt(sum(c.weights^2) / n) 
    
    Current.width <- abs(psi.limit.upper - psi.limit.lower)
    Diff.width.Full <- Current.width - width
    }
    return(n)
    }

if(!is.null(assurance))
    { psi.sign <- psi
    psi <- abs(psi)
    
    if((assurance <= 0) | (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call.=FALSE)
    if(assurance <= .50) stop("The 'assurance' should be larger than 0.5 (but less than 1).", call.=FALSE)
    
    n0 <- ss.aipe.sc(psi=psi, conf.level=conf.level, width=width, c.weights=c.weights, assurance=NULL, ...)
    
    Limit.2.Sided <- sqrt( sum(c.weights^2) / n0)*conf.limits.nct(ncp=psi/sqrt( sum(c.weights^2) / n0), df=n0*J-J, 
    alpha.upper=(1-assurance)/2, alpha.lower=(1-assurance)/2)$Upper.Limit

    Limit.1.Sided <- sqrt( sum(c.weights^2) / n0) * conf.limits.nct(ncp=psi/sqrt( sum(c.weights^2) / n0), df=n0*J-J, 
    alpha.upper=1-assurance, alpha.lower=0)$Upper.Limit
      
    determine.limit <- function(current.psi.limit=current.psi.limit, sample.size=n0, psi=psi, 
    assurance=assurance)
        { Less <- pt(q=-current.psi.limit/sqrt( sum(c.weights^2) / sample.size),
        df=J*sample.size-J, ncp=psi/sqrt( sum(c.weights^2) / sample.size))
        
        Greater <- 1 - pt(q=current.psi.limit/sqrt( sum(c.weights^2) / sample.size), 
        df=J*sample.size-J, ncp=psi/sqrt( sum(c.weights^2) / sample.size))
        
        Expected.Widths.Too.Large <- Less + Greater
        return((Expected.Widths.Too.Large - (1-assurance))^2)
        }
    Optimize.Result <- optimize(f=determine.limit, interval=c(Limit.1.Sided, Limit.2.Sided), 
    psi=psi, assurance=assurance)
    
    n <- ss.aipe.sc(psi=Optimize.Result$minimum, conf.level=1-alpha, width=width, c.weights=c.weights,
    assurance=NULL, ...)
    return(n)
    }

}
