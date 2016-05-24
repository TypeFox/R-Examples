ss.aipe.sc.ancova <-function(Psi=NULL, sigma.anova=NULL, sigma.ancova=NULL, psi=NULL, ratio=NULL, rho=NULL, divisor="s.ancova", c.weights, width, conf.level=.95, assurance=NULL, ...)
  
  
{options(warn=-1)
if(divisor!="s.ancova" && divisor!="s.anova") stop("The argument 'divisor' must be either 's.ancova' or 's.anova'")
if(!is.null(ratio))
    {if(ratio>1 | ratio <0) stop("'ratio' must be larger than 0 and smaller than or equal to 1.")
    if(!is.null(rho)) stop("You just need to specify either 'ratio' or 'rho', not both.")
    }
if(!is.null(rho))
    {if(rho>1 | rho< -1) stop("'rho' must be no smaller than -1 and no larger than 1.")
    ratio<-sqrt(1-rho^2)
    }

if(divisor=="s.ancova")
    {options(warn=-1)
    
    if(sum(c.weights)!=0) stop("The sum of the contrast weights must be zero")
    if(sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")
    if(is.null(psi))
        {if(is.null(Psi) | is.null(sigma.ancova)) stop("You must specifiy either 'psi', or both 'Psi' and 'sigma.ancova'.")
        psi<- Psi/sigma.ancova
        }
    
    alpha <- 1-conf.level
    J <- length(c.weights)
    
    if(is.null(assurance))
        {# Initial starting value for n using the z distribution.
        n.0 <- sum(c.weights^2)* 4 *(qnorm(1-alpha/2))^2 /(width)^2 
        
        # Second starting value for n using the central t distribiton. 
        n <- sum(c.weights^2)* 4 *(qt(1-alpha/2, n.0*J-J-1))^2 / (width)^2
        
        # measures the discrepency between the inital and second starting values.
        Difference <- abs(n-n.0)
        
        while(Difference > .000001) 
            { n.p <- n
            n <- sum(c.weights^2)* 4 *(qt(1-alpha/2, n*J-J-1))^2 / (width)^2
            Difference <- abs(n - n.p) 
            }
        n <- ceiling(n)
        
        # Initial estimate of noncentral value.
        lambda.0 <- psi/ sqrt( sum(c.weights^2) / n)
        
        lambda.limits.0 <- conf.limits.nct(ncp=lambda.0, df=n*J-J-1, conf.level=1-alpha)
        psi.limit.upper.0 <- lambda.limits.0$Upper.Limit* sqrt(sum(c.weights^2) / n)
        psi.limit.lower.0 <- lambda.limits.0$Lower.Limit* sqrt(sum(c.weights^2) / n) 
            
        # Initial full-width for confidence interval.
        Diff.width.Full <- abs(psi.limit.upper.0 - psi.limit.lower.0) - width
        
        while(Diff.width.Full > 0)
            { n <- n + 1
            lambda <- psi/ sqrt( sum(c.weights^2) / n )
            lambda.limits <- conf.limits.nct(ncp=lambda, df=n*J-J-1, conf.level=1-alpha)
            psi.limit.upper <- lambda.limits$Upper.Limit * sqrt(sum(c.weights^2) / n)
            psi.limit.lower <- lambda.limits$Lower.Limit * sqrt(sum(c.weights^2) / n)  
            Current.width <- abs(psi.limit.upper - psi.limit.lower)
            Diff.width.Full <- Current.width - width
            }
        return(n)
        }
    
    if(!is.null(assurance))
        {psi <- abs(psi)
        
        if((assurance <= 0) | (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call.=FALSE)
        if(assurance <= .50) stop("The 'assurance' should be larger than 0.5 (but less than 1).", call.=FALSE)
        
        n0 <- ss.aipe.sc.ancova(psi=psi, conf.level=1-alpha, width=width, c.weights=c.weights, assurance=NULL, ...)
        
        lambda.2.sided<-conf.limits.nct(ncp=psi/sqrt( sum(c.weights^2) / n0), df=n0*J-J-1, alpha.upper=(1-assurance)/2, alpha.lower=(1-assurance)/2)
        
        Limit.2.Sided <- sqrt(sum(c.weights^2) / n0) * lambda.2.sided$Upper.Limit
        
        lambda.1.sided<- conf.limits.nct(ncp=psi/sqrt( sum(c.weights^2) / n0), df=n0*J-J-1, alpha.upper=1-assurance, alpha.lower=0)
        
        Limit.1.Sided <- sqrt( sum(c.weights^2) / n0) * lambda.1.sided$Upper.Limit
        
        determine.limit <- function(current.psi.limit=current.psi.limit, sample.size=n0, psi=psi, assurance=assurance)
            { Less <- pt(q=-current.psi.limit/sqrt( sum(c.weights^2) / sample.size),
            df=J*sample.size-J-1, ncp=psi/sqrt( sum(c.weights^2) / sample.size) )
            
            Greater <- 1 - pt(q=current.psi.limit/ sqrt( sum(c.weights^2) / sample.size), 
            df=J*sample.size-J-1, ncp=psi/ sqrt( sum(c.weights^2) / sample.size) )
            
            Expected.Widths.Too.Large <- Less + Greater
            return((Expected.Widths.Too.Large - (1-assurance))^2)
            }
        Optimize.Result <- optimize(f=determine.limit, interval=c(Limit.1.Sided, Limit.2.Sided), psi=psi, assurance=assurance)
        
        n <- ss.aipe.sc.ancova(psi=Optimize.Result$minimum, conf.level=1-alpha, width=width, c.weights=c.weights, assurance=NULL, ...)
        
        return(n)
        }
    }
####################################################################################
####################################################################################

if(divisor=="s.anova")
    {options(warn=-1)
    
    if(sum(c.weights)!=0) stop("The sum of the contrast weights must be zero")
    if(sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")
    
    if(!is.null(Psi)| !is.null(sigma.anova) | !is.null(sigma.ancova) )
        {if(is.null(Psi)| is.null(sigma.anova) | is.null(sigma.ancova)) stop("'Psi', 'sigma.anova', and 'sigma.ancova' must be all specified.") }
    
    if(!is.null(psi) | !is.null(ratio))
        {if(is.null(psi)| is.null(ratio)) stop("Both 'psi' and 'ratio' (or 'rho') must be specified")}
    
    alpha <- 1-conf.level
    J <- length(c.weights)
    
    if(is.null(ratio)) 
        {ratio<- sigma.ancova/sigma.anova
        psi<- Psi/sigma.anova    
        }
    
    if(is.null(assurance))
        {# Initial starting value for n using the z distribution.
        n.0 <- 4*sum(c.weights^2)* ratio^2 *(qnorm(1-alpha/2))^2 / (width)^2 
        
        # Second starting value for n using the central t distribiton. 
        n <- 4*sum(c.weights^2) * ratio^2 *(qt(1-alpha/2, n.0*J-J-1))^2 / (width)^2
        
        # measures the discrepency between the inital and second starting values.
        Difference <- abs(n-n.0)
        
        while(Difference > .000001) 
            { n.p <- n
            n <- 4*sum(c.weights^2) * ratio^2 *(qt(1-alpha/2, n*J-J-1))^2 / (width)^2
            Difference <- abs(n - n.p)
            }
        n <- ceiling(n)
        
        # Initial estimate of noncentral value.
        lambda.0 <- psi/ ( ratio* sqrt( sum(c.weights^2)/n ) )
        
        lambda.limits.0 <- conf.limits.nct(ncp=lambda.0, df=n*J-J-1, conf.level=1-alpha)
        psi.limit.upper.0 <- lambda.limits.0$Upper.Limit *ratio*sqrt(sum(c.weights^2)/n)
        psi.limit.lower.0 <- lambda.limits.0$Lower.Limit *ratio*sqrt(sum(c.weights^2)/n) 
        # Initial full-width for confidence interval.
        Diff.width.Full <- abs(psi.limit.upper.0 - psi.limit.lower.0) - width
        
        while(Diff.width.Full > 0)
            { n <- n + 1
            lambda <- psi/( ratio* sqrt( sum(c.weights^2) / n ) )
            lambda.limits <- conf.limits.nct(ncp=lambda, df=n*J-J-1, conf.level=1-alpha)
            psi.limit.upper <- lambda.limits$Upper.Limit* ratio*sqrt(sum(c.weights^2)/n)
            psi.limit.lower <- lambda.limits$Lower.Limit* ratio*sqrt(sum(c.weights^2)/n)  
            Current.width <- abs(psi.limit.upper - psi.limit.lower)
            Diff.width.Full <- Current.width - width
            }
        return(n)
        }
    
    if(!is.null(assurance))
        {psi <- abs(psi)
        
        if((assurance <= 0) | (assurance >= 1)) stop("The 'assurance' must either be NULL or some value greater than zero and less than unity.", call.=FALSE)
        if(assurance <= .50) stop("The 'assurance' should be larger than 0.5 (but less than 1).", call.=FALSE)
        
        n0 <- ss.aipe.sc.ancova(psi=psi, ratio=ratio, divisor="s.anova", conf.level=conf.level, width=width, c.weights=c.weights, assurance=NULL, ...)
        
        lambda.2.sided<-conf.limits.nct(ncp=psi/(ratio* sqrt( sum(c.weights^2)/n0 )), df=n0*J-J-1, alpha.upper=(1-assurance)/2, alpha.lower=(1-assurance)/2 )
        
        Limit.2.Sided <- sqrt( sum(c.weights^2)/n0 )* ratio* lambda.2.sided$Upper.Limit
        
        lambda.1.sided<- conf.limits.nct(ncp=psi/(ratio*sqrt( sum(c.weights^2) / n0)), df=n0*J-J-1, alpha.upper=1-assurance, alpha.lower=0)
        
        Limit.1.Sided <- sqrt( sum(c.weights^2)/n0)* ratio* lambda.1.sided$Upper.Limit
        
        determine.limit <- function(current.psi.limit=current.psi.limit, sample.size=n0, psi=psi, assurance=assurance)
            { Less <-pt(q=-current.psi.limit/(ratio*sqrt(sum(c.weights^2)/sample.size)), df=J*sample.size-J-1, ncp=psi/(ratio*sqrt( sum(c.weights^2) / sample.size)))
            
            Greater<-1-pt(q=current.psi.limit/(ratio*sqrt(sum(c.weights^2)/sample.size)), df=J*sample.size-J-1, ncp=psi/(ratio*sqrt( sum(c.weights^2)/sample.size)))
            
            Expected.Widths.Too.Large <- Less + Greater
            return((Expected.Widths.Too.Large - (1-assurance))^2)
            }
        Optimize.Result <- optimize(f=determine.limit, interval=c(Limit.1.Sided, Limit.2.Sided), psi=psi, assurance=assurance)
        
        n <- ss.aipe.sc.ancova(psi=Optimize.Result$minimum, ratio=ratio, divisor="s.anova", conf.level=1-alpha, width=width, c.weights=c.weights, assurance=NULL, ...)
        
        return(n)
        }
    }

}

