# Mandel-Paule algorithm
# Based on code by S Cowen
#
# 2011-03-26 amended tolerance handling
# to prevent 'small' variances being ignored
#
# 2011-04-15 amended to include df based on count(x)
# 
mandel.paule <- function(x, ..., tol=.Machine$double.eps^0.25, maxiter=25) {
        mpaule(x, ..., tol=tol, maxiter=maxiter)
}

mpaule <- function(x, ..., tol=.Machine$double.eps^0.25, maxiter=25) {
        UseMethod("mpaule")
}

mpaule.default <- function(x, u=NULL, n=NULL, groups=NULL, 
        tol=.Machine$double.eps^0.25, maxiter=25, ...)          
{
#If n present, u is interpreted as sd of n observations
        
        if(tol >= 1.0 ) stop("Tolerance must be <1.0")
                
        count<-function(x) sum(!is.na(x))
        
        if(!is.null(groups)) {
                groups <- factor(groups)
                x.i <- as.vector(tapply(x, groups, mean, na.rm=TRUE))
                u.i <- as.vector(tapply(x, groups, sd, na.rm=TRUE))
                u.i <- u.i/sqrt(as.vector(tapply(x, groups, count)))
        } else {
                x.i <- x
                if(is.null(n)) 
                        u.i <- u
                else
                        u.i <- u/sqrt(rep(n, length(x.i)))  #Recycles n

        }

        cons.mean <- NA
                #Guarantees cons.mean exists in case no iterations are run
                
        v <- v.x <- var(x.i)
        iter <- 0
        dv <- v
        converged <- 0L 
        while(iter < maxiter && abs(dv) > tol*v.x )     
        {
                iter<- iter + 1
                wt <- 1 / ( u.i^2 + v )
                cons.mean <- sum(wt * x.i) / sum(wt)
                F <- sum( wt * (x.i - cons.mean)^2 ) - (length(x.i) - 1)        
                dv <- F / sum( wt^2 * (x.i - cons.mean)^2 )
                v <- v + dv
                if(v < 0) {
                        v <- 0.0
                        dv <- 0.0
                        wt <- 1 / ( u.i^2 + v )
                        cons.mean <- sum(wt * x.i) / sum(wt)
                        converged <- 2L
                }
                
        }
        
        if(abs(dv) >= tol*v.x) {
                warning("Maximum iterations reached; M-P may not have converged", call.=TRUE)
                
        } else {
                if(converged == 0L) converged<-1L
                        #Not changed if already set to 2L
        }
        

        rv <- .construct.loc.est( x=cons.mean, u=1/sqrt(sum(wt)), df=count(x.i)-1, 
                xi=x.i, ui=u.i, u.eff=sqrt(v+u.i^2), w=rep(1, length(x.i)), method="Mandel-Paule", 
                method.details=list(var.between=v, iter=iter, converged=converged, tol=tol*v.x))
        
        return(rv)
}
