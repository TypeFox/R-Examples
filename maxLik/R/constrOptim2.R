#  This file is a modified copy of src/library/stats/R/constrOptim.R
#  Part of the R package, http://www.R-project.org

### This foutine is not intended for end-user use.
### API is subject to change.
constrOptim2<-function(theta,
                       f,grad=NULL,
                       ineqA,ineqB,
                       mu=0.0001,control=list(),
                  method=if(is.null(grad)) "Nelder-Mead" else "BFGS",
                  outer.iterations=100,outer.eps=0.00001,
                       ...){
   ## Optimize with inequality constraint using SUMT/logarithmic
   ## barrier
   ## 
   ## start      initial value of parameters, included the fixed ones
   ##
   ## This function has to operate with free parameter components
   ## only as 'optim' cannot handle
   ## fixed parameters.  However, for computing constraints in
   ## 'R' and 'dR' we have to use the complete parameter vector.
   ## 
   R <- function(thetaFree, thetaFree.old, ...) {
      ## Wrapper for the function.  As this will be feed to the
      ## 'optim', we have to call it with free parameters only
      ## (thetaFree) and internally expand it to the full (theta)
      ## 
      ## Were we called with 'fixed' argument in ... ?
      dotdotdot <- list(...)
                           # can this be made better?
      fixed <- dotdotdot[["fixed"]]
      theta <- addFixedPar( theta = thetaFree, start = theta0, fixed = fixed)
      theta.old <- addFixedPar( theta = thetaFree.old, start = theta0, fixed = fixed)
       ineqA.theta<-ineqA%*%theta
       gi<- ineqA.theta + ineqB
      if(any(gi < 0))
           ## at least one of the constraints not fulfilled
           return(NaN)
       gi.old <- ineqA%*%theta.old + ineqB
      bar <- sum(gi.old*log(gi) - ineqA.theta)
                           # logarithmic barrier value: sum over
                           # components
      if(!is.finite(bar))
          bar<- -Inf
      result <- f(thetaFree, ...)-mu*bar
                           # do not send 'fixed' and 'start' to the
                           # function here -- we have already
                           # expanded theta to the full parameter
      result
    }
    dR<-function(thetaFree, thetaFree.old, ...){
       ## Wrapper for the function.  As this will be feed to the 'optim',
       ## we have to call it with free parameters only (thetaFree) and
       ## internally expand it to the full (theta)
       ## 
       ## Were we called with 'fixed' argument in ... ?
       dotdotdot <- list(...)
                           # can this be made better?
       fixed <- dotdotdot[["fixed"]]
       theta <- addFixedPar( theta = thetaFree, start = theta0, fixed = fixed)
       theta.old <- addFixedPar( theta = thetaFree.old, start = theta0, fixed = fixed)
       ineqA.theta<-ineqA%*%theta
        gi<-drop(ineqA.theta + ineqB)
        gi.old<-drop(ineqA%*%theta.old + ineqB)
        dbar<-colSums( ineqA*gi.old/gi-ineqA)
       if(!is.null(fixed))
           gr <- grad(thetaFree,...)- (mu*dbar)[!fixed]
                           # grad only gives gradient for the free parameters in order to maintain
                           # compatibility with 'optim'.  Hence we compute barrier gradient
                           # for the free parameters only as well.
       else
           gr <- grad(thetaFree,...)- (mu*dbar)
       return(gr)
    }
    if (!is.null(control$fnscale) && control$fnscale<0)
      mu <- -mu ##maximizing
    if(any(ineqA%*%theta + ineqB < 0))
        stop("initial value not the feasible region")
    theta0 <- theta
                           # inital value, for keeping the fixed params
    ## Were we called with 'fixed' argument in ... ?
    fixed <- list(...)[["fixed"]]
    if(!is.null(fixed))
        thetaFree <- theta[!fixed]
    else
        thetaFree <- theta
    ##
    obj<-f(thetaFree, ...)
    r<-R(thetaFree,thetaFree,...)
    for(i in 1L:outer.iterations){
        obj.old<-obj
        r.old<-r
        thetaFree.old<-thetaFree
        fun<-function(thetaFree,...){ R(thetaFree,thetaFree.old,...)}
        
        if( method == "SANN" ) {
          if( is.null( grad ) ) {
            gradient <- NULL
          } else {
            gradient <- grad
          }
        } else {
          gradient <- function(thetaFree, ...) {
            dR(thetaFree, thetaFree.old, ...)
          }
        }
        ## As 'optim' does not directly support fixed parameters, 
        a<-optim(par=thetaFree.old,fn=fun,gr=gradient,control=control,method=method,...)
        r<-a$value
        if (is.finite(r) && is.finite(r.old) && abs(r-r.old)/(outer.eps+abs(r-r.old))<outer.eps)
            break
        thetaFree<-a$par
        obj<-f(thetaFree, ...)
        if (obj>obj.old) break
    }
    if (i==outer.iterations){
        a$convergence<-7
        a$message<-"Barrier algorithm ran out of iterations and did not converge"
    }
    if (mu>0 && obj>obj.old){
        a$convergence<-11
        a$message<-paste("Objective function increased at outer iteration",i)
    }
    if (mu<0 && obj<obj.old){
        a$convergence<-11
        a$message<-paste("Objective function decreased at outer iteration",i)
    }

        
    a$outer.iterations<-i
    a$barrier.value<-a$value
    a$value<-f(a$par, ...)
    a$barrier.value<-a$barrier.value-a$value
    a
    
}
