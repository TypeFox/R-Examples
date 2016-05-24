## Function which optimizes the location of the new mixture component
## __input__
## FUN      : [function] to be optimized (either 'KERNEL' or 'fn.w')
## mu0      : [kx1 vector] starting vector
## control  : [list] of optimization parameters (see optim)
## ...      : additional parameters used by 'KERNEL'
## __output__
## [list] with the following components:
## $mu      : [kx1 vector] of location parameters
## $Sigma   : [k^2x1 vector] of scale matrices (in vector form)
## $value   : [double] value of the function at optimum
## $method  : [character] method used in the optimization ("BFGS", "Nelder-Mead", "IS")
## $time    : [double] time of the optimization
## __20080502__
'fn.optmu' <- function(FUN, mu0, control, ...)
  {
    'fn.optmu_sub' <- function(method)
      tmp <- optim(par=mu0, fn=FUN, method=method,
                   control=list(
                     trace=control$trace,
                     maxit=control$maxit,
                     reltol=control$reltol,
                     fnscale=-1), ## maximization
                   hessian=TRUE, log=TRUE, ...)

    ptm <- proc.time()[3]
    k <- length(mu0)
    method <- "BFGS"
    Sigma <- NA

    ## first pass optimization
    tmp <- fn.optmu_sub(method)
    
    if (tmp$convergence>1)
      { ## if bad convergence
        if (k>1)
          { ## and multivariate optimization, then use Nelder-Mead optimization
            method <- "Nelder-Mead"
            tmp <- fn.optmu_sub(method)
          }
      }

    if (tmp$convergence>1)
      { ## if still bad convervence        
        method <- "IS"
      }
    
    if (fn.isSingular(tmp$hessian))
      { ## or if the Hessian matrix is singular
        method <- "IS"
      }
    else
      { 
        Sigma <- -solve(tmp$hessian)
        if (!fn.isPD(Sigma))
          { ## or if the scale matrix is not positive definite
            method <- "IS"
          }
      }
    
    list(mu=as.vector(tmp$par),
         Sigma=matrix(Sigma, nrow=1),
         value=as.numeric(tmp$value),
         method=as.character(method),
         time=as.numeric(proc.time()[3]-ptm))
  }

