## Function which optimizes the probabilities by minimizing
## the squared coefficient of variation
## __input__
## ph       : [(H-1)x1 vector] of past probabilities
## lnK      : [NpxH matrix] of kernel values
## lnD      : [NpxH^2 matrix] of Student-t densities components
## control  : [list] of optimization parameters
## __output__
## [list] with the following components
## $ph      : [Hx1 vector] of optimized probabilities
## $method  : [character] indicating the optimization method used ("NLMIN", "BFGS", "Nelder-Mead", "NONE")
## $time    : [double] indicating the time of the optimization
## __20080427__
'fn.optp' <- function(ph, lnK, lnD, control)
  {    
    ## function which transforms the probabilities (positivity and summability)
    'fn.lambdap' <- function(lambda)
      {
        e <- c(exp(lambda),1)
        as.vector(e/sum(e))
      }

     ## objective function
    'fn.lnf' <- function(lambda)
      {
        r <- .C('fnlnf_C',
                lnp = as.double(log(fn.lambdap(lambda))),
                lnk = as.double(as.vector(t(lnK))),
                lnd = as.double(as.vector(t(lnD))),
                Np = as.integer(Np),
                H = as.integer(H),
                f = as.double(0),
                grad = vector('double',H),
                PACKAGE = 'AdMit',
                NAOK = TRUE,
                DUP = FALSE)

        assign('grad', r$grad, envir = memo)
        as.numeric(r$f)
      }
    
    ## gradient of 'fn.lnf' 
    'fn.gradlnf' <- function(lambda)
      {
        e <- c(exp(lambda),1)
        s <- sum(e)
        tmp <- -e %*% t(e) / s^2
        diag(tmp) <- (e*s-e^2) / s^2
        gradlambda <- as.matrix(tmp[1:H,1:(H-1)])
        grad <- get('grad', envir = memo) ## take gradient from fn.lnf
        
        as.vector(t(gradlambda)%*%grad)
      }
    
    ptm <- proc.time()[3]
    Np <- nrow(lnK)
    H <- ncol(lnK)
    ph <- c((1-control$weightNC)*ph, control$weightNC)
    lambda0 <- log(ph[1:(H-1)]/ph[H])
    memo <- new.env(hash=FALSE)

    method <- "NLMINB"
    tmp <- nlminb(start=lambda0,
                  objective=fn.lnf,
                  gradient=fn.gradlnf,
                  control=list(
                    trace=control$trace,
                    iter.max=control$iter.max,
                    rel.tol=control$rel.tol,
                    x.tol=1e-10))
    
    if (tmp$convergence>0 & length(lambda0)==1)
      { ## if no convergence and univariate optimization
        method <- "BFGS"
      }
    
    if (tmp$convergence>0 & length(lambda0)>1)
      { ## if no convergence and multivariate optimization        
        method <- "Nelder-Mead"
      }
    
    if (tmp$convergence>0)
      { ## run the optimization again with the changed methods
        tmp <- optim(par=lambda0,
                     fn=fn.lnf,
                     gr=fn.gradlnf,
                     method=method,
                     control=list(
                       trace=control$trace,
                       maxit=control$iter.max,
                       reltol=control$rel.tol))
      }
    
    if (tmp$convergence>0)
      { ## if still no convergence, keep past values
        method <- "NONE"
        tmp$par <- ph
      }

    list(p=as.vector(fn.lambdap(tmp$par)),
         method=as.character(method),
         time=as.numeric(proc.time()[3]-ptm))
  }
