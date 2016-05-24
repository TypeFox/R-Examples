## Function which computes mu and Sigma from the residuals weights
## See Hoogerheide (2006,Sect.2.3.1,p.48)
## __input__
## theta    : [Nxk matrix] of draws
## w        : [Nx1 vector] of weights
## ISctild  : [double] tuning parameter
## ISmax    : [integer>0] maximum iteration allowed
## ISfactor : [double] factor by which ctild is multiplied
## __output__
## [list] with the following components:
## $mu      : [kx1 vector] of locations
## $Sigma   : [k^2x1 vector] of scale matrices (in vector form)
## $time    : [double] time of the optimization
## __20080429__
'fn.wRes' <- function(theta, w, control)
  {
    ptm <- proc.time()[3]
    k <- ncol(theta)
    N <- nrow(theta)
    ctild <- 100*mean(w)
    ISstop <- FALSE
    while (ISstop==FALSE)
      { 
        wres <- w-ctild
        wres[wres<0] <- 0
        wrestild <- wres / sum(wres)
        wrestild <- matrix(wrestild, N, k)
        mu <- apply( wrestild * theta, 2, sum)
        tmp <- theta-matrix(mu, N, k, byrow=TRUE)
        tmpSigma <- t(tmp) %*% (wrestild*tmp)

        if (any(is.na(tmpSigma)) | any(is.nan(tmpSigma)))
          { ## if 'NA' of 'NaN' detected, scale ctild
            ctild <- .5*ctild
          }
        else if (!fn.isPD(tmpSigma))
          { ## if the matrix is not positive definite, scale ctild
            ctild <- .5*ctild
          }
        else
          { ## otherwise, stop the iteration
            ISstop <- TRUE
          }
      }

    Sigma <- NULL
    for (j in control$scale)
      { ## iterate of the scaling factors
        Sigma <- rbind(Sigma, as.vector(j*tmpSigma))
      }
    
    list(mu=as.vector(mu),
         Sigma=matrix(Sigma, ncol=k*k, byrow=TRUE),
         time=as.numeric(proc.time()[3]-ptm))
  }
