## Function which computes mu and Sigma from the largest weights
## __input__
## theta      : [Nxk matrix] of draws
## w          : [Nx1 vector] of weights
## control    : [list] containing the following components:
## $ISpercent : [k1x1 vector] of percentages
## $ISscale   : [k2x1 vector] of scaling factors
## __output__
## [list] with the following components:
## mu         : [kx1 vector] of location
## Sigma      : [(k1*k2)xk^2 matrix] of scale matrics (in rows and in vector form)
## method     : [(k1*k2)x1 character vector] of IS computations
## time       : [double] computation time
## __20080429__
'fn.wIS' <- function(theta, w, control)
  {
    ptm <- proc.time()[3]
    k <- ncol(theta)
    n <- nrow(theta)
    pos <- order(w, decreasing=TRUE)
    Sigma <- method <- NULL
    for (i in control$percent)
      { ## iterate over percentages of importance weights
        tmppos <- pos[1:floor(n*i)]
        tmp <- fn.muSigma(w[tmppos], theta[tmppos,], mu=theta[pos[1],])
        for (j in control$scale)
          { ## iterate over scaling factors
            Sigma <- rbind(Sigma, j*tmp$Sigma)
            method <- c(method, paste("IS", paste(i,j,sep="-")))
          }
      }
    
    list(mu=tmp$mu,
         Sigma=matrix(Sigma, ncol=k*k),
         method=as.character(method),
         time=as.numeric(proc.time()[3]-ptm))
  }

         
