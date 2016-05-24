## Function which generates draws from a mixture of Student-t densities
## __input__
## N        : [integer>0] number of draws (default: 1)
## mit      : [list] of mixture information (default: univariate Cauchy)
##  $p      : [Hx1 vector] of probabilities
##  $mu     : [Hxk matrix] of mean vectors (in row)
##  $Sigma  : [Hxk^2 matrix] of scale matrices (in row)
##  $df     : [integer>0] degrees of freedom parameter
## __output__
## [Nxk matrix] of draws
## __20080429__
## changes marked by # date 20120816 (updated mit definition: mit$df can be a double or vector)
'rMit' <- function(N=1, mit=list())
{
  H <- length(mit$p)
  if (H==0)
  { ## default
    warning ("'mit' not well defined in 'rMit'; set to default")
    mit <- list(p=1, mu=as.matrix(0), Sigma=as.matrix(1), df=1)
    H <- 1
  }
  comp <- sample(1:H, N, prob=mit$p, replace=TRUE)
  r <- NULL
  # date 20120816: df can be a vector of size 1 or H, replicated in the former case
  if(length(mit$df)==1 & H >1)
    mit$df = rep(mit$df,H)
  for (h in 1:H) {
    nh <- length(comp[comp == h])
    tmp <- NULL
    if (nh > 0) 
      # date 20120816: simulation from mixture with possibly different df parameters
      #        tmp <- fn.rmvt(nh, mit$mu[h, ], mit$Sigma[h, ], mit$df)
      tmp <- fn.rmvt(nh, mit$mu[h, ], mit$Sigma[h, ], mit$df[h])
    r <- rbind(r, tmp)
  }
  r <- as.matrix(r[sample(1:N), ])
  if (ncol(r) == 1) 
    r <- as.vector(r)
  return(r)
}



