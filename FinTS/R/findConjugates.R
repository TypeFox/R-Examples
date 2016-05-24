findConjugates <- function(x, complex.eps=.Machine[["double.eps"]]){
##
##  1.  compute normalization 
##
  if(length(x)<1)return(complex(0))
  ax <- abs(x)
  m2 <- outer(ax, ax, pmax)
##
##  2.  Compute complex differences 
##
  c2 <- (abs(outer(x, Conj(x), "-") / m2) < complex.eps)
  c2[m2==0] <- FALSE
  c2 <- (c2 & lower.tri(c2))
##
## 3.  Any differences exceed complex.eps?  
##
  if(any(c2)){
#     check standard differences
    d2 <- (abs(outer(x, x, "-") / m2) > complex.eps)
    d2[m2==0] <- FALSE
#
    cd2 <- (c2 & d2)
    if(any(cd2)){
      ic <- sort(unique(row(cd2)[cd2]))
      return(x[ic])
    }
  }
  complex(0)
}
