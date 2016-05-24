# 'Robust' sampling from mixture of (multivariate) t densities
# samples from 'mit' are re-drawn if they correspond to a bad parameter region with '0' kernel density
#
# input: 
#  robustify : [logical] if TRUE, gets N draws from non-zero part of the kernel,
#              if FALSE, gets N draws in total (effectively, may lead to less draws 
#              than N in the relevant region)
#  N         : [double] number of draws
#  mit       : [list] mixture density with H components, each component is a d-variate t dens.
#              Following components must be provided:
#    p       : [H-vector] probability of each component
#    mu      : [Hxd matrix] modes (in row) of each component
#    Sigma   : [Hx(d*d)  matrix] scales (in row, vectorized) of each component
#    df      : [H-vector] degree of freedom for each component
#  KERNEL    : [function] which computes the kernel. Must be vectorized
#  ...       : additional parameters used by 'KERNEL'
# outputs    : [list] with the following components:
#  theta     : [Nxd] matrix with draws in each row
#  lnk       : [N-vector] of log-kernel evaluations at draws
#
# author : Nalan Basturk
# date   : 20120912
fn.rmvgt_robust <- function(robustify=TRUE,N,mit,KERNEL,...){
   theta <- rmvgt(N,mit)    
   lnk = KERNEL(theta,...)
   # get more draws from relevant regions, if robustness is specified
   if(robustify){
     while(any(lnk == -Inf)){
        ind.remove        <- which(lnk == -Inf)
      draw.new          <- rmvgt(length(ind.remove),mit=mit)
      theta[ind.remove,] =  draw.new 
      lnk[ind.remove]   <- KERNEL(draw.new,...)
     }
   }
   return(list(theta=theta,lnk=lnk))
}
