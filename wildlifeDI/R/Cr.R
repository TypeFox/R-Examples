# ---- roxygen documentation ----
#
#' @title Movement Correlation Coefficient
#'
#' @description
#'   The function \code{Cr} computes the correlation statistic for movement data as presented 
#'   in the paper by Shirabe (2006). The statistic is essentially a Pearson product-moment 
#'   correlation statistic formulated for use with movement data.
#'   
#' @details
#'   The function \code{Cr} can be used to measure the level of dynamic interaction (termed correlation)
#'   between a pair of simultaneously moving objects. The statistic is sensitive to
#'   interaction in both movement direction (azimuth) and displacement, but is unable to
#'   disentangle the effects of these components. 
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#'
#' @return
#'   This function returns the Shirabe (2006) correlation statistic for two moving objects.
#'
#' @references
#' Shirabe, T. 2006. Correlation analysis of discrete motions. In: Raubal, M.,
#' Miller, HJ, Frank, AU, and Goodchild, M. eds. GIScience 2006, LNCS 4197. Berlin: Springer-Verlag;
#' 370-382.
#'
#' @keywords indices
#' @seealso GetSimultaneous, DI
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes
#' Cr(deer37, deer38, tc = 7.5*60)
#' 
#' @export
#
# ---- End of roxygen documentation ----
Cr <- function(traj1, traj2, tc = 0){
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])[,1:2]
  tr2 <- ld(trajs[2])[,1:2]
  
  n <- nrow(tr1)
  
  #compute vectors for each movement segment
  for (j in n:2)
  {
    tr1[j,] <- tr1[j,] - tr1[j-1,]
    tr2[j,] <- tr2[j,] - tr2[j-1,]
  }
  #Remove the beginning fix
  tr1 <- tr1[2:n,]
  tr2 <- tr2[2:n,]
  
  #compute path means for x,y coords
  tr1.bar <- apply(tr1,2,mean)
  tr2.bar <- apply(tr2,2,mean)
  #compute mean subtracted data matrix
  tr1.m <- t(t(as.matrix(tr1)) - tr1.bar)
  tr2.m <- t(t(as.matrix(tr2)) - tr2.bar)
  #function for numerator
  numer <- function(a,b)
  {
    numer <- sum(a*b)
    return(numer)
  }
  #compute numerator
  R.numer <- rep(0,(n-1))
  for (i in 1:(n-1))
  {
    R.numer[i] <- numer(tr1.m[i,],tr2.m[i,])
  }
  #compute denominator
  len <- function(v)
  {
    l <- sqrt(sum(v*v))
    return(l)
  }
  tr1.denom <- sqrt(sum(apply(tr1.m,1,len)^2))
  tr2.denom <- sqrt(sum(apply(tr2.m,1,len)^2))
  R.denom <- tr1.denom*tr2.denom
  #compute and return the statistic
  R <- sum(R.numer) / R.denom
  return(R)
}
#======================== End of shirabe Function ==============================