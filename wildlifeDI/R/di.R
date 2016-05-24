# ---- roxygen documentation ----
#
#' @title Dynamic interaction index
#'
#' @description
#' The function \code{DI} measures dynamic interaction between two moving objects. It
#' calculates the local level di statistic for movement displacement, direction, 
#' and overall. DI can compute time- and/or distance-based weighting schemes
#' following Long and Nelson (2013).
#' 
#' @details
#' This function can be used for calculating the dynamic interaction (DI) statistic 
#' as described in Long and Nelson (2013). The DI statistic can be used to measure 
#' the local level of dynamic interaction between two moving objects. Specifically, 
#' it measures dynamic interaction in movement direction and displacement. 
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped 
#' movement fixes of the first object. Note this object must be a \code{type II 
#' traj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param alpha value for the \eqn{\alpha} parameter in the formula for di\eqn{_d}.
#' @param TimeWeight whether or not to compute the weights for time-based weighting.
#' @param DistWeight whether or not to compute the weights for distance-based weighting.
#' @param local logical, whether or not to perform local analysis.
#'
#' @return
#' The function DI returns the DI index (along with DI\eqn{_{theta}}{_theta} and DI\eqn{_d}), or if \code{local = TRUE}, 
#' a dataframe that contains the localized \code{di} values (see Long and Nelson 2013). The columns for \code{di},
#' \code{di.theta}, and \code{di.d} representdynamic interaction overall, in 
#' direction (azimuth), and in displacement respectively. If time- and/or distance-based weighting 
#' is selected, the corresponding weights are included as additional columns \code{t.weight} and \code{d.weight},
#' respectively. Please see Long and Nelson (2013) for a more detailed description.
#'
#' @references
#' Long, J.A., Nelson, T.A. 2013. Measuring dynamic interaction in movement
#' data. \emph{Transactions in GIS}. 17(1): 62-77.
#'
#' @keywords indices
#' @seealso GetSimultaneous, Cr
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes
#' DI(deer37, deer38, tc = 7.5*60)
#' df <- DI(deer37, deer38, tc = 7.5*60, local = TRUE)
#' 
#' @export
#
# ---- End of roxygen documentation ----
DI <- function(traj1, traj2, tc = 0, alpha = 1, local=FALSE, TimeWeight = FALSE, DistWeight = FALSE){
  #convert ltraj objects to dataframes
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  n <- nrow(tr1)
  
  #create some empty vectors
  f.theta <- rep(0,n-1)
  g.displ <- rep(0,n-1)
  
  #loop through the records
  for (i in 1:(n-1)){
    #compute interaction in azimuth
    if (is.na(tr1$abs.angle[i]) == FALSE){
      if (is.na(tr2$abs.angle[i])== FALSE){
        f.theta[i] <- cos(tr1$abs.angle[i] - tr2$abs.angle[i])
      } else {f.theta[i] <- 0}
    } else {
      if (is.na(tr2$abs.angle[i]) == TRUE){
        f.theta[i] <- 1
      } else {f.theta[i] <- 0}
    }
    
    #compute interaction in displacement
    if (tr1$dist[i] + tr2$dist[i] > 0){
      g.displ[i] <- 1 - (abs(tr1$dist[i] - tr2$dist[i])/(tr1$dist[i] + tr2$dist[i]))^alpha
    } else {g.displ[i] <- 1}
  }
  #compute overall interaction
  di.total <- f.theta*g.displ
  
  outdf <- data.frame(date = tr1$date, di = c(di.total,NA), di.theta = c(f.theta,NA), di.d = c(g.displ,NA))
  
  weights <- rep(1,(n-1))
  #Time Weighting
  if (TimeWeight == TRUE){
    TT <- sum(tr1$dt[which(!is.na(tr1$dt))])
    tWeights <- (tr1$dt[which(!is.na(tr1$dt))] / TT)*(n-1)
    weights <- weights*tWeights
    outdf$t.Weight <- c(tWeights,NA)
  }
  
  #Distance Weighting
  if (DistWeight == TRUE){
    #I have implemented the average as the distance based weight, other forms
    # might be useful and could be implemented here
    d.avg <- (tr1$dist[which(!is.na(tr1$dist))] + tr2$dist[which(!is.na(tr2$dist))])/2
    #-------------------------------------------------------------------------
    DD <- sum(d.avg)
    dWeights <- (d.avg / DD)*(n-1)
    weights <- weights*dWeights
    outdf$d.weight <- c(dWeights,NA)
  } 
  
  DI.out <- list(DI = mean(di.total*weights), DI.theta = mean(f.theta*weights), DI.d = mean(g.displ*weights))
  
  if (local == TRUE) {return(outdf)}
  else {return(DI.out)}
}
#==================== End of di Function =======================================
