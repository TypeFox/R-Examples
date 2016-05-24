# ---- roxygen documentation ----
#
#' @title Proximity Index
#'
#' @description
#' The function \code{Prox} simply computes the proportion of fixes that are proximal, based on some spatial
#' threshold -- \code{dc} (Bertrand et al. 1996). It also facilitates local-level proximity analysis
#'
#' @details
#' The function \code{Prox} can be used to test for the presence of attraction (via proximity) in
#' wildlife telemetry data. Prox is simply the proportion of simultaneous fixes within
#' the threshold distance -- \code{dc}. The local output (dataframe) can be useful for examining variation
#' in proximity through time.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when 
#'    two fixes are spatially together.
#' @param local logical value indicating whether or not a dataframe, containing the
#'    distance between each simultaneous fix, should be returned.
#'
#' @return
#' If \code{local=FALSE} (the default) Prox returns the numeric value of the Prox index.
#' If \code{local=TRUE} Prox returns a dataframe (containing the date/times of \emph{all} simultaneous fixes, 
#' along with the distance between fixes at each time.
#'
#' @references
#' Bertrand, M.R., DeNicola, A.J., Beissinger, S.R, Swihart, R.K. (1996) Effects of parturition
#' on home ranges and social affiliations of female white-tailed deer.
#' \emph{Journal of Wildlife Management}, \bold{60}: 899-909.
#'
#' @keywords indices
#' @seealso GetSimultaneous, DI, contacts
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes, dc = 50 meters
#' Prox(deer37, deer38, tc=7.5*60, dc=50)
#' df <- Prox(deer37, deer38, tc=7.5*60, dc=50, local=TRUE)
#' 
#' @export
#
# ---- End of roxygen documentation ----
Prox <- function(traj1,traj2,tc=0,dc=50,local=FALSE){
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  n <- nrow(tr1)
  #Calculate the observed distances
  prox.df <- data.frame(date=tr1$date,prox=sqrt(((tr1$x - tr2$x)^2) + ((tr1$y - tr2$y)^2)))
  #compute the Prox index
  nprox <- length(which(prox.df$prox < dc))
  val <- nprox/n
  #Return a list of the Prox index value, and a dataframe for plotting if set to true
  if (local == TRUE){
    return(prox.df)
  }
  else {return(val)}
}
#================= end of Prox function ========================================