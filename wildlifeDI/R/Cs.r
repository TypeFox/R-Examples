# ---- roxygen documentation ----
#
#' @title Coefficient of Sociality
#'
#' @description
#'    The function \code{Cs} computes the coefficient of sociality between two moving
#'    objects following the methods outlined by Kenward et al. (1993). It also uses a
#'    signed Wilcoxon-rank test to test for significance.
#'
#' @details
#'  This function can be used to calculate the Kenward et al. (1993) coefficient of sociality (Cs)
#'  between two animals. The Cs statistic tests the observed mean
#'  distance between simultaneous fixes against that expected by the overall
#'  distribution of distances between all fixes.
#'  \deqn{Cs=\frac{D_E-D_O}{D_O+D_E}}{Cs=(De-Do)/(De+Do)}
#'  Where \eqn{D_O}{Do} is the mean observed distance between simultaneous fixes, and \eqn{D_E}{De}
#'  is the mean expected distance between all fixes. Kenward et al. (1993) propose Cs
#'  as a useful metric for exploring attraction or avoidance behaviour.
#'  Values for Cs closer to 1 indicate
#'  attraction, while values for Cs closer to -1 indicate avoidance. Values of Cs
#'  near 0 indicate that the two animals' movements have no influence on one another.
#'  \cr \cr
#'  Further, the difference between the observed and expected distances are compared
#'  using a paired signed-rank test (both one-sided tests, indicative of attraction
#'  or avoidance). See the function \code{GetSimultaneous} for details on how
#'  simultaneous fixes are determined from two trajectories.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#'  @param traj2 same as \code{traj1}.
#'  @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#'
#' @return
#'  This function returns a list of objects representing the calculated values from the
#'  Cs statistic and associated \emph{p}-values from the signed rank test.
#'  \itemize{
#'    \item Do -- The mean distance of simultaneous fixes.
#'    \item De -- The mean expected distance, from all fixes.
#'    \item Cs -- The coefficient of sociality, see \bold{Details}.
#'    \item p.Attract -- One sided \emph{p}-value from signed rank test, testing for attraction.
#'    \item p.Avoid -- One sided \emph{p}-value from signed rank test, testing for avoidance.
#'    }
#'
#' @references
#' Kenward, R.E., Marcstrom, V. and Karlbom, M. (1993) Post-nestling behaviour in
#'  goshawks, \emph{Accipiter gentilis: II}. Sex differences in sociality and nest-switching.
#'  \emph{Animal Behaviour}. \bold{46}, 371--378.
#'
#' @keywords indices
#' @seealso GetSimultaneous
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes
#' Cs(deer37, deer38, tc = 7.5*60) 
#' 
#' @export
#
# ---- End of roxygen documentation ----
Cs <- function(traj1,traj2, tc=0){
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  n <- nrow(tr1)

  #calculate the observed distances
  Do <- sqrt(((tr1$x - tr2$x)^2) + ((tr1$y - tr2$y)^2))


  #calculate the expected distances
  De <- rep(0,n)
  for (j in 1:n){
    De[j] <- sum(sqrt(((tr1$x[j] - tr2$x)^2) + ((tr1$y[j] - tr2$y)^2)))/n
    }

  #calculate the significance of differences b/w Do and De using a Wilcoxon signed rank test
  p.less <- wilcox.test(Do,De,paired=T,alternative="less")$p.value
  p.great <- wilcox.test(Do,De,paired=T,alternative="greater")$p.value

  #Compute the Coefficient of sociality (Kenward et al. 1993)
  DO <- sum(Do)/n
  DE <- sum(De)/n
  Cs <- (DE - DO) / (DO + DE)

  #return output
  output <- list(Do=DO, De=DE, Cs=Cs,p.Attract=p.less,p.Avoid=p.great)
  return(output)
  }


