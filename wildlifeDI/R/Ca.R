# ---- roxygen documentation ----
#
#' @title Coefficient of Association
#'
#' @description
#'  This function measures the dynamic interaction between two moving objects following
#'  the methods first described by Cole (1949), and more recently employed by Bauman (1998).
#'
#' @details
#'  This function can be used to calculate the Cole (1949) measure of dynamic
#'  interaction between two animals. Termed a coefficient of association, the Ca
#'  statistic tests the number of fixes the animals are observed together against the
#'  total number of fixes following:
#'  \deqn{Ca = \frac{2AB}{A+B}}{2AB/(A+B)}
#'  where \eqn{A} (respectively \eqn{B}) is the number of times animal 1 (resp. 2) are
#'  observed, and \eqn{AB} is the number of times the two animals are observed together.
#'  Several works, including Bauman (1998) have suggested that Ca > 0.5 indicates
#'  affiliation or fidelity, while Ca < 0.5 indicates no association between the
#'  two animals. Note that this function calls \code{GetSimultaneous} to identify the temporal
#'  component of identifying when fixes together.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped 
#'         movement fixes of the first object. Note this object must be a \code{type II 
#'         ltraj} object. For more information on objects of this type see \code{
#'         help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc temporal tolerance limit (in seconds) for defining when two fixes
#'         are simultaneous or together. Parameter passed to function \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when 
#'         two fixes are spatially together.
#'
#' @return
#'  This function returns a numeric result of the Ca statistic.
#'
#' @references
#'  Bauman, P.J. (1998) The Wind Cave National Park elk herd: home ranges, seasonal movements, and alternative control methods.
#'    M.S. Thesis. South Dakota State University, Brookings, South Dakota, USA. \cr\cr
#'  Cole, L.C. (1949) The measurement of interspecific association. \emph{Ecology}. \bold{30}, 411--424.
#'
#' @keywords indices
#' @seealso GetSimultaneous, Prox, HAI
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes, dc = 50 meters
#' Ca(deer37, deer38, tc = 7.5*60, dc = 50)
#' 
#' @export
#
# ---- End of roxygen documentation ----

Ca <- function(traj1,traj2,tc=0,dc=50){
  #convert ltraj objects to dataframes
  A <- dim(ld(traj1))[1]
  B <- dim(ld(traj2))[1]
  
  trajs <- GetSimultaneous(traj1,traj2,tc)
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  
  trDist <- sqrt(((tr1$x - tr2$x)^2) + ((tr1$y - tr2$y)^2))
  
  AB <- length(which(trDist <= dc))
  
  #Compute coefficent of association
  Ca <- 2*AB/(A+B)
  return(Ca)
}
#=============== End of Ca Function =======================================
