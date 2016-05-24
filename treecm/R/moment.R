#' @title Returns the coordinates of centre of mass of branches and logs
#'
#' @description Computes the cartesian coordinates of centre of mass of branches and 
#' logs along with their \eqn{x}, \eqn{y}, \eqn{z} moments
#'
#' The \eqn{x} and \eqn{y} coordinates are computed from the polar coordinates (angle and distance, 
#' defined as the length of its projection on ground), measured in the field. 
#' The \eqn{z} coordinate is computed by adding the height of branch insertion on the stem 
#' (measured in the field) to the height of the branch (calculated through its 
#' mean tilt, in case it was measured in the filed).
#' The \eqn{x}, \eqn{y}, \eqn{z} coordinates are corrected to take into account where the actual
#' centre of mass lies on the branches themselves by multiplying them by branchesCM,
#' a real number from 0.01 (CM at branch base) to 1.00 (CM at branch tip). 
#' As a rule of thumb, average live branches, with an average amount of foliage, 
#' have CM approx. \eqn{2/3} of their length, ie. branchesCM = 0.66.
#' \eqn{x}, \eqn{y}, \eqn{z} moments are computed by multiplying the corresponding cartesian coordinates by 
#' branch or log mass, e.g. \eqn{m_x = F \cdot l_x}, where \eqn{F} is branch or log mass, \eqn{l_x} is the \eqn{x} component of the lever arm (e.g. the \eqn{x} component of the branch or log projection on the ground).
#'
#' @param azimuth Branch compass heading
#' @param dBase unused argument
#' @param dTip unused argument
#' @param length Branch length
#' @param tipD unused argument
#' @param height Height of branch insertion on the stem or the height of log lower section
#' @param tilt Inclination of the branch or log in degrees
#' @param toBePruned unused argument
#' @param biomass Mass of the branch or log
#' @param branchesCM a real number varying from 0.01 to 1 proportional to the centre of
#' mass position along the branch (0.01 branch base, 1 branch tip)
#' @return a vector holding 5 reals:
#' \itemize{
#'  \item{the \eqn{x} coordinate of branch CM}
#'  \item{the \eqn{y} coordinate of branch CM}
#'  \item{the \eqn{x} moment of the branch}
#'  \item{the \eqn{y} moment of the branch}
#'  \item{the \eqn{z} moment of the branch}}
#' @note BranchCM is assumed to have same value in branches and logs. This is not the case in the real world. As a measure of safety one should use the highest value possible, eg branchesCM = 1.
#'
#' @note \eqn{z} coordinate of CM is not returned because it would be useless in a 2D plot. It is computed using \eqn{mz}, which is, as a matter of facts, returned
getCoordinatesAndMoment <- function (azimuth, dBase, dTip, length, tipD, height, tilt, toBePruned, biomass, branchesCM) {
  # height (h) to be added to branch height (z), as a function of the 
  # angle of its tilt (0 degrees = horiz., 90 degrees = vert.), its distance (length of its 
  # projection on the ground, 
  # from tree base to branch tip), and the estimated position of the centre of mass
  h  <- tipD * sin(tilt * pi / 180) * branchesCM
  ## computes cartesian coordinates of centre of mass of branches and their moments (mx, my, mz).
  ## When branchesCM = 1, x and y are coordinates of branch tip
  xy  <- toCartesianXY(azimuth, (tipD * branchesCM))
  z   <- height + h
  
  mx  <- biomass * xy[1]
  my  <- biomass * xy[2]
  mz  <- biomass * z
  c(xy, mx, my, mz)
}
