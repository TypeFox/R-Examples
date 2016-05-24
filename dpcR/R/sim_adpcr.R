#' Simulate Array Digital PCR
#' 
#' A function that simulates results of an array digital PCR.
#' 
#' The array digital PCR is performed on plates containing many microfluidic
#' chambers with a randomly distributed DNA template, fluorescence labels and
#' standard PCR reagents. After the amplification reaction, performed
#' independently in each chamber, the chambers with the fluorescence level
#' below certain threshold are treated as negative. From differences between
#' amplification curves of positive chambers it is possible to calculate both
#' total number of template molecules and their approximate number in a
#' single chamber.
#' 
#' The function contains two implementations of the array digital PCR
#' simulation. First one was described in Dube at. al (2008). This method is
#' based on random distributing \eqn{m \times times}{m * times} molecules
#' between \eqn{n \times times}{n * times} chambers.  After this step, the
#' required number of plates is created by the random sampling of chambers
#' without replacement. The above method is used, when the \code{dube} argument
#' has value \code{TRUE}.
#' 
#' The second method treats the total number of template molecules as random
#' variable with a normal distribution \eqn{\mathcal{N}(n, 0.05n)}{N(n,
#' 0.05*n)}. The exact sum of total molecules per plate is calculated and
#' randomly adjusted to the value of \eqn{m \times times}{m * times}. The above
#' method is used, when the \code{dube} argument has value \code{FALSE}. This
#' implementation is much faster than previous one, especially for big
#' simulations. The higher the value of the argument \code{times}, the
#' simulation result is closer to theoretical calculations.
#' 
#' @aliases sim_adpcr
#' @param m the total number of template molecules added to the plate. Must be
#' a positive integer.
#' @param n the number of chambers per plate. Must be a positive integer.
#' @param times number of repetitions (see Details).
#' @param n_panels the number of panels that are simulated by the function.
#' Cannot have higher value than the \code{times} argument.
#' @param dube if \code{TRUE}, the function is strict implementation of array
#' digital PCR simulation (as in Dube et al., 2008). If \code{FALSE}, the
#' function calculates only approximation of Dube's experiment. See Details and
#' References.
#' @param pos_sums if \code{TRUE}, function returns only the total number of
#' positive (containing at least one molecule) chamber per panel. If
#' \code{FALSE}, the functions returns a vector of length equal to the number
#' of chambers. Each element of the vector represents the number of template
#' molecules in a given chamber.
#' @return If the \code{pos_sums} argument has value \code{FALSE}, the function
#' returns a matrix with \eqn{n} rows and \eqn{n_panels} columns. Each column
#' represents one plate. The type of such simulation would be \code{"nm"}. If the
#' \code{pos_sums} argument has value \code{TRUE}, the function returns a
#' matrix with one row and \eqn{n_panels} columns. Each column contains the
#' total number of positive chambers in each plate and type of simulation would
#' be set as \code{"tnp"}.
#' 
#' In each case the value is an object of the \code{\linkS4class{adpcr}} class.
#' @author Michal Burdukiewicz.
#' @seealso \code{\link{sim_ddpcr}}.
#' @references Dube S, Qin J, Ramakrishnan R, \emph{Mathematical Analysis of
#' Copy Number Variation in a DNA Sample Using Digital PCR on a Nanofluidic
#' Device}.  PLoS ONE 3(8), 2008.
#' @keywords datagen
#' @examples
#' 
#' # Simulation of a digital PCR experiment with a chamber based technology.
#' # The parameter pos_sums was altered to change how the total number of positive 
#' # chamber per panel are returned. An alteration of the parameter has an impact 
#' # in the system performance.
#' adpcr_big <- sim_adpcr(m = 10, n = 40, times = 1000, pos_sums = FALSE, n_panels = 1000)
#' adpcr_small <- sim_adpcr(m = 10, n = 40, times = 1000, pos_sums = TRUE, n_panels = 1000)
#' # with pos_sums = TRUE, output allocates less memory
#' object.size(adpcr_big)
#' object.size(adpcr_small)
#' 
#' # Mini version of Dube et al. 2008 experiment, full requires replic <- 70000
#' # The number of replicates was reduced by a factor of 100 to lower the computation time.
#' replic <- 700
#' dube <- sim_adpcr(400, 765, times = replic, dube = TRUE, 
#'     pos_sums = TRUE, n_panels = replic)
#' mean(dube) # 311.5616
#' sd(dube) # 13.64159
#' 
#' # Create a barplot from the simulated data similar to Dube et al. 2008
#' bp <- barplot(table(factor(dube, levels = min(dube):max(dube))), 
#' 	      space = 0)
#' lines(bp, dnorm(min(dube):max(dube), mean = 311.5, sd = 13.59)*replic, 
#'       col = "green", lwd = 3)
#' 
#' # Exact Dube's method is a bit slower than other one, but more accurate
#' system.time(dub <- sim_adpcr(m = 400, n = 765, times = 500, n_panels = 500, 
#'   pos_sums = TRUE))
#' system.time(mul <- sim_adpcr(m = 400, n = 765, times = 500, n_panels = 500, 
#'   pos_sums = FALSE))
#' 
#' @export sim_adpcr
sim_adpcr <- function(m, n, times, n_panels = 1, dube = FALSE, pos_sums = FALSE) {
  n <- num2int(n)
  res <- sim_dpcr(m, n, times, dube, pos_sums, n_panels)
  create_adpcr(res, n = rep(n, n_panels), breaks = 0L:max(res), 
               type = ifelse(pos_sums, "tnp", "nm"))
}
