#' Dispersion Test for Spatial Point Pattern in Array dPCR Based on Quadrat
#' Counts
#' 
#' Performs a test of Complete Spatial Randomness for each
#' plate. This function is a wrapper around \code{\link[spatstat]{quadrat.test}}
#' function working directly on the objects of \code{\linkS4class{adpcr}}.
#' 
#' @details 
#' Under optimal conditions, the point pattern of dPCR events (e.g., positive droplet 
#' & negative droplets) should be randomly distrubuted over a planar chip. 
#' This function verifies this assumption using chi-square or Monte Carlo test.
#' Arrays with non-random patterns should be checked for integrity.
#' 
#' @param X Object of the \code{\linkS4class{adpcr}} class containing data from
#' one or more panels.
#' @param nx Number of quadrats in the x direction.
#' @param ny Number of quadrats in the y direction.
#' @param alternative \code{character} string (partially matched) specifying the
#' alternative hypothesis.
#' @param method \code{character} string (partially matched) specifying the test to
#' use: either \code{"Chisq"} for the chi-squared test (the default), or
#' \code{"MonteCarlo"} for a Monte Carlo test.
#' @param conditional \code{logical}. Should the Monte Carlo test be conducted
#' conditionally upon the observed number of points of the pattern? Ignored if
#' method="Chisq".
#' @param nsim The number of simulated samples to generate when
#' method="MonteCarlo".
#' @return A \code{list} of objects of class \code{"htest"} with the length equal to the
#' number of plates (minimum 1).
#' @note A similar result can be achived by using \code{\link{adpcr2ppp}} and
#' \code{\link[spatstat]{quadrat.test}}. See Examples.
#' @author Adrian Baddeley, Rolf Turner, Michal Burdukiewcz, Stefan Roediger.
#' @seealso \code{\link[spatstat]{quadrat.test}}.
#' @references http://www.spatstat.org/
#' @keywords pattern quadrat spatial dPCR
#' @examples
#' 
#' many_panels <- sim_adpcr(m = 400, n = 765, times = 1000, pos_sums = FALSE, 
#'                    n_panels = 5)
#' test_panel(many_panels)
#' 
#' #test only one plate
#' test_panel(extract_dpcr(many_panels, 3))
#' 
#' #do test_panel manually
#' require(spatstat)
#' ppp_data <- adpcr2ppp(many_panels)
#' lapply(ppp_data, function(single_panel) quadrat.test(single_panel))
#' 
#' 
#' @export test_panel
test_panel <- function(X, nx = 5, ny = 5, alternative = c("two.sided", "regular", "clustered"), 
                       method = c("Chisq", "MonteCarlo"), conditional = TRUE, nsim = 1999) {
  ppp_data <- adpcr2ppp(X)
  lapply(ppp_data, function(single_panel)
    quadrat.test(single_panel, nx, ny, alternative, method, conditional, nsim = 1999))
}
