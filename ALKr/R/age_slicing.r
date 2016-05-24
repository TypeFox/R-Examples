#' Age slicing
#' 
#' Estimation of a population's age distribution from its length distribution
#' based on von Bertalanffy's growth curve, as described by Kell and Kell (2011)
#' 
#' Age distribution is calculated by using the inverse of the von Bertalanffy
#' growth curve, whose parameters must be known. Limits for the minimum and
#' maximum ages must also be given. 
#'
#' @param fi A vector with the number of individuals per length class.
#' @param li A vector containing the length value for each length class. If
#' blank, will use \code{as.numeric(names(fi))}.
#' @param vb_params A named vector with the parameters of the von Bertalanffy
#' growth equation, \code{Linf}, \code{K} and \code{t0}.
#' @param age_limits A vector with two elements, containing the lowest and
#' highest age classes.
#' @param timing Correction for the offset between the data collection and
#' recruitment. Defaults to \code{0.5}, i.e. half a year
#' 
#' @return A vector containing the number of individuals in each age class.
#' 
#' @references
#' Kell, L., Kell, A. (2011). A comparison of age slicing and statistical age
#' estimation for mediterranean swordfish (Xiphias gladius). \emph{Collect. Vol.
#' Sci. Pap. ICCAT}. \strong{66}/4, 1522-1534
#' 
#' @examples
#' data(hom)
#' age_slicing(fi = hom$F1992,
#'   vb_params = c(Linf = 54.98, K = 0.064, t0 = -4.68),
#'   age_limits = c(0,5))
#' 
#' @export
#' 
age_slicing <- function(fi, li = as.numeric(names(fi)), vb_params, age_limits,
                        timing = 0.5) {
  
  if (!all(names(vb_params) %in% c("t0", "K", "Linf")))
    stop("vb_params must contain 't0', 'K', 'Linf'")
  if (!length(age_limits) == 2)
    stop("age_limits must be a vector with two elements")
    
  age <- floor(vb_params["t0"] - log(1 - pmin(li / vb_params["Linf"], .999999999)) / vb_params["K"] + timing)
  age <- pmax(pmin(age, age_limits[2]), age_limits[1]) # age can't be lower than 0 nor greater that last_age
  age <- factor(age, levels = age_limits[1]:age_limits[2])
  
  tapply(fi, age, sum)
}

