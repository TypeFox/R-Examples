#' Fitting routines for the Vitality family of mortality models.
#'
#' This package provides support for fitting models in the vitality family of mortality models.
#' Currently, the 4-parameter and 6paramter 2-process models are included, but planned updates will
#' include all published versions of the models, as well as Bayesian parameter estimation routines.
#' 
#' Support developing this package  was provided to J. Anderson by Bonneville Power Administration
#' and the University of Washington Center for Statistics and the Social Sciences
#' and to G. Passolt by the University of Washington Center for studies in Demography and Ecology.
#' 
#' @examples \dontrun{
#' data(swedish_females)
#' head(swedish_females)
#' initial_age <- 0 # (Could be adjusted up)
#' time <- initial_age:max(swedish_females$age)
#' survival_fraction <- swe$lx / swe$lx[swe$age == initial_age]
#' sample_size <- swe$Lx[swe$age == initial_age] #sample size
#' results <- vitality.2ps(time = time,
#'     sdata = survival_fraction,
#'     init.params=c(0.012, 0.01, 0.1, 0.1),
#'     se = sample_size,
#'     Mplot=F)
#' }
#' 
#' @references
#' \itemize{
#'   \item Li, T. and J.J. Anderson (in press).
#'   "Shaping human mortality patterns through intrinsic and extrinsiv vitality processes."
#'   Demographic Research.
#'   \item Salinger, D.H., J.J. Anderson, and O.S. Hamel. 2003.
#'   "A parameter estimation routine for the vitality-based survival model."
#'   Ecological Modelling 166 (3): 287-29
#'   \item Li, T. and J.J. Anderson. 2009.
#'   "The vitality model: A way to understand population survival and demographic heterogeneity."
#'   Theoretical Population Biology 76: 118-131.
#'   \item Anderson, J.J., Molly C. Gildea, Drew W. Williams, and Ting Li. 2008.
#'   "Linking Growth, Survival, and Heterogeneity through Vitality.
#'   The American Naturalist 171 (1): E20-E43.
#'   }
#'   
#' @import IMIS
#' @docType package
#' @name vitality
NULL
