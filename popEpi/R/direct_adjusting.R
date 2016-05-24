


#' @title Direct Adjusting in \pkg{popEpi} Using Weights
#' @author Joonas Miettinen
#' @name direct_standardization
#' @aliases direct_adjusting
#' @description 
#' 
#' Several functions in \pkg{popEpi} have support for direct standardization
#' of estimates. This document explains the usage of weighting with those
#' functions.
#' 
#' @details 
#' 
#' Direct standardization is performed by computing estimates of
#' \code{E}
#' by the set of adjusting variables \code{A}, to which a set of weights 
#' \code{W} is applicable. The weighted average over \code{A} is then the
#' direct-adjusted estimate of \code{E} (\code{E*}).
#' 
#' To enable both quick and easy as well as more rigorous usage of direct
#' standardization with weights, the weights arguments in \pkg{popEpi}
#' can be supplied in several ways. Ability to use the different
#' ways depends on the number of adjusting variables.
#' 
#' The weights are always handled internally to sum to 1, so they do not
#' need to be scaled in this manner when they are supplied. E.g.
#' counts of subjects in strata may be passed.
#' 
#' @section Basic usage - one adjusting variable:
#' 
#' In the simple case where we are adjusting by only one variable 
#' (e.g. by age group), one can simply supply a vector of weights:
#' 
#' \code{FUN(weights = c(0.1, 0.25, 0.25, 0.2, 0.2))}
#' 
#' which may be stored in advance:
#' 
#' \code{w <- c(0.1, 0.25, 0.25, 0.2, 0.2)}
#' 
#' \code{FUN(weights = w)}
#' 
#' The order of the weights matters. \pkg{popEpi} functions with direct
#' adjusting enabled match the supplied weights to the adjusting variables
#' as follows: If the adjusting variable is a \code{factor}, the order
#' of the levels is used. Otherwise, the alphabetic order of the unique
#' values is used (try \code{sort} to see how it works). For clarity
#' and certainty we recommend using \code{factor} or \code{numeric} variables
#' when possible. \code{character} variables should be avoided: to see why,
#' try \code{sort(15:9)} and \code{sort(as.character(15:9))}.
#' 
#' It is also possible to supply a \code{character} string corresponding
#' to one of the age group standardization schemes integrated into \pkg{popEpi}:
#' 
#' \itemize{
#' \item 'europe_1976_18of5' - european std. popupulation (1976), 18 age groups
#' \item 'nordic_2000_18of5' - nordic std. popupulation (2000), 18 age groups
#' \item 'world_1966_18of5' - world standard (1966), 18 age groups
#' \item 'world_2000_18of5' - world standard (2000), 18 agegroups
#' \item 'world_2000_20of5' - world standard (2000), 20 agegroups 
#' \item 'world_2000_101of1' - world standard (2000), 101 agegroups
#' }
#' 
#' You may also supply \code{weights = "internal"} to use internally
#' computed weights, i.e. usually simply the counts of subjects / person-time
#' experienced in each stratum. E.g.
#' 
#' \code{FUN(weights = "world_2000_18of5")}
#' 
#' will use the world standard population from 2000 as 
#' weights for 18 age groups, that your adjusting variable is 
#' assumed to contain. The adjusting variable must be coded in this case as 
#' a numeric variable containing \code{1:18} or as a \code{factor} with
#' 18 levels (coded from the youngest to the oldest age group). 
#' 
#' @section More than one adjusting variable:
#' 
#' In the case that you employ more than one adjusting variable, separate
#' weights should be passed to match to the levels of the different adjusting
#' variables. When supplied correctly, "grand" weights are formed based on
#' the variable-specific weights by multiplying over the variable-specific
#' weights (e.g. if men have \code{w = 0.5} and the age group 0-4 has 
#' \code{w = 0.1}, the "grand" weight for men aged 0-4 is \code{0.5*0.1}).
#' The "grand" weights are then used for adjusting after ensuring they
#' sum to one.
#' 
#' When using multiple adjusting variables, you
#' are allowed to pass either a named \code{list} of 
#' weights or a \code{data.frame} of weights. E.g.
#' 
#' \code{WL <- list(agegroup = age_w, sex = sex_w)}
#' 
#' \code{FUN(weights = WL)}
#' 
#' where \code{age_w} and \code{sex_w} are numeric vectors. Given the 
#' conditions explained in the previous section are satisfied, you may also do
#' e.g.
#' 
#' \code{WL <- list(agegroup = "world_2000_18of", sex = sex_w)}
#' 
#' \code{FUN(weights = WL)}
#' 
#' and the world standard pop is used as weights for the age groups as outlined
#' in the previous section.
#' 
#' Sometimes using a \code{data.frame} can be clearer (and it is fool-proof
#' as well). To do this, form a \code{data.frame} that repeats the levels
#' of your adjusting variables by each level of every other adjusting variable,
#' and assign the weights as a column named \code{"weights"}. E.g.
#' 
#' \code{wdf <- data.frame(sex = rep(0:1, each = 18), agegroup = rep(1:18, 2))}
#' 
#' \code{wdf$weights <- rbinom(36, size = 100, prob = 0.25)}
#' 
#' \code{FUN(weights = wdf)}
#' 
#' If you want to use the counts of subjects in strata as the weights,
#' one way to do this is by e.g.
#' 
#' \code{wdf <- as.data.frame(x$V1, x$V2, x$V3)}
#' \code{names(wdf) <- c("V1", "V2", "V3", "weights")}
#' 
#' @seealso \code{\link{flexible_arguments}}
#' 
#' @references
#' Source of the Nordic standard population in 5-year age groups 
#' (also contains European & 1966 world standards):
#' \url{http://www-dep.iarc.fr/NORDCAN/english/glossary.htm}
#' 
#' Source of the 1976 European standard population: 
#' 
#' Waterhouse, J.,Muir, C.S.,Correa, P.,Powell, J., eds (1976). 
#' Cancer Incidence in Five Continents, Vol. III. 
#' IARC Scientific Publications, No. 15, Lyon, IARC
#' 
#' A comparison of the 1966 vs. 2000 world standard populations in 5-year age groups:
#' \url{http://www3.ha.org.hk/cancereg/e_asr.asp}
#' 
#' Source of 2000 world standard population in 1-year age groups:
#' \url{http://seer.cancer.gov/stdpopulations/stdpop.singleages.html}

NULL









