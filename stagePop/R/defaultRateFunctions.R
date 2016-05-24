#' @include reproductionFunction.R
#' @include deathFunction.R
#' @include durationFunction.R
#' @include immigrationFunction.R
#' @include developmentFunction.R
#' @include emigrationFunction.R
NULL

#' defaultRateFunctions
#'
#' These default implementations will simply generate errors when
#' run. To create implementations, please see the documentation linked
#' to below.
#'
#' The list should contain the following names, each mapped to a function of the correct signature.
#'
#' \itemize{
#' \item \code{reproFunc}
#' \item \code{deathFunc}
#' \item \code{durationFunc}
#' \item \code{immigrationFunc}
#' \item \code{develFunc}. Note that by defaul develFunc is NULL as it is not required for all simulation types
#' }
#'
#' @seealso \code{\link{RateFunctions}} 
defaultRateFunctions <- list(
                             reproFunc=reproFuncDefault,
                                                          
                             deathFunc=deathFuncDefault,                     
                             
                             durationFunc=durationFuncDefault,

                             immigrationFunc=immigrationFuncDefault,
                             
                             develFunc=NULL,

                             emigrationFunc=emigrationFuncDefault
                             )
