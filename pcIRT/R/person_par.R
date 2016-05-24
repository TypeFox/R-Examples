#' Estimation of person parameters
#'
#' This function performs the estimation of person parameters for the
#' multidimensional polytomous Rasch model or the continuous Rating Scale
#' model.
#'
#' The estimation is performed by Maximum Likelihood Estimation. Thus,
#' parameters for extreme scores are not calculated!
#'
#' @aliases person_par person_par.MPRM person_par.CRSM
#' @param object Object of class \code{MPRM} or \code{CRSM}
#' @param \dots \dots{}
#' @return \item{ptable}{table showing for each (observed) raw score the
#' corresponding estimated person parameter and standard error}
#' \item{pparList}{for each person raw score, estimated person parameter and
#' the standard error is displayed} \item{fun_calls}{number of function calls}
#' \item{call}{function call}
#' @author Christine Hohensinn
#' @seealso \code{\link{CRSM}}
#' @references Fischer, G. H. (1974). Einfuehrung in die Theorie
#' psychologischer Tests [Introduction to test theory]. Bern: Huber.
#'
#' Mueller, H. (1999). Probabilistische Testmodelle fuer diskrete und
#' kontinuierliche Ratingskalen. [Probabilistic models for discrete and
#' continuous rating scales]. Bern: Huber.
#' @keywords person parameter
#'
#' @rdname perspar
#'
#' @examples
#'
#' #estimate CRSM for the first four items
#' data(analog)
#' res_cr <- CRSM(extraversion, low=-10, high=10)
#'
#' #estimate person parameters for CRSM
#' pp <- person_par(res_cr)
#'
#'
#' @export person_par
person_par <-
function(object,...)UseMethod("person_par")
