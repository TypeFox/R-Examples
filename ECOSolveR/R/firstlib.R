#' Return the default optimization parameters for ECOS
#'
#' This is used to control the behavior of the undelying optimization code.
#'
#' @param feastol the tolerance on the primal and dual residual, default 1e-8
#' @param abstol the absolute tolerance on the duality gap, default 1e-8
#' @param reltol the relative tolerance on the duality gap, default 1e-8
#' @param feastol_inacc the tolerance on the primal and dual residual if reduced precisions, default 1e-4
#' @param abstol_inacc the absolute tolerance on the duality gap if reduced precision, default 5e-5
#' @param reltol_inacc the relative tolerance on the duality gap if reduced precision, default 5e-5
#' @param maxit the maximum number of iterations for ecos, default 100L
#' @param mi_max_iters the maximum number of branch and bound iterations (mixed integer problems only), default 1000L
#' @param mi_int_tol the integer tolerence (mixed integer problems only), default 1e-4
#' @param mi_abs_eps the absolute tolerance between upper and lower bounds (mixed integer problems only), default 1e-6
#' @param mi_rel_eps the relative tolerance, \eqn{(U-L)/L}, between upper and lower bounds (mixed integer problems only), default 1e-6
#' @param verbose verbosity level, default 0L. A verbosity level of 1L will show more detail, but clutter session transcript.
#'
#' @return a list with the following elements:
#'  \describe{
#'    \item{FEASTOL}{ the tolerance on the primal and dual residual, parameter \code{feastol}}
#'    \item{ABSTOL}{ the absolute tolerance on the duality gap, parameter \code{abstol}}
#'    \item{RELTOL}{ the relative tolerance on the duality gap, parameter \code{reltol}}
#'    \item{FEASTOL_INACC}{ the tolerance on the primal and dual residual if reduced precisions, parameter \code{feastol_inacc}}
#'    \item{ABSTOL_INACC}{ the absolute tolerance on the duality gap if reduced precision, parameter \code{abstol_inacc}}
#'    \item{RELTOL_INACC}{ the relative tolerance on the duality gap if reduced precision, parameter \code{reltol_inacc}}
#'    \item{MAXIT}{ the maximum number of iterations for ecos, parameter \code{maxit}}
#'    \item{MI_MAX_ITERS}{ the maximum number of branch and bound iterations (mixed integer problems only), parameter \code{mi_max_iters}}
#'    \item{MI_INT_TOL}{ the integer tolerence (mixed integer problems only), parameter \code{mi_int_tol}}
#'    \item{MI_ABS_EPS}{ the absolute tolerance between upper and lower bounds (mixed integer problems only), parameter \code{mi_abs_eps}}
#'    \item{MI_REL_EPS}{ the relative tolerance, \eqn{(U-L)/L}, between upper and lower bounds (mixed integer problems only), parameter \code{mi_rel_eps}}
#'    \item{VERBOSE}{ verbosity level, parameter \code{verbose}}
#'   }
#'
#' @export
ecos.control <- function(maxit = 100L,
                         feastol = 1e-8,
                         reltol = 1e-8,
                         abstol = 1e-8,
                         feastol_inacc = 1e-4,
                         abstol_inacc = 5e-5,
                         reltol_inacc = 5e-5,
                         verbose = 0L,
                         mi_max_iters = 1000L,
                         mi_int_tol = 1e-4,
                         mi_abs_eps = 1e-6,
                         mi_rel_eps = 1e-6) {
    list(MAXIT = maxit,
         FEASTOL = feastol,
         RELTOL = reltol,
         ABSTOL = abstol,
         FEASTOL_INACC = feastol_inacc,
         ABSTOL_INACC = abstol_inacc,
         RELTOL_INACC = reltol_inacc,
         VERBOSE = verbose,
         MI_MAX_ITERS = mi_max_iters,
         MI_INT_TOL = mi_int_tol,
         MI_ABS_EPS = mi_abs_eps,
         MI_REL_EPS = mi_rel_eps)
}

isNonnegativeInt <- function(x) {
    (!is.null(x)) && (typeof(x) == "integer") && (length(x) == 1) && ! (x < 0)
}

isNonnegativeFloat <- function(x) {
    (!is.null(x)) && (typeof(x) == "double") && (length(x) == 1) && ! (x < 0)
}

isNontrivialNumericVector <- function(x) {
    (typeof(x) == "double") && (length(x) > 0)
}

isNontrivialIntegerVector <- function(x) {
    (typeof(x) == "integer") && (length(x) > 0)
}

checkOptions <- function(control)  {

    ## Check options
    if (!isNonnegativeInt(control$VERBOSE)) {
        return("Expected non-negative integer for option VERBOSE")
    }

    if (!isNonnegativeInt(control$MAXIT)) {
        return("Expected non-negative integer for option MAXIT")
    }

    if (!isNonnegativeFloat(control$FEASTOL)) {
        return("Expected non-negative numeric for option FEASTOL")
    }

    if (!isNonnegativeFloat(control$ABSTOL)) {
        return("Expected non-negative numeric for option ABSTOL")
    }

    if (!isNonnegativeFloat(control$RELTOL)) {
        return("Expected non-negative numeric for option RELTOL")
    }

    if (!isNonnegativeFloat(control$ABSTOL_INACC)) {
        return("Expected non-negative numeric for option ABSTOL_INACC")
    }

    if (!isNonnegativeFloat(control$FEASTOL_INACC)) {
        return("Expected non-negative numeric for option FEASTOL_INACC")
    }

    if (!isNonnegativeFloat(control$RELTOL_INACC)) {
        return("Expected non-negative numeric for option RELTOL_INACC")
    }

    if (!isNonnegativeInt(control$MI_MAX_ITERS)) {
        return("Expected non-negative integer for option MI_MAX_ITERS")
    }

    if (!isNonnegativeFloat(control$MI_ABS_EPS)) {
        return("Expected non-negative numeric for option MI_ABS_EPS")
    }

    if (!isNonnegativeFloat(control$MI_REL_EPS)) {
        return("Expected non-negative numeric for option MI_REL_EPS")
    }

    if (!isNonnegativeFloat(control$MI_INT_TOL)) {
        return("Expected non-negative numeric for option MI_INT_TOL")
    }
    ## else return NULL
    return(NULL)

}


.onLoad <- function(lib, pkg) {
    ## if (is.null(getOption("ecos.control"))) {
    ##     options(ecos.control = .ecos.control)
    ## }
}

.onUnload <- function(libpath) {
    library.dynam.unload("ECOSolveR", libpath)
}
