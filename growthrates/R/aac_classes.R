### ============================================================================
###  The S4 classes
### ============================================================================

### ----------------------------------------------------------------------------
### Helper Classes
### ----------------------------------------------------------------------------

#' @include aaa_growthmodel-class.R aab_growthmodel-constructor.R

setOldClass(c("modFit", "summary.modFit")) # classes from package FME

setOldClass("smooth.spline")


### ----------------------------------------------------------------------------
### S4 Classes Representing Single Fits
### ----------------------------------------------------------------------------

#' S4 Classes of Package \pkg{growthrates}
#'
#'
#' \code{growthrates_fit}: top-level class representing a growthrates fit with
#' any method.
#'
#' @slot FUN model function used.
#' @slot fit results of the model fit.
#' @slot obs observation data used for model fitting.
#' @slot rsquared coefficient of determination.
#'
#' @rdname growthrates-classes
#' @exportClass growthrates_fit
#'
setClass("growthrates_fit",
         representation(
           FUN = "function_growthmodel",
           fit = "ANY",
           obs =  "data.frame",
           rsquared = "numeric"
        )
)


#' \code{nonlinear_fit}: single nonlinear growthrates fit with package FME.
#'
#' @slot par parameters of the fit.
#' @rdname growthrates-classes
#' @exportClass nonlinear_fit
#'
setClass("nonlinear_fit",
         representation(
           fit = "modFit",
           par = "numeric"   ## fitted and fixed parms
         ),
         contains = "growthrates_fit"
)


#' \code{easylinear_fit}: single fit from the ``growthrates made easy''-method.
#'
#' @slot ndx index values of the time points used (for \code{easylinear_fit}).
#' @rdname growthrates-classes
#' @exportClass easylinear_fit
#'
setClass("easylinear_fit",
         representation(
           fit = "lm",
           par = "numeric",
           ndx = "numeric"
         ),
         contains = "growthrates_fit"
)

#' \code{spline_fit}: single fit with (optionally cross-validated) smoothing splines.
#'
#' @slot xy x and y values at the maximum of the spline.
#' @rdname growthrates-classes
#' @exportClass smooth.spline_fit
#'
setClass("smooth.spline_fit",
         representation(
           fit = "smooth.spline",
           xy  = "numeric",
           par = "numeric"
         ),
         contains = "growthrates_fit"
)



### ----------------------------------------------------------------------------
### S4 Classes Representing Multiple Fits
### ----------------------------------------------------------------------------


#' \code{multiple_fits}: top-level class representing multiple fits with
#' any method.
#'
#' @rdname growthrates-classes
#' @aliases multiple_fits
#' @exportClass multiple_fits
#'
setClass("multiple_fits",
         representation(
           fits = "list",
           grouping = "character"
         )
)


#' \code{multiple_easylinear_fits}: class representing multiple fits with
#' the ``growthrates made easy''-method.
#'
#' @rdname growthrates-classes
#' @exportClass multiple_easylinear_fits
#'
setClass("multiple_easylinear_fits",
         contains = "multiple_fits"
)


#' \code{multiple_nonlinear_fits}: class representing multiple nonlinear fits.
#'
#' @rdname growthrates-classes
#' @exportClass multiple_nonlinear_fits
#'
setClass("multiple_nonlinear_fits",
         contains = "multiple_fits"
)

#' \code{multiple_smooth.spline_fits}: class representing multiple smooth.spline fits.
#'
#' @rdname growthrates-classes
#' @exportClass multiple_smooth.spline_fits
#'
setClass("multiple_smooth.spline_fits",
         contains = "multiple_fits"
)

