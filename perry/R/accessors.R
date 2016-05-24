# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Access or set information on resampling-based prediction error results
#' 
#' Retrieve or set the names of resampling-based prediction error results, 
#' retrieve or set the identifiers of the models, or retrieve the number of 
#' prediction error results or included models.
#' 
#' @rdname accessors
#' @name accessors
#' 
#' @param x  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param value  a vector of replacement values.
#' 
#' @return 
#' \code{peNames} returns the names of the prediction error results.  The 
#' replacement function thereby returns them invisibly.
#' 
#' \code{fits} returns the identifiers of the models for objects inheriting 
#' from class \code{"perrySelect"} and \code{NULL} for objects inheriting from 
#' class \code{"perry"}.  The replacement function thereby returns those values 
#' invisibly.
#' 
#' \code{npe} returns the number of prediction error results.
#' 
#' \code{nfits} returns the number of models included in objects inheriting 
#' from class \code{"perrySelect"} and \code{NULL} for objects inheriting from 
#' class \code{"perry"}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perrySelect}}, 
#' \code{\link{perryTuning}}
#' 
#' @example inst/doc/examples/example-accessors.R
#' 
#' @keywords utilities

NULL


#' @rdname accessors
#' @export
peNames <- function(x) UseMethod("peNames")

#' @S3method peNames perry
peNames.perry <- function(x) names(x$pe)

#' @S3method peNames perrySelect
peNames.perrySelect <- function(x) names(x$pe)[-1]


#' @rdname accessors
#' @usage peNames(x) <- value
#' @export
"peNames<-" <- function(x, value) UseMethod("peNames<-")

#' @S3method peNames<- perry
"peNames<-.perry" <- function(x, value) {
  object <- x
  names(object$pe) <- names(object$se) <- value
  if(hasComponent(object, "reps")) colnames(object$reps) <- value
  eval.parent(substitute(x <- object))
}

#' @S3method peNames<- perrySelect
"peNames<-.perrySelect" <- function(x, value) {
  object <- x
  names(object$best) <- value
  value <- c("Fit", value)
  names(object$pe) <- names(object$se) <- value
  if(hasComponent(object, "reps")) names(object$reps) <- value
  eval.parent(substitute(x <- object))
}


#' @rdname accessors
#' @export
fits <- function(x) UseMethod("fits")

#' @S3method fits perry
fits.perry <- function(x) NULL

#' @S3method fits perrySelect
fits.perrySelect <- function(x) x$pe$Fit


#' @rdname accessors
#' @usage fits(x) <- value
#' @export
"fits<-" <- function(x, value) UseMethod("fits<-")

#' @S3method fits<- perry
"fits<-.perry" <- function(x, value) eval.parent(substitute(x))

#' @S3method fits<- perrySelect
"fits<-.perrySelect" <- function(x, value) {
  object <- x
  if(is.factor(value)) value <- factor(as.character(value), levels=value)
  object$pe$Fit <- object$se$Fit <- value
  if(!is.null(reps <- x$reps)) {
    indices <- match(reps$Fit, x$pe$Fit, nomatch=0)
    object$reps$Fit <- value[indices]
  }
  names(object$yHat) <- value
  eval.parent(substitute(x <- object))
}


#' @rdname accessors
#' @export
npe <- function(x) UseMethod("npe")

#' @S3method npe perry
npe.perry <- function(x) length(x$pe)

#' @S3method npe perrySelect
npe.perrySelect <- function(x) ncol(x$pe) - 1


#' @rdname accessors
#' @export
nfits <- function(x) UseMethod("nfits")

#' @S3method nfits perry
#' @S3method nfits perrySelect
nfits.perry <- nfits.perrySelect <- function(x) nrow(x$pe)


# ## @rdname accessors
# ## @export
# getR <- function(x) UseMethod("getR")
# 
# ## @S3method getR cvFolds
# ## @S3method getR randomSplits
# ## @S3method getR bootSamples
# getR.cvFolds <- getR.randomSplits <- getR.bootSamples <- function(x) x$R
# 
# ## @S3method getR perry
# ## @S3method getR perrySelect
# getR.perry <- getR.perrySelect <- function(x) getR(x$splits)

# -------------------------

#' @S3method coef perryTuning
coef.perryTuning <- function(object, ...) {
  finalModel <- object$finalModel
  if(is.null(finalModel)) stop("final model not available")
  coef(finalModel, ...)
}

#' @S3method fitted perryTuning
fitted.perryTuning <- function(object, ...) {
  finalModel <- object$finalModel
  if(is.null(finalModel)) stop("final model not available")
  fitted(finalModel, ...)
}

#' @S3method residuals perryTuning
residuals.perryTuning <- function(object, ...) {
  finalModel <- object$finalModel
  if(is.null(finalModel)) stop("final model not available")
  residuals(finalModel, ...)
}
