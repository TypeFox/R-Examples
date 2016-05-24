# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Access or set information on cross-validation results
#' 
#' Retrieve or set the names of cross-validation results, retrieve or set the 
#' identifiers of the models, or retrieve the number of cross-validation 
#' results or included models.
#' 
#' @rdname accessors
#' @name accessors
#' 
#' @param x  an object inheriting from class \code{"cv"} or \code{"cvSelect"} 
#' that contains cross-validation results.
#' @param value  a vector of replacement values.
#' 
#' @return 
#' \code{cvNames} returns the names of the cross-validation results.  The 
#' replacement function thereby returns them invisibly.
#' 
#' \code{fits} returns the identifiers of the models for objects inheriting 
#' from class \code{"cvSelect"} and \code{NULL} for objects inheriting from 
#' class \code{"cv"}.  The replacement function thereby returns those values 
#' invisibly.
#' 
#' \code{ncv} returns the number of cross-validation results.
#' 
#' \code{nfits} returns the number of models included in objects inheriting 
#' from class \code{"cvSelect"} and \code{NULL} for objects inheriting from 
#' class \code{"cv"}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvSelect}}, \code{\link{cvTuning}}
#' 
#' @example inst/doc/examples/example-accessors.R
#' 
#' @keywords utilities

NULL


#' @rdname accessors
#' @export

cvNames <- function(x) UseMethod("cvNames")

#' @S3method cvNames cv
cvNames.cv <- function(x) names(x$cv)

#' @S3method cvNames cvSelect
cvNames.cvSelect <- function(x) names(x$cv)[-1]


#' @rdname accessors
#' @usage cvNames(x) <- value
#' @export

"cvNames<-" <- function(x, value) UseMethod("cvNames<-")

#' @S3method cvNames<- cv
"cvNames<-.cv" <- function(x, value) {
    object <- x
    names(object$cv) <- names(object$se) <- value
    if(!is.null(x$reps)) colnames(object$reps) <- value
    eval.parent(substitute(x <- object))
}

#' @S3method cvNames<- cvSelect
"cvNames<-.cvSelect" <- function(x, value) {
    object <- x
    names(object$best) <- value
    value <- c("Fit", value)
    names(object$cv) <- names(object$se) <- value
    if(!is.null(x$reps)) names(object$reps) <- value
    eval.parent(substitute(x <- object))
}


#' @rdname accessors
#' @export

fits <- function(x) UseMethod("fits")

#' @S3method fits cv
fits.cv <- function(x) NULL

#' @S3method fits cvSelect
fits.cvSelect <- function(x) x$cv$Fit


#' @rdname accessors
#' @usage fits(x) <- value
#' @export

"fits<-" <- function(x, value) UseMethod("fits<-")

#' @S3method fits<- cv
"fits<-.cv" <- function(x, value) eval.parent(substitute(x))

#' @S3method fits<- cvSelect
"fits<-.cvSelect" <- function(x, value) {
    object <- x
    if(is.factor(value)) value <- factor(as.character(value), levels=value)
    object$cv$Fit <- object$se$Fit <- value
    if(!is.null(reps <- x$reps)) {
        indices <- match(reps$Fit, x$cv$Fit, nomatch=0)
        object$reps$Fit <- value[indices]
    }
    eval.parent(substitute(x <- object))
}


#' @rdname accessors
#' @export

ncv <- function(x) UseMethod("ncv")

#' @S3method ncv cv
ncv.cv <- function(x) length(x$cv)

#' @S3method ncv cvSelect
ncv.cvSelect <- function(x) ncol(x$cv) - 1


#' @rdname accessors
#' @export

nfits <- function(x) UseMethod("nfits")

#' @S3method nfits cv
#' @S3method nfits cvSelect
nfits.cv <- nfits.cvSelect <- function(x) nrow(x$cv)
