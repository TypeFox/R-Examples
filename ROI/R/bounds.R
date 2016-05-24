################################################################################
## 'bounds'

################################################################################
## 'V_bound' constructor

##' Constructs a variable bounds object.
##'
##' This function returns a sparse representation of objective
##' variable bounds.
##' @title Objective Variable Bounds
##' @param li an integer vector specifying the indices of non-standard
##' (i.e., values != 0) lower bounds.
##' @param ui an integer vector specifying the indices of non-standard
##' (i.e., values != Inf) upper bounds.
##' @param lb a numeric vector with lower bounds.
##' @param ub a numeric vector with upper bounds.
##' @param nobj an integer representing the number of objective variables
##' @return An S3 object of class \code{"V_bound"} containing lower and
##' upper bounds of the objective variables.
##' @author Stefan Theussl
##' @export
V_bound <- function( li, ui, lb, ub, nobj = max(li, ui) ) {
    if(missing(li))
        li <- integer()
    if(missing(ui))
        ui <- integer()
    if(missing(lb))
        lb <- double()
    if(missing(ub))
        ub <- double()
    li <- as.integer(li)
    ui <- as.integer(ui)
    lb <- as.double(lb)
    ub <- as.double(ub)
    if(length(lb)){
        zero <- lb == 0
        lb <- lb[!zero]
        li <- li[!zero]
    }
        if(length(ub)){
        inf <- ub == Inf
        ub <- ub[!inf]
        ui <- ui[!inf]
    }

    ## Sanity checking
    if( (length(li) != length(lb)) || (length(ui) != length(ub)) )
        stop("length of indices must be equal to the length of the corresponding values.")
    if( any(duplicated(li)) || any(duplicated(ui)) )
        stop("duplicated entries in indices.")
    if(length(li))
        if( max(li) > nobj )
            stop("indices must not exceed number of objective coefficients.")
    if(length(ui))
        if( max(ui) > nobj )
            stop("indices must not exceed number of objective coefficients.")
    if( any(lb >= Inf) )
      stop("lower bound cannot be 'Inf'.")
    if( any(ub <= -Inf) )
        stop("upper bounds cannot be '-Inf'.")
    ## FIXME: lower bounds vs. upper bounds -> lb cannot be higher than ub and
    ##        the other way round
    structure( list(lower = list(ind = li, val = lb),
                    upper = list(ind = ui, val = ub),
                    nobj = as.integer(nobj)),
              class = "V_bound" )
}

as.V_bound <- function( x ){
    UseMethod( "as.V_bound" )
}

##' @noRd
##' @S3method as.V_bound V_bound
as.V_bound.V_bound <- identity

##' @noRd
##' @S3method as.V_bound NULL
as.V_bound.NULL <- function( x )
    .make_standard_bounds()

##' @noRd
##' @method as.list V_bound
##' @S3method as.list V_bound
as.list.V_bound <- function( x, ... )
  unclass( x )

##' @noRd
##' @method print V_bound
##' @S3method print V_bound
print.V_bound <- function(x, ...){
    writeLines( "ROI Variable Bounds:\n" )

    writeLines( sprintf("%d lower and %d upper non-standard variable bounds.",
                        length(x$lower$ind), length(x$upper$ind)) )
}


################################################################################
## 'bounds' extractor functions

##' Extract bounds from its argument (typically \pkg{ROI} objects) and
##' return them.
##'
##' Currently, there is no default method. See \code{\link{bounds.OP}}
##' for extracting bounds from \pkg{ROI} objects of class \code{"OP"}.
##' @title Extract Objective Variable Bounds
##' @param x an object used to select the method.
##' @return the extracted bounds object.
##' @author Stefan Theussl
##' @export
bounds <- function( x )
  UseMethod("bounds")

##' Extract bounds from ROI objects of class \code{"OP"} and return
##' them.
##'
##' @title Extract Objective Variable Bounds
##' @param x an object of class \code{"OP"}.
##' @return an object of class \code{"V_bound"}.
##' @author Stefan Theussl
##' @method bounds OP
##' @S3method bounds OP
bounds.OP <- function( x )
   x$bounds

################################################################################
## 'bounds' replacement functions

##' Replaces the (variable) bounds in R objects (typically ROI
##' objects of class \code{"OP"}).
##'
##' Currently, there is no default method. Bounds in ROI objects of
##' class \code{"OP"} given by the argument \code{x} are replaced with
##' \code{value}, either being an object of class \code{"V_bound"} or
##' \code{NULL} (standard variable bound). The updated \code{"OP"}
##' object will be returned.
##' @title Replacement of Variable Bounds
##' @name bounds-replace
##' @aliases bounds<- bounds<-.OP
##' @usage bounds(x) <- value
##' @param x an R object.
##' @param value an R object.
##' @return the updated object.
##' @author Stefan Theussl
##' @export bounds<-
'bounds<-' <- function( x, value )
  UseMethod("bounds<-")


##' @noRd
##' @S3method bounds<- OP
'bounds<-.OP' <- function( x, value ) {
   x$bounds <- as.V_bound(value)
   x
}

.make_standard_bounds <- function( x )
  NULL
