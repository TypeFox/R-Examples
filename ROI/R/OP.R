################################################################################
## Package: ROI
## File:    OP.R
## Author:  Stefan Theussl
## Changed: 2013-11-25
################################################################################



##' Optimization problem constructor
##'
##' @title Optimization Problem Constructor
##' @param objective an object inheriting from class \code{"objective"}.
##' @param constraints an object inheriting from class \code{"constraints"}.
##' @param bounds \code{NULL} (default) or a list with elements
##' \code{upper} and \code{lower} containing the indices and
##' corresponding bounds of the objective variables.  The default for
##' each variable is a bound between 0 and \code{Inf}.
##' @param types a character vector giving the types of the objective
##' variables, with \code{"C"}, \code{"I"}, and \code{"B"}
##' corresponding to continuous, integer, and binary, respectively, or
##' \code{NULL} (default), taken as all-continuous.  Recycled as
##' needed.
##' @param maximum a logical giving the direction of the
##' optimization. \code{TRUE} means that the objective is to maximize
##' the objective function, \code{FALSE} (default) means to minimize
##' it.
##' @return A list containing the optimal solution, with the following
##' components.
##' \item{solution}{the vector of optimal coefficients}
##' \item{objval}{the value of the objective function at the optimum}
##' \item{status}{an integer with status information about the
##' solution returned: 0 if the optimal solution was found, a non-zero
##' value otherwise}
##' \item{msg}{the status code and additional
##' information about the solution provided by the solver.}
##' @author Stefan Theussl
##' @export
OP <- function( objective, constraints = NULL, types = NULL, bounds = NULL,
  maximum = FALSE )
    .check_OP_for_sanity( structure(list(objective = as.objective(objective),
                                         constraints = as.constraint(constraints),
                                         bounds = bounds,
                                         types = as.types(types),
                                         maximum = as.logical(maximum)), class = "OP")
                         )

.check_OP_for_sanity <- function( x ){
    if( length(objective(x)) != dim(constraints(x))[2] )
        stop( "dimensions of 'objective' and 'constraints' not conformable." )
    len_types <- length(types)
    if( len_types && (len_types > 1L) )
        if( length(objective(x)) != len_types )
            stop( "dimensions of 'objective' and 'types' not conformable." )
    if( !is.null(bounds(x)) )
        if( length(objective(x)) != bounds(x)$nobj )
            stop( "dimensions of 'objective' and 'bounds' not conformable." )
    x
}

## FIXME: also consider objective function

##' @noRd
##' @method print OP
##' @S3method print OP
print.OP <- function(x, ...){
    writeLines( "ROI Optimization Problem:\n" )
    ## objective
    op_type <- switch( class(objective(x))[2],
                       "L_objective" = "linear",
                       "Q_objective" = "quadratic",
                       "F_objective" = "nonlinear",
                       "abstract" )
    writeLines( sprintf("%s a %s objective function with",
                        ifelse(x$maximum, "Maximize", "Minimize"), op_type) )
    writeLines( sprintf("- %d objective variables,", length(objective(x))) )
    writeLines( "\nsubject to" )
    ## constraints
    types <- c( L_constraint = "linear",
                Q_constraint = "quadratic",
                F_constraint = "nonlinear" )
    writeLines( sprintf("- %d constraints of type %s.",
                        length(constraints(x)),
                        paste(na.omit(types[class(constraints(x))])[1],
                              collapse = ", ")
                        ) )
        ## calculate if types have to be written to stdout
    writetypes <- FALSE
    if( !is.null(types(x)) )
        if( any(types(x) %in% available_types()[2:3]) )
            writetypes <- TRUE
    if( writetypes ){
        writeLines( "" )
        writeLines( "Some of the objective variables are of type binary or integer." )
    }
}

##' Coerces objects of type \code{"OP"}.
##'
##' Objects from the following classes can be coerced to \code{"OP"}:
##' \code{"NULL"}, and \code{"numeric"}. The former represents an
##' empty optimization problem, the latter an unconstrained linear
##' programming problem where the elements of a \code{"numeric"}
##' vector \eqn{c} are treated as being objective variable
##' coefficients in \eqn{c^\top x}).  inherits from class
##' \code{"objective"}.
##' @title Optimization Problem Object
##' @param x an R object.
##' @return an object of class \code{"OP"}.
##' @author Stefan Theussl
##' @export
as.OP <- function(x)
    UseMethod("as.OP")

##' @noRd
##' @method as.OP OP
##' @S3method as.OP OP
as.OP.OP <- identity

##' @noRd
##' @method as.OP numeric
##' @S3method as.OP numeric
as.OP.numeric <- function(x){
    OP( objective = x, constraints = NULL, bounds = NULL, types = NULL,
        maximum = FALSE )

##' @noRd
##' @method as.OP default
##' @S3method as.OP default
as.OP.default <- function(x, ...)
    stop("Method not implemented.")

}

## OP_class <- function( x ){
##     x <- as.OP( x )
##     uniq_types <- if( is.null(types(x)) )
##         available_types()[1]
##     else unique(types(x))
##     signature <- list(

##                       )
##     c(sapply(available_objective_classes(),
##                           function(what) inherits(objective(x), what) ),
##                    sapply(available_constraint_classes(),
##                           function(what) inherits(constraints(x), what)),
##                    sapply(available_types(),
##                           function(what) what %in% uniq_types),
##                    bounds  = !is.null(bounds(x)),
##                    maximum = x$maximum
##     )
##     signature
## }

## NOTE: objective(x) returns something which inherits from function and class(x).
##       this is why we need to derive the type of objective by taking the 2nd element.
OP_signature <- function( x ){
    x <- as.OP( x )
    uniq_types <- if( is.null(types(x)) )
        available_types()[1]
    else paste(unique(types(x)), collapse = "")
    ROI_make_signature( objective = names( available_objective_classes() )[ available_objective_classes() %in% class(objective(x))[2] ],
                        constraints = names( available_constraint_classes() )[ available_constraint_classes() %in% class(constraints(x))[1] ],
                        types = uniq_types,
                        bounds  = !is.null(bounds(x)),
                        maximum = x$maximum
                       )
}

