
## Solver functions for a given OP signature

## Signature: objective_constraints_types<_bounds_maximize>

## Signature: L_L_C<_FALSE_FALSE> (LP)
L_L_C <- function(x)
    UseMethod( "L_L_C" )

## Signature: L_L_I<_FALSE_FALSE> (LIP)
L_L_I <- function(x)
    UseMethod( "L_L_I" )

## Signature: L_L_C+I<_FALSE_FALSE> (MILP)
L_L_CI <- function(x)
    UseMethod( "L_L_CI" )

all_signatures <- function(){
    sigs <- expand.grid( objective = names(available_objective_classes()),
                         constraints = names(available_objective_classes()),
                         types = c("C", "I", "B", "CI", "CB", "IB", "CIB"),
                         bounds = c("TRUE", "FALSE"),
                         maximum = c("TRUE", "FALSE") )
    stopifnot( ncol(sigs) == length(formals(OP)) )
    stopifnot( identical(colnames(sigs), names(formals(OP))) )
    sigs
}

#make_generics <- function(){
#    lapply( apply( all_signatures(), 1, .make_signature), function(x) sprintf('%s <- function(x) UseMethod( "%s" )', x, x))
#}
