linearPredictors <- function( x, ... ) {
    UseMethod("linearPredictors")
}

linearPredictors.probit <- function( x, ... ) {
   mm <- naresid(na.action(x), model.matrix(x))
                           # naresid works (for now).  Should we keep it or replace it?
   result <- drop( mm %*% x$estimate )
   names( result ) <- rownames( mm )
   return( result )
}
