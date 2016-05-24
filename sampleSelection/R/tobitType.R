### Type of the corresponding tobit model.
### It should be done by classes instead..

tobitType <- function( x, ... ) {
    UseMethod("tobitType")
}

tobitType.selection <- function( x, ... ) {
   x$tobitType
}
