### treatReg: which is the name of the selection outcome variable
### used for determining if it is included in the outcome model

selectionVariableName <- function( x, ... ) {
    UseMethod("selectionVariableName")
}

selectionVariableName.selection <- function( x, ... ) {
   if(tobitType(x) == "treatment")
      return(x$param$selectionVariableName)
   return(integer(0))
}
