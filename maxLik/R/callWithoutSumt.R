## strip possible SUMT parameters and call the function thereafter
callWithoutSumt <- function(theta, fName, ...) {
   return( callWithoutArgs( theta, fName = fName, 
      args = names(formals(sumt)), ... ) )
}
