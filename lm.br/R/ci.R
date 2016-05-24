
# this function is internal, not meant for the user

# R code to process user input, then call
# the corresponding function of the C++ object


.ci <- function( z, CL =0.95, method ="clr" )  {
  method <- toupper(method)
  met <- integer(1)
  if( method=="CLR" )  met <- 1  else  {
    if( method=="AF" )  met <- 2  else
      stop( "'method' must be \"CLR\" or \"AF\"" )
  }
  (z$CppObj)$ci( CL, met )
}
