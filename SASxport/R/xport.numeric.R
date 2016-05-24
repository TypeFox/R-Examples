xport.numeric <- function( value )
  {
    if(length(value)!=1) stop("Only a single numeric value is permitted.")
    if(is.na(value))
      return( xport.NA() )

    .C("fill_numeric_field",
       value = as.double(value),
       NAOK=TRUE,
       PACKAGE="SASxport"
       )

    .Call("getRawBuffer", PACKAGE="SASxport")
  }

