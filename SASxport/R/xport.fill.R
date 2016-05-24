xport.fill <- function( useSpace, width )
  {
    
    
    .C("fill_space",
       value = as.integer(useSpace),
       width = as.integer(width),
       PACKAGE="SASxport"
       )

    .Call("getRawBuffer", PACKAGE="SASxport")
  }
