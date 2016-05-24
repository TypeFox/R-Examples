
xport.NA <- function()
  {
    .C("fill_numeric_NA", PACKAGE="SASxport")
    .Call("getRawBuffer", PACKAGE="SASxport")
  }
