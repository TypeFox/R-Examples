
##==============================================================================
## shadepalette    : creates a palette that is suited for shading
##==============================================================================

shadepalette <- function(n=100, endcol = "red", inicol = "white",
  interval = c(0.0,1.0)) {

  x.to <- seq ( interval[1], interval[2], length.out=n)
  return( intpalette( rbind(inicol,endcol), n, x.to=x.to))

}
