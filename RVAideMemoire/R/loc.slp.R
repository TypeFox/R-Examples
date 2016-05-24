loc.slp <- function() {
  coord <- locator(2)
  slp <- (coord$y[2]-coord$y[1])/(coord$x[2]-coord$x[1])
  return(slp)
}
