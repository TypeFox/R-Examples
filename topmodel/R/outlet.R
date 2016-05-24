outlet <- function(map,out,radius = 2) {

  ## TODO: plot this as values with coordinates

  i <- out[1]
  j <- out[2]

  return(map[(i-radius):(i+radius),(j-radius):(j+radius)])

}
