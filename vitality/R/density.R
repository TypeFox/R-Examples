####################################################################################################
## Density functions

#' Density function for 3-param r, s, u
#' 
#' None
#' @param xx age
#' @param r r value
#' @param s s value
#' @param u u value
#' @return density
ft.4p <- function(xx, r, s, u) {
  temp1 <- s^2 * xx + u^2
  temp2 <- u^2 * r + s^2
  if (xx==0) value = 0 
  else value <- temp2 / sqrt(2 * pi * temp1^3) * exp(-(1-r * xx)^2/(2 * temp1))
  return(value)
}


#' Vectorized density function
#' 
#' None
#' 
#' @param xx vector of ages
#' @param r r value
#' @param s s value
#' @param u u value
#' @return vector of densities
vft.4p <- Vectorize(FUN = ft.4p, vectorize.args = "xx")

####################################################################################################
## Density functions for 6 par

#' Density function for intrinsic
#' 
#' None
#' @param xx age
#' @param r r value
#' @param s s value
#' @return density
ft.6p <- function(xx, r, s) {
  temp1 <- s^2 * xx
  temp2 <- s^2
  if (xx==0) value = 0 
  else value <- ((xx^-(3/2)) / (s*sqrt(2 * pi))) * exp(-(1-r * xx)^2/(2 * temp1))
  return(value)
}


#' Vectorized density function
#' 
#' None
#' 
#' @param xx vector of ages
#' @param r r value
#' @param s s value
#' @return vector of densities
vft.6p <- Vectorize(FUN = ft.6p, vectorize.args = "xx")
