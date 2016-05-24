#' Get the apex for a gamma detection function
#'
#' @param ddfobj ddf object
#'
#' @return the distance at which the gamma peaks
#' @author Jeff Laake
apex.gamma <- function(ddfobj){
  key.scale <- scalevalue(ddfobj$scale$parameters,ddfobj$scale$dm)
  key.shape <- scalevalue(ddfobj$shape$parameters,ddfobj$shape$dm)

  fr <- (1/gamma(key.shape)) * (((key.shape - 1)/exp(1))^(key.shape - 1))
  return(key.scale*fr*(key.shape-1))
}
