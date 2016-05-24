#' Check if Points as Inside a Polygon
#'
#' @param pontos A vector
#' @param limites A vector
#' @return Vector of points
#' @export
setDentroFora <- function(pontos, limites) {
  pontos@data$DENTRO<-!is.na(over(pontos, as(limites, "SpatialPolygons")))
  plot(limites)
  plot(pontos[pontos@data$DENTRO==T,],add=T, col="blue")
  plot(pontos[pontos@data$DENTRO==F,],add=T, col="red")
  return(pontos)
}
