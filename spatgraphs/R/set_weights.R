#' Set weights to edges of sg
#'
#' For each edge e(i,j) between points i,j, set the weight f(||x_i-x_j||)
#'
#' @param g sg object
#' @param x point pattern used in g
#' @param f function for the weight
#' @param scale additional scale parameter for the default f
#' @param ... ignored
#'
#' @details
#' Default f(x) = exp(-x^2/scale)
#' @export

weight.sg<-function(g, x, f=function(x)exp(-x^2/scale), scale=1, ...) {
  if(!is(g, "sg")) stop("g not an sg object.")
  D<-as.matrix(dist(sg_parse_coordinates(x), upper=T, diag=T))
  W<-f(D)
  weights<-list()
  for(i in 1:g$N)	weights[[i]]<-W[i, g$edges[[i]]]
  g$weights<-weights
  g$note <- c(g$note, "weighted edges")
  g
}
