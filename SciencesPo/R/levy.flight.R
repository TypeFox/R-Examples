#' @title Simulates a Levy Walk
#'
#' @description This function simulates a Levy walk
#' @param n The lenght of walk.
#' @param alpha The exponent of the Levy distribution.
#' @param min.lenght The minimum length of a step.
#' @param max.lenght The maximum length of a step.
#' @param plot A logical, TRUE will make a plot.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' #' @examples
#' levy.flight(n=100, alpha=2)
#' @export
`levy.flight`<- function (n=500, alpha = 3, min.lenght=1, max.lenght=5, plot=TRUE){
  theta <- runif(n-1)*2*pi
  f <- runif(n-1, max.lenght^(-alpha), min.lenght^(-alpha))^(-1/alpha)
  x <- c(0, cumsum(f*cos(theta)))
  y <- c(0, cumsum(f*sin(theta)))
  out<- data.frame(x,y)
  if(plot){
    plot(x~y, type="l")
  }
  return(out)
}
NULL
