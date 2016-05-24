#' f_helmert
#'
#' Function to set up a Helmert transformation of a (d-1)-dimensional 
#' shape in p-space down into its (p-1)-space. Mostly internally used, 
#' but could be useful for setting up new compositional data.
#'
#' @param d object
#' @return
#'  \item{helmert }{transformation matrix}
#' @references \url{http://schloerke.github.io/geozoo/simplices/}
#' @author Di Cook
#' @export
f_helmert <- function(d){
  helmert <- rep(1 / sqrt(d), d)
  for (i in 1:(d - 1)){
    x <- rep(1 / sqrt(i * (i + 1)), i)
    x <- c(x, -i / sqrt(i * (i + 1)))
    x <- c(x, rep(0, d - i - 1))
    helmert <- rbind(helmert, x)
  }
  return(helmert)
}

#' f_composition
#'
#' Function to take a d-dimensional compositional data set and 
#' transform it using a Helmert transformation into (p-1)-space, 
#' where it lives. Mostly internally used, but could be useful 
#' for setting up new compositional data.
#'
#' @param data object
#' @return
#'  \item{data }{points in (d-1)-dimensional space}
#' @references \url{http://schloerke.github.io/geozoo/simplices/}
#' @author Di Cook
#' @export
f_composition <- function(data){
  d <- dim(data)[2]
  hm <- f_helmert(d)
  x <- data - matrix(1 / d, dim(data)[1], d)
  return( (x %*% t(hm))[, -1])
}


#' Simplex
#'
#' A function to generate a simplex
#'
#' @param p dimension of object
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/simplices/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a simplex
#' simplex(p = 3)
#'
#' @keywords dynamic
#' @export
simplex <- function(p = 3){
  vert <- f_composition(diag(p + 1))

  wires <- do.call(expand.grid, list(c(1:nrow(vert)), c(1:nrow(vert))))

  structure(
    list(
      points = vert,
      edges = wires[! (wires[, 1] == wires[, 2]), ]
    ),
    class = c("geozooNoScale", "geozoo")
  )
}
