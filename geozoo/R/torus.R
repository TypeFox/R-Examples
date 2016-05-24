
#' Torus
#'
#' A function to generate a torus in any dimension
#'
#' @param p dimension of object
#' @param n number of points
#' @param radius radiuses of the torus, set from largest to smallest
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/torus/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a torus
#' torus(p = 3, n = 1000)
#'
#' @keywords dynamic
#' @export
torus <- function(p = 3, n = 10000, radius = 2 ^ ( (p - 2):0) ){
  if (length(radius) != (p - 1)) {
    stop("'radius' length does not equal p")
  }

  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, torus_row(radius, p))
      )
    ),
    ncol = p, byrow = TRUE
  )
  wires <- NULL

  structure(
    list(points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
torus_row <- function(radius, p) {
  ##Generates Angles
  pm1 <- p - 1
  t <- runif(pm1, min = 0, max = 2 * pi)

  ##Generates Row of Data
  torus <- c(
    rep(cos(t[pm1]) * radius[pm1], pm1),
    sin(t[pm1]) * radius[pm1]
  )

  if (p > 2) {
    for (i in (pm1):2) {
      for (j in (i - 1):1) {
        torus[j] <- (torus[j] + radius[i - 1]) * cos(t[i - 1])
      }
      torus[i] <- (torus[i] + radius[i - 1]) * sin(t[i - 1])
    }
  }
  torus
}


#' Flat Torus
#'
#' A function to generate a flat torus in any dimension
#'
#' @param p dimension of object (number of circles x2)
#' @param n number of points
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/torus/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a Flat Torus
#' torus.flat(p = 4, n = 1000)
#'
#' @keywords dynamic
#' @export
torus.flat <- function(p = 4, n = 10000){
  p <- floor(p / 2)
  vert <- do.call("rbind", replicate(n, torus_flat_row(p), simplify = FALSE))
  wires <- NULL
  structure(
    list(points = vert, edges = wires),
    class = "geozoo"
  )
}

#' @keywords internal
torus_flat_row <- function(p){
  a <- runif(p, min = 0, max = 2 * pi)
  as.vector(rbind(cos(a), sin(a)))
}
