##Parametric Eqn's

#' Boy Surface
#'
#' A function to produce a Boy Surface.
#'
#' @param n number of points
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a Boy Surface
#' boy.surface(n = 1000)
#'
#' @keywords dynamic
#' @export
boy.surface <- function(n = 10000) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, boy_surface_row())
      )
    ),
    ncol = 3, byrow = TRUE
  )
  wires <- NULL

  structure(
    list(points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
boy_surface_row <- function( ){
  u <- runif( 1, min = 0, max = pi)
  v <- runif( 1, min = 0, max = pi)
  a <- cos( v) * sin( u)
  b <- sin( v) * sin( u)
  c <- cos( u)
  x <- 1 / 2 * (
      ( 2 * a ^ 2 - b ^ 2 - c ^ 2) +
      2 * b * c * ( b ^ 2 - c ^ 2) +
      c * a * ( a ^ 2 - c ^ 2) +
      a * b * ( b ^ 2 - a ^ 2)
    )
  y <- sqrt(3) / 2 * (
    ( b ^ 2 - c ^ 2) +
    c * a * ( c ^ 2 - a ^ 2) +
    a * b * ( b ^ 2 - a ^ 2)
  )
  z <- ( a + b + c) * (
    ( a + b + c) ^ 3 +
    4 * ( b - a) * ( c - b) * ( a - c)
  )
  z <- z / 8
  return( cbind( x, y, z))
}


#' Conic Spiral
#'
#' A function to produce a conic spiral
#'
#' @param n number of points
#' @param a final radius of cone
#' @param b height of object
#' @param c inner radius
#' @param w number of spirals
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a Conic Spiral
#' conic.spiral(n = 1000)
#'
#' @keywords dynamic
#' @export
conic.spiral <- function(n = 10000, a = .2, b = 1, c = .1, w = 2) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, conic_spiral_row(a, b, c, w))
      )
    ),
    ncol = 3, byrow = TRUE
  )
  wires <- NULL

  structure(
    list(points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
conic_spiral_row <- function(a, b, c, w) {
  u <- runif( 1, min = 0, max = 2 * pi)
  v <- runif( 1, min = 0, max = 2 * pi)

  x <- a * ( 1 - v / ( 2 * pi)) * cos( w * v) * ( 1 + cos( u)) +
    c * cos( w * v)
  y <- a * ( 1 - v / ( 2 * pi)) * sin( w * v) * ( 1 + cos( u)) +
    c * sin( w * v)
  z <- ( b * v + a * ( 1 - v / ( 2 * pi)) * sin( u)) / ( 2 * pi)
  return( cbind( x, y, z))
}



#' Conic Spiral (Nautilus Shape)
#'
#' A function to produce a Conic Spiral in a nautilus shape
#'
#' @param n number of points
#' @param a final radius of cone
#' @param b height of object
#' @param c inner radius
#' @param w number of spirals
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a Nautilus Conic Spiral
#' conic.spiral.nautilus( n = 1000 )
#'
#' @keywords dynamic
#' @export
conic.spiral.nautilus <- function(n = 10000, a = .2, b = .1, c = 0, w = 2) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, conic_spiral_row(a, b, c, w))
      )
    ),
    ncol = 3, byrow = TRUE
  )
  wires <- NULL

  structure(
    list(points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}


#' Cross Cap
#'
#' A function to generate a cross cap
#'
#' @param n number of points
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a Cross Cap
#' cross.cap( n = 1000 )
#'
#' @keywords dynamic
#' @export
cross.cap <- function(n = 10000) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, cross_cap_row())
      )
    ),
    ncol = 3, byrow = TRUE
  )
  wires <- NULL

  structure(
    list(points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
cross_cap_row <- function( ) {
  u <- runif( 1, min = 0, max = pi)
  v <- runif( 1, min = 0, max = pi)
  x <- cos( u) * sin( 2 * v)
  y <- sin( u) * sin( 2 * v)
  z <- cos( v) * cos( v) - cos( u) * cos( u) * sin( v) * sin( v)
  return( cbind( x, y, z))
}


#' Dini Surface
#'
#' A function to generate a dini surface.
#'
#' @param n number of points
#' @param a outer radius of object
#' @param b space between loops
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a Dini Surface
#' dini.surface(n = 1000, a = 1, b = 1)
#'
#' @keywords dynamic
#' @export
dini.surface <- function(n = 10000, a = 1, b = 1) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, dini_surface_row( a, b))
      )
    ),
    ncol = 3, byrow = TRUE
  )
  wires <- NULL

  structure(
    list( points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
dini_surface_row <- function(a = 1, b = 1) {
  u <- runif( 1, min = 0, max = 4 * pi)
  v <- runif( 1, min = 0.0000000001, max = 2)
  x <- a * cos( u) * sin( v)
  y <- a * sin( u) * sin( v)
  z <- a * ( cos(v) + log(tan(v / 2))) + (b * u)
  return( cbind( x, y, z))
}


#' Ellipsoid
#'
#' A function to generate an ellipsoid
#'
#' @param n number of points
#' @param a radius in x direction
#' @param b radius in y direction
#' @param c radius in z direction
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates an ellipsoid
#' ellipsoid(n = 1000, a = 1, b = 1, c = 3)
#'
#' @keywords dynamic
#' @export
ellipsoid <- function(n = 10000, a = 1, b = 1, c = 3) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, ellipsoid_row(a, b, c))
      )
    ),
    ncol = 3, byrow = TRUE
  )
  wires <- NULL

  structure(
    list(points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
ellipsoid_row <- function(a, b, c) {
  u <- runif( 1, min = 0, max = 2 * pi)
  v <- runif( 1, min = 0, max = 2 * pi)
  x <- a * cos( u) * sin( v)
  y <- b * sin( u) * sin( v)
  z <- c * cos( v)
  return( cbind( x, y, z))
}


#' Enneper's Surface
#'
#' A function to generate Enneper's surface
#'
#' @param n number of points
#' @param a angle, radians, minimum and maximum.  -a < angle < a
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates an Enneper Surface
#' enneper.surface(n = 1000, a = 4)
#'
#' @keywords dynamic
#' @export
enneper.surface <- function(n = 10000, a = 4) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, enneper_surface_row(a))
      )
    ),
    ncol = 3, byrow = TRUE
  )
  wires <- NULL

  structure(
    list( points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
enneper_surface_row <- function( a = 4) {
  u <- runif( 1, min = - a, max = a)
  v <- runif( 1, min = - a, max = a)
  x <- u - u ^ 3 / 3 + u * v ^ 2
  y <- v - v ^ 3 / 3 + v * u ^ 2
  z <- u ^ 2 - v ^ 2
  return( cbind( x, y, z))
}



#' Figure Eight Klein Bottle
#'
#' A function to generate a figure eight Klein bottle
#'
#' @param n number of points
#' @param a radius of outer radius
#' @param b radius of inner radius
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a figure eight Klein bottle.
#' klein.fig.eight(n = 1000, a = 3, b = 1)
#'
#' @keywords dynamic
#' @export
klein.fig.eight <- function(n = 10000, a = 3, b = 1) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, klein_fig_eight_row(a, b))
      )
    ),
    ncol = 4,
    byrow = TRUE
  )
  wires <- NULL

  structure(
    list( points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
klein_fig_eight_row <- function( a = 3, b = 1) {
  u <- runif( 1, min = - pi, max = pi)
  v <- runif( 1, min = - pi, max = pi)
  x <- ( b * cos( v) + a) * cos( u)
  y <- ( b * cos( v) + a) * sin( u)
  z <- b * sin( v) * cos( u / 2)
  w <- b * sin( v) * sin( u / 2)
  return( cbind( x, y, z, w))
}



#' Roman Surface
#'
#' A function to generate a Roman surface, also known as a Steiner surface
#'
#' @param n number of points
#' @param a maximum radius of object
#' @return
#'  \item{points }{location of points}
#'  \item{edges }{edges of the object (null)}
#' @references \url{http://schloerke.github.io/geozoo/mobius/other/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a Roman surface.
#' roman.surface(n = 1000, a = 1)
#'
#' @keywords dynamic
#' @export
roman.surface <- function(n = 10000, a = 1) {
  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, roman_surface_row(a))
      )
    ),
    ncol = 3, byrow = TRUE
  )
  wires <- NULL

  structure(
    list(points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
roman_surface_row <- function(a = 1) {
  u <- runif( 1, min = 0, max = pi)
  v <- runif( 1, min = 0, max = pi)
  x <- a ^ 2 * cos( v) * cos( v) * sin( 2 * u) / 2
  y <- a ^ 2 * sin( u) * sin( 2 * v) / 2
  z <- a ^ 2 * cos( u) * sin( 2 * v) / 2
  return( cbind( x, y, z))
}
