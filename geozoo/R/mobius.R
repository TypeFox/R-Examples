#' Mobius
#'
#' A function to generate a mobius strip in the third or fourth dimension.
#'
#' @param p dimension of object.  (3)
#' @param n number of points
#' @references \url{http://schloerke.github.io/geozoo/mobius/mobius/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a mobius strip.
#' mobius(3, n = 1000)
#'
#' @keywords dynamic
#' @export
mobius <- function(p = 3, n = 10000){

  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, mobius_row(3))
      )
    ),
    ncol = 3,
    byrow = TRUE
  )
  wires <- NULL

  structure(
    list( points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
mobius_row <- function(p = 3) {

  ##Generates Angles
  a <- runif(1, min = 0, max = 2 * pi)
  a <- c(a, a / 2)

  ##Generates Small Radius
  radius <- c(1, runif(1, min = -.4, max = .4))

  ##Generates Row of Data
  mobius <- c(
    (cos(a[2]) * radius[2] + radius[1]) * cos(a[1]),
    (cos(a[2]) * radius[2] + radius[1]) * sin(a[1]),
    sin(a[2]) * radius[2]
  )

  mobius
}

#' Mobius Experiment
#'
#' A function to generate a 5-D mobius strip in the third dimension.
#'
#' @param p dimension of object.  (5)
#' @param n number of points
#' @references \url{http://schloerke.github.io/geozoo/mobius/mobius/}
#' @author Barret Schloerke
#' @examples
#' ## Generates a mobius strip.
#' mobius.experiment(5, n = 1000)
#'
#' @keywords dynamic
#' @export
mobius.experiment <- function(p = 5, n = 10000){

  p <- 5

  vert <- matrix(
    do.call(
      "rbind",
      as.list(
        replicate(n, mobius_experiment_row())
      )
    ),
    ncol = 3,
    byrow = TRUE
  )
  wires <- NULL

  structure(
    list( points = vert, edges = wires),
    class = c("geozooNoScale", "geozoo")
  )
}

#' @keywords internal
mobius_experiment_row <- function(){

  ##Generates Angles
  a <- runif(1, min = 0, max = 2 * pi)
  a <- c(a, a / 2)

  ##Generates Small Radius
  radius <- c(1, runif(1, min = -.4, max = .4))

  ##Generates Row of Data
  mobius <- c(
    (cos(a[2]) * radius[2] + radius[1]) * cos(a[1]),
    (cos(a[2]) * radius[2] + radius[1]) * sin(a[1]),
    sin(a[2]) * radius[2]
  )

  k <- runif(1, min = 0, max = pi)
  ## Rot over x axis
  rot_1 <- matrix(
    c(
      0,
      cos(k),
      -sin(k),
      1, 0, 0, 0,
      sin(k),
      cos(k)
    ),
    ncol = 3, byrow = TRUE
  )
  ## Rot over z axis
  rot_2 <- matrix(
    c(
      cos(2 * k),
      -sin(2 * k),
      0,
      sin(2 * k),
      cos(2 * k),
      0, 0, 0, 1
    ),
    ncol = 3, byrow = TRUE
  )
  ## Trans perpendicular to z axis
  trans <- matrix(c(4 * cos(2 * k), 4 * sin(2 * k), 0), ncol = 1)

  mobius <- rot_2 %*% rot_1 %*% mobius + trans

  mobius
}
