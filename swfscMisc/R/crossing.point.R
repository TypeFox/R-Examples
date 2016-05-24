#' @title Crossing Point
#' @description Return point where two lines cross
#'
#' @param l1,l2 matrices representing two lines, where first two columns are
#'   x and y values respectively
#'
#' @return a data.frame of x and y values of points where lines cross
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' x <- 1:100
#' line1 <- cbind(x, 3 + 3 * x)
#' line2 <- cbind(x, 10 - 3 * x)
#' plot(line1[, 1], line1[, 2], type = "l", col = "red")
#' lines(line2[, 1], line2[, 2], col = "blue")
#' cr.pt <- crossing.point(line1, line2)
#' print(cr.pt)
#' 
#' @importFrom spatstat crossing.psp as.owin as.data.frame.ppp psp
#' @export
#' 
crossing.point <- function(l1, l2) {
  i1.1 <- 1:(nrow(l1) - 1)
  i1.2 <- 2:nrow(l1)
  i2.1 <- 1:(nrow(l2) - 1)
  i2.2 <- 2:nrow(l2)
  xlim <- range(c(l1[, 1], l2[, 1]))
  ylim <- range(c(l1[, 2], l2[, 2]))
  w <- as.owin(list(xrange = xlim, yrange = ylim))
  psp1 <- psp(l1[i1.1, 1], l1[i1.1, 2], l1[i1.2, 1], l1[i1.2, 2], window = w)
  psp2 <- psp(l2[i2.1, 1], l2[i2.1, 2], l2[i2.2, 1], l2[i2.2, 2], window = w)
  cross <- crossing.psp(psp1, psp2, fatal = FALSE)
  if(is.null(cross)) return(NULL) else as.data.frame(cross)
}
