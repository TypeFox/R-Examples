#' na.interp1
#'
#' This function combines pracma's \code{\link[pracma]{interp1}}
#' constant interpolation method with zoo's \code{\link[zoo]{na.approx}} linear
#' interpolation method. Here, \code{x = x} rather than
#' \code{x = index(object)} in na.approx. Here, \code{y = y} rather than
#' \code{y = object} in na.approx. Also, here, \code{xi} is used instead
#' of \code{xout} in na.approx. The Arguments list was obtained from both
#' interp1 and na.approx.
#'
#' @param x Numeric vector; points on the x-axis; at least two points
#' 	required; will be sorted if necessary.
#' @param y Numeric vector; values of the assumed underlying function;
#' 	\code{x} and \code{y} must be of the same length.
#' @param xi Numeric vector; points at which to compute the
#' 	interpolation; all points must lie between \code{min(x)} and
#' 	\code{max(x)}.
#' @param na.rm logical. If the result of the (\code{spline})
#' 	interpolation still results in \code{NA}s, should these be removed?
#' @param maxgap maximum number of consecutive \code{NA}s to fill. Any
#' 	longer gaps will be left unchanged. Note that all methods listed
#' 	above can accept \code{maxgap} as it is ultimately passed to the
#' 	default method.
#' @param ...	further arguments passed to methods. The \code{n}
#' 	argument of \code{approx} is currently not supported.
#'
#' @return Numeric vector representing values at points \code{xi}.
#'
#'
#' @author Hans Werner Borchers (pracma interp1), Felix Andrews
#'         (zoo na.approx), Irucka Embry
#'
#'
#' @source
#' \enumerate{
#'	\item zoo's na.approx.R - modified on Fri Aug 6 00:26:22 2010 UTC by felix. See \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/zoo/R/na.approx.R?view=markup&revision=781&root=zoo}.
#'	\item pracma interp1 function definition - R package pracma created and maintained by Hans Werner Borchers. See \code{\link[pracma]{interp1}}.
#' }
#'
#' @encoding UTF-8
#'
#'
#'
#' @seealso \code{\link[zoo]{na.approx}}, \code{\link[pracma]{interp1}}
#'
#'
#'
#'
#' @examples
#' library(iemisc)
#'
#' # zoo time series example
#' zoo1 <- structure(c(1.6, 1.7, 1.7, 1.7, 1.7, 1.7, 1.6, 1.7, 1.7, 1.7,
#' 1.7, 1.7, 2, 2.1, 2.1, NA, NA, 2.1, 2.1, NA, 2.3, NA, 2, 2.1), .Dim = c(12L,
#' 2L), .Dimnames = list(NULL, c("V1", "V2")), index = structure(c(1395242100,
#' 1395243000, 1395243900, 1395244800, 1395245700, 1395256500, 1395257400,
#' 1395258300, 1395259200, 1395260100, 1395261000, 1395261900), class =
#' c("POSIXct", "POSIXt"), tzone = "GMT"), class = "zoo")
#'
#' zoo1 <- as.data.frame(zoo1) # to data.frame from zoo
#' zoo1[, "Time"] <- as.POSIXct(rownames(zoo1)) # create column named Time as a
#' # POSIXct class
#' zoo1 <- setDT(zoo1) # create data.table out of data.frame
#' setcolorder(zoo1, c(3, 1, 2)) # set the column order as the 3rd column
#' # followed by the 2nd and 1st columns
#' zoo1 <- setDF(zoo1) # return to data.frame
#'
#' rowsinterps1 <- which(is.na(zoo1$V2 == TRUE))
#' # index of rows of zoo1 that have NA (to be interpolated)
#' xi <- as.numeric(zoo1[which(is.na(zoo1$V2 == TRUE)), 1])
#' # the Date-Times for V2 to be interpolated in numeric format
#' interps1 <- na.interp1(as.numeric(zoo1$Time), zoo1$V2, xi = xi,
#' na.rm = FALSE, maxgap = 1)
#' # the interpolated values where only gap sizes of 1 are filled
#' zoo1[rowsinterps1, 3] <- interps1
#' # replace the NAs in V2 with the interpolated V2 values
#' zoo1
#'
#'
#'
#' # data frame time series example
#' df1 <- structure(list(Time = structure(c(1395242100, 1395243000, 1395243900,
#'  1395244800, 1395245700, 1395256500, 1395257400, 1395258300, 1395259200,
#'  1395260100, 1395261000, 1395261900), class = c("POSIXct", "POSIXt"),
#'  tzone = "GMT"), V1 = c(1.6, 1.7, 1.7, 1.7, 1.7, 1.7, 1.6, 1.7, 1.7, 1.7,
#'  1.7, 1.7), V2 = c(2, 2.1, 2.1, NA, NA, 2.1, 2.1, NA, 2.3, NA, 2, 2.1)),
#'  .Names = c("Time", "V1", "V2"), row.names = c(NA, -12L),
#'  class = "data.frame")
#'
#' rowsinterps1 <- which(is.na(df1$V2 == TRUE))
#' # index of rows of df1 that have NA (to be interpolated)
#' xi <- as.numeric(df1[which(is.na(df1$V2 == TRUE)), 1])
#' # the Date-Times for V2 to be interpolated in numeric format
#' interps1 <- na.interp1(as.numeric(df1$Time), df1$V2, xi = xi,
#'  na.rm = FALSE, maxgap = 1)
#' # the interpolated values where only gap sizes of 1 are filled
#' df1[rowsinterps1, 3] <- interps1
#' # replace the NAs in V2 with the interpolated V2 values
#' df1
#'
#'
#'
#' # data.table time series example
#' dt1 <- structure(list(Time = structure(c(1395242100, 1395243000, 1395243900,
#'  1395244800, 1395245700, 1395256500, 1395257400, 1395258300, 1395259200,
#'  1395260100, 1395261000, 1395261900), class = c("POSIXct", "POSIXt"),
#'  tzone = "GMT"), V1 = c(1.6, 1.7, 1.7, 1.7, 1.7, 1.7, 1.6, 1.7, 1.7, 1.7,
#'  1.7, 1.7), V2 = c(2, 2.1, 2.1, NA, NA, 2.1, 2.1, NA, 2.3, NA, 2, 2.1)),
#'  .Names = c("Time", "V1", "V2"), row.names = c(NA, -12L), class =
#'  c("data.table", "data.frame"), sorted = "Time")
#'
#' rowsinterps2 <- which(is.na(dt1[, 3, with = FALSE] == TRUE))
#' # index of rows of x that have NA (to be interpolated)
#' xi <- as.numeric(dt1[rowsinterps2, Time])
#' # the Date-Times for V2 to be interpolated in numeric format
#' interps2 <- dt1[, na.interp1(as.numeric(Time), V2, xi = xi,
#'  na.rm = FALSE, maxgap = 1)]
#' # the interpolated values where only gap sizes of 1 are filled
#' dt1[rowsinterps2, `:=` (V2 = interps2)]
#' # replace the NAs in V2 with the interpolated V2 values
#' dt1
#'
#'
#' @import stats
#' @import zoo
#'
#' @export
# Sources 1 and 2 begin
na.interp1 <- function (x, y, xi = x, ..., na.rm = TRUE, maxgap = Inf) {

  na.interp1.vec <- function (x, y, xi = x, ...) {

  na <- is.na(y)

  yi <- approx(x[!na], y[!na], xi, ...)$y

if (maxgap < length(y)) {

# construct a series like y but with only gaps > maxgap
# (actual values don't matter as we only use is.na(ygap) below)

  ygap <- .fill_short_gaps(y, seq_along(y), maxgap = maxgap)
# construct y values at 'xi', keeping NAs from ygap
# (using indexing, as approx() does not allow NAs to be propagated)

  ix <- approx(x, seq_along(y), xi, ...)$y

  yx <- ifelse (is.na(ygap[floor(ix)] + ygap[ceiling(ix)]), NA, yi)

  yx

 } else {

  yi
}
}

if (!identical(length(x), length(y))) {

  stop("x and y must have the same length")
}

  x. <- as.numeric(x)

if (missing(xi) || is.null(xi)) xi <- x.

  xi. <- as.numeric(xi)

  y. <- y

  result <- if (length(dim(y.)) < 2) {

  na.interp1.vec(x., y., xi = xi., ...)

} else {

  apply(y., 2, na.interp1.vec, x = x., xi = xi., ...)

}

if (na.rm) {

  result <- na.trim(result, is.na = "all")

}

  result

}

# x = series with gaps
# fill = same series with filled gaps
.fill_short_gaps <- function(x, fill, maxgap) {

if (maxgap <= 0)

  return(x)

if (maxgap >= length(x))

  return(fill)

  naruns <- rle(is.na(x))

  naruns$values[naruns$lengths > maxgap] <- FALSE

  naok <- inverse.rle(naruns)

ifelse (naok, fill, x)
# Sources 1 and 2 end
}
