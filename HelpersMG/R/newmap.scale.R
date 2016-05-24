#' newmap.scale Add Scale to Existing Unprojected Map
#' @title Add Scale to Existing Unprojected Map
#' @author See map.scale {maps}
#' @return The exact calculated scale is returned.
#' @param x Location of left end of distance gauge.
#' @param y Location of left end of distance gauge.
#' @param relwidth Proportion of width of display to be used for the scale. The default is 0.15.
#' @param metric If TRUE, the distance gauge will be in km, otherwise miles.
#' @param ratio If FALSE, the scale ratio of the map is not displayed.
#' @param col.line The color of lines for the gauge.
#' @param ... Further plotting parameters may be specified as for the command text().
#' @description Adds a scale to an existing map, both as a ratio and a distance gauge.
#' If x or y are not specified, this will be taken to be near the lower left corner of the map.
#' @examples
#' \dontrun{
#' library("maps")
#' library("HelpersMG")
#' map("world", "China")
#' newmap.scale(col.line = "red", col="blue")
#' }
#' @export


newmap.scale <- function (x, y, relwidth = 0.15, metric = TRUE, ratio = TRUE, 
                        col.line = "black", ...) 
{
  format.pretty <- function(x, digits = 2) {
    x = signif(x, 2)
    prettyNum(formatC(x, format = "fg", digits = digits), 
              big.mark = ",")
  }
  usr <- par("usr")
  if (missing(y)) 
    y <- (9 * usr[3] + usr[4])/10
  if (abs(y) >= 90) 
    warning("location of scale out of this world!")
  if (missing(x)) 
    x <- (9 * usr[1] + usr[2])/10
  cosy <- cos((2 * pi * y)/360)
  perdeg <- (2 * pi * (6356.78 + 21.38 * cosy) * cosy)/360
  scale <- (perdeg * 1e+05)/(2.54 * (par("pin")/diff(par("usr"))[-2])[1])
  if (metric) 
    unit <- "km"
  else {
    perdeg <- perdeg * 0.6213712
    unit <- "mi"
  }
  len <- perdeg * relwidth * (usr[2] - usr[1])
  ats <- pretty(c(0, len), n = 2)
  nats <- length(ats)
  labs <- as.character(ats)
  labs[nats] <- paste(labs[nats], unit)
  linexy <- matrix(NA, ncol = 2, nrow = 3 * nats)
  colnames(linexy) <- c("x", "y")
  cxy <- par("cxy")
  dy <- cxy[2] * par("tcl")
  dx <- ats[nats]/perdeg/(nats - 1)
  linexy[1, ] <- c(x, y)
  linexy[2, ] <- c(x, y + dy)
  for (i in 1:(nats - 1)) {
    linexy[3 * i, ] <- c(x + (i - 1) * dx, y)
    linexy[3 * i + 1, ] <- c(x + i * dx, y)
    linexy[3 * i + 2, ] <- c(x + i * dx, y + dy)
  }
  lines(linexy, col = col.line)
  text(x + ats/perdeg, y + dy - 0.5 * cxy[2], labs, adj = c(0.4, 
                                                            0.5), ...)
  if (ratio) 
    text(x, y + 0.5 * cxy[2], paste("scale approx 1:", format.pretty(scale), 
                                    sep = ""), adj = 0, ...)
  invisible(scale)
}
