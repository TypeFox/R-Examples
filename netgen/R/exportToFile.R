#' Exports a network to an proprietary format.
#'
#' The format used is similar to the TSPlib format (see \code{\link{exportToTSPlibFormat}}),
#' but it saves not only the point coordinates. It also saves the arrival times
#' of dynamic customers.
#'
#' @param x [\code{Network}]\cr
#'   Network to export.
#' @param filename [\code{character(1)}]\cr
#'   File name.
#' @param digits [\code{integer(1)}]\cr
#'   Round coordinates to this number of digits. Default is 2.
#' @return Nothing
#' @export
exportToFile = function(x, filename, digits = 2L) {
  assertClass(x, c("Network"))
  assertCharacter(filename, len = 1L, any.missing = FALSE)
  assertInt(digits, lower = 0L, upper = 10L, na.ok = FALSE)

  df = as.data.frame(x, include.extra = TRUE)

  out = paste("NAME: ", x$name, "\n", sep = "")
  out = paste(out, "COMMENT: ", x$comment, "\n", sep = "")
  out = paste(out, "DIMENSION: ", ncol(x$coordinates), "\n", sep = "")
  out = paste(out, "N_NODES: ", getNumberOfNodes(x), "\n", sep = "")
  out = paste(out, "N_CLUSTERS: ", getNumberOfClusters(x), "\n", sep = "")
  out = paste(out, "N_DEPOTS: ", getNumberOfDepots(x), "\n", sep = "")
  out = paste(out, "LOWER: ", x$lower, "\n", sep = "")
  out = paste(out, "UPPER: ", x$upper, "\n", sep = "")
  out = paste(out, "DATA_SECTION", sep = "")

  # round coordinates
  df$x1 = round(df$x1, digits)
  df$x2 = round(df$x2, digits)

  if (!is.null(x$arrival.times)) {
    arrival.times = round(x$arrival.times, digits)
    if (hasDepots(x)) {
      arrival.times = c(rep(0, getNumberOfDepots(x)), arrival.times)
    }
    df = cbind(df, data.frame("arrival.time" = arrival.times))
  }
  write(x = out, file = filename)
  suppressWarnings({write.table(df, file = filename, sep = ",",
    append = TRUE, # append to header stuff
    row.names = FALSE,
    quote = FALSE
  )})
}
