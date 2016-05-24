# ------------------------------------------------------------------------------
# Function 'dissim'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# Last update: 24 May 2014
# Depends: -
# ------------------------------------------------------------------------------
dissim <- function(x, data, nb, adjust = FALSE, p2n.args, n2m.args, 
                   wVECT.args, v2n.args, verbose = FALSE) {

  if (!missing(x)) {
    if (inherits(x, "SpatialPolygons"))
      data <- suppressMessages(chksegdata(x, data))$data
    else {
      msg <- paste("'x' is neither \"SpatialPolygons\"",
                   "nor \"SpatialPolygonsDataFrame\"")
      stop(msg, call. = FALSE)
    }
  }

  if (ncol(data) > 2) {
    warning("'data' has more than two columns; only the first two are used",
            call. = FALSE)
    data <- data[,1:2]
  }
  if (any(data < 0))
    stop("negative value(s) in 'data'", call. = FALSE)
  colsum <- apply(data, 2, sum)
  if (any(colsum <= 0))
    stop("the sum of each column in 'data' must be > 0", call. = FALSE)

  out <- list(d = NA, dm = NA, dw = NA, ds = NA, user = NA)
  
  # Duncan and Duncan's index of dissimilarity
  b <- data[,1]/sum(data[,1])     # Blacks
  w <- data[,2]/sum(data[,2])     # Whites
  out$d <- as.vector(sum(abs(b-w))/2)
  
  if (!missing(x) & adjust) {
    userpkg <- .packages(all.available = TRUE)
    if ("spdep" %in% userpkg) {
      tmp <- tryCatch(.use.spdep(x, data, p2n.args, n2m.args, verbose),
                      error = function(e) print(e))
      if (is.numeric(tmp))
        out$dm <- out$d - tmp
      else if (verbose)
        message("failed to calculate D(adj)")
      
      if (all(c("spgrass6", "rgdal") %in% userpkg)) {
        tmp <- tryCatch(.use.spgrass6(x, data, wVECT.args, v2n.args, verbose),
                        error = function(e) print(e))
        if (is.numeric(tmp[1]))
          out$dw <- out$d - tmp[1]
        else if (verbose)
          message("failed to calculate D(w)")
        if (is.numeric(tmp[2]))
          out$ds <- out$d - tmp[2]
        else if (verbose)
          message("failed to calculate D(s)")
      }
    }
  }
  
  if (!missing(nb))
    out$user <- out$d - .d.adjust(data, nb)
  out
}
