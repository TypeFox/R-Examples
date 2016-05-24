#' Lag-distance classes for variogram estimation
#' 
#' Computation of lag-distance classes for variogram estimation.
#' 
#' @param coords Data frame or matrix with the projected x- and y-coordinates.
#' 
#' @param n.lags Integer value defining the number of lag-distance classes that 
#' should be computed. Defaults to \code{n = 7}.
#' 
#' @param type Character value defining the type of lag-distance classes that
#' should be computed, with options \code{"equi"} (equidistant) and \code{"exp"}
#' (exponential). Defaults to \code{type = "exp"}.
#' 
#' @param cutoff Numeric value defining the fraction of the diagonal of the
#' rectangle that spans the data (bounding box) that should be used to set the
#' maximum distance up to which lag-distance classes should be computed. 
#' Defaults to \code{cutoff = 0.5}, i.e. half the diagonal of the bounding box.
#' 
#' @param base Numeric value defining the base of the exponential expression 
#' used to create exponentially spaced lag-distance classes. Used only when 
#' \code{type = "exp"}. Defaults to \code{base = 2}, i.e. the width of the 
#' rightmost lag-distance classes is equal to half the diagonal of 
#' \code{cutoff}, and so on.
#' 
#' @param zero Numeric value setting the minimum pair-wise separation distance 
#' that should be used to compute the lag-distance classes. Defaults to 
#' \code{zero = 0.0001}.
#' 
#' @param count Should the number of points (\code{"points"}) or point-pairs 
#' (\code{"pairs"}) per lag-distance class be computed? Defaults to 
#' \code{count = "pairs"}.
#' 
#' @return Vector of numeric values with the lower and upper boundaries of the 
#' lag-distance classes. The number of points or point-pairs per lag-distance
#' class is returned as an attribute.
#' 
#' @author Alessandro Samuel-Rosa <\email{alessandrosamuelrosa@@gmail.com}>
#' @seealso \code{\link[spsann]{optimPPL}}
#' @concept variogram
#' @references
#' Truong, P. N.; Heuvelink, G. B. M.; Gosling, J. P. Web-based tool for expert
#' elicitation of the variogram. \emph{Computers and Geosciences}. v. 51, p.
#' 390-399, 2013.
#' @export
#' @examples
#' data(meuse, package = "sp")
#' lags_points <- vgmLags(coords = meuse[, 1:2], count = "points")
#' lags_pairs <- vgmLags(coords = meuse[, 1:2], count = "pairs")
# FUNCTION - MAIN ##############################################################
vgmLags <-
  function (coords, n.lags = 7, type = "exp", cutoff = 0.5, base = 2, 
            zero = 0.001, count = "pairs") {
    
    # Check if suggested packages are installed
    pkg <- c("SpatialTools")
    id <- !sapply(pkg, requireNamespace, quietly = TRUE)
    if (any(id)) {
      pkg <- paste(pkg[which(id)], collapse = " ")
      stop(paste("Package(s) needed for this function to work but not",
                 "installed: ", pkg, sep = ""), call. = FALSE)
    }
    
    # Compute cutoff
    if (class(coords) %in% c("matrix", "data.frame", "array")) {
      cutoff <- sqrt(
        sum(apply(apply(coords[, 1:2], 2, range), 2, diff) ^ 2)) * cutoff
    } else {
      message("'coords' should be a data frame with the projected coordinates")
    }
    
    # Compute the boundaries of the lag-distance classes
    n_pts <- nrow(coords)
    lags <- switch(
      type,
      equi = { # Equidistant
        seq(zero, cutoff, length.out = n.lags + 1)
      },
      exp = { # Exponential
        idx <- base ^ c(1:n.lags - 1)
        c(zero, rev(cutoff / idx))
      }
    )
    
    # Count the number of points or point-pairs per lag-distance class
    dm <- SpatialTools::dist1(as.matrix(coords))
    ppl <- switch (
      count,
      pairs = { # Point-pairs per lag-distance class
        diff(sapply(
          1:length(lags), function (i) 
            length(which(dm <= lags[i]))) - n_pts) * 0.5
      },
      points = { # Points per lag-distance class
        sapply(1:n.lags, function (i)
          length(unique(c(
            which(dm > lags[i] & dm <= lags[i + 1], arr.ind = TRUE)))))
      })
    
    # Output with attributes
    a <- attributes(lags)
    a$count <- ppl
    attributes(lags) <- a
    return (lags)
  }
