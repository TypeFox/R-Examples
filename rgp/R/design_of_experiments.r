## design_of_experiments.R
##   - Utility functions for creating experiment designs
##     (i.e. samples of n-dimensional spaces)
##
## RGP - a GP system for R
## 2010-2012 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Patrick Koch, Olaf Mersmann and Joerg Stork
## Latin Hypercube Sampling Code by Chrisitian Lasarczyk
## released under the GPL v2
##

##' Create a regular grid design matrix
##'
##' Returns a n = length(points)**dimension times m = dimension matrix containing
##' the coordinates of sample points from a hypervolume of the given dimension.
##' Points are sampled in a grid defined by the vector points.
##'
##' @param dimension The number of columns in the design matrix to create.
##' @param points A vector of points to sample at in each dimension.
##' @return The regular grid design matrix.
##'
##' @export
gridDesign <- function(dimension, points = seq(from = 0.0, to = 1.0, length.out = 10)) {
  s = length(points)
  D <- matrix(0.0, nrow = s ** dimension, ncol = dimension)

  for (i in 1:nrow(D)) {
    for (j in 1:ncol(D)) {
      index <- (ceiling(i / s ** (j - 1)) - 1) %% s + 1
      D[i, j] <- points[index] 
    }
  }

  return (D)
}

##' Create a normalized design matrix
##'
##' Produces a normalized design and calculates the minimal distance 
##' if required. Returns a design is a matrix with \code{dim} columns and
##' \code{size} rows.
##' 
##' @param dimension Dimension of the problem (will be no. of columns of the result matrix).
##' @param size Number of points with that dimension needed. (will be no. of rows of the result matrix).
##' @param calcMinDistance Indicates whether a minimal distance should be calculated.
##' @return List \code{L} consists of a matrix and nd (if required) a minimal distance.
##'
##' @export
normalizedDesign <- function(dimension, size, calcMinDistance = FALSE) {
	step <- 1 / size
	D <- replicate(dimension, sample(0:(size - 1), size) * step + runif(size) * step)

	if (calcMinDistance) {
		minDistance <- min(dist(D))
	} else {
		minDistance <- NA
  }

	list(design = D,  minDistance = minDistance)
}

##' Create a latin hypercube design (LHD)
##'
##' Produces a LHD matrix with \code{dimension} columns and \code{size}
##' rows.
##'
##' @param dimension Dimension of the problem (will be no. of columns of the result matrix).
##' @param size Number of design points, defaults to \code{max(11 * dimension,
##'   1 + 3 * dimension + dimension * (dimension - 1) / 2 + 1)}.
##' @param lowerBounds Numeric vector of length \code{dimension} giving lower bounds
##'   for sampling, defaults to \code{c(0.0, ...)}.
##' @param upperBounds Numeric vector of length \code{dimension} giving upper bounds
##'   for sampling, defaults to \code{c(1.0, ...)}.
##' @param retries Number of retries, which is the number of trials to find a design
##'   with the lowest distance, default is \code{2 * dimension}.
##' @return A LHD matrix. 
##'
##' @export
latinHypercubeDesign <- function(dimension,
                                 size = max(11 * dimension, 1 + 3 * dimension + dimension * (dimension - 1) / 2 + 1),
                                 lowerBounds = replicate(dimension, 0.0),
                                 upperBounds = replicate(dimension, 1.0),
                                 retries = 2 * dimension) {
	# min distance does not have to be calculated if there is only one try 
	bestNormalizedD <- normalizedDesign(dimension, size, calcMinDistance = retries > 0)
	if (retries > 0) {
		for (i in 1:retries) {
			candidateD <- normalizedDesign(dimension, size, calcMinDistance = TRUE)
			## maximize minimal distance
			if (candidateD$minDistance > bestNormalizedD$minDistance) {
				bestNormalizedD <- candidateD 
      }
		}
	}

  # scale normalized matrix to bounds
  D <- bestNormalizedD$design
  for (i in 1:dimension) {
    D[, i] <- lowerBounds[i] + bestNormalizedD$design[, i] * (upperBounds[i] - lowerBounds[i])
  }

	return (D)
}
