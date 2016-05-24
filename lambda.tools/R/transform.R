# :vim set filetype=R
#' Force values into a set of bins
#'
#' This function quantizes data into a set of bins based on a 
#' metric function. Each value in the input is evaluated with each
#' quantization level (the bin), and the level with the smallest
#' distance is assigned to the input value.
#'
#' @section Usage:
#' quantize(x, bins=c(-1,0,1), metric=function(a,b) abs(a-b))
#'
#' @section Details:
#' When converting analog signals to digital signals, 
#' quantization is a natural phenomenon. This concept can be extended
#' to contexts outside of DSP. More generally it can be thought of
#' as a way to classify a sequence of numbers according to some
#' arbitrary distance function.
#' 
#' The default distance function is the Euclidean distance in 1 dimension.
#' For the default set of bins, values from (-infty, -.5] will map to -1.
#' The values from (-.5, .5] map to 0, and the segment (.5, infty) map to 1.
#' Regardless of the ordering of the bins, this behavior is
#' guaranteed. Hence for a collection of boundary points k and bins b,
#' where |b| = |k| + 1, the mapping will always have the form
#' (-infty, k_1] => b_1, (k_1, k_2] => b_2, ... (k_n, infty) => b_n.
#'
#' @name quantize
#' @param x A sequence
#' @param bins The bins to quantize into
#' @param metric The method to attract values to the bins
#' @return A vector containing quantized data
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{confine}}
#'
#' @examples
#' x <- seq(-2, 2, by=.1)  
#' quantize(x)
#'
#' quantize(x, bins=-1.5:1.5)
quantize(x, bins=c(-1,0,1), metric=function(a,b) abs(a-b)) %as% {
  bins <- sort(bins)
  ds <- sapply(bins, function(b) metric(x,b))
  ds <- onlyif(is.null(dim(ds)), t, ds)
  apply(ds,1, function(d) item(bins, which.min(d)))
}


#' Confine values to the given bounds
#'
#' Given a sequence this function confines the sequence values to within
#' the specified bounds. This behavior is equivalent to clipping in
#' digital signal processing.
#'
#' @section Usage:
#' confine(x, min.level, max.level) %::% numeric : numeric : numeric : numeric
#' confine(x, min.level=-1, max.level=1)
#'
#' @section Details:
#' The confine function can be thought of a transform that limits the
#' range of a sequence. Any values outside the range [min.level, max.level]
#' are adjusted to be exactly min.level or max.level.
#' 
#' Care should be taken when using this function as it is not always
#' a good idea to change the value of outliers. Sometimes it is better
#' to remove these values from a data set instead.
#'
#' @name confine
#' @param x A numeric vector
#' @param min.level The lower bound
#' @param max.level The upper bound
#' @return A sequence with values outside of [min.level, max.level]
#' clipped to those values
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{quantize}}
#'
#' @examples
#' confine(seq(-2,2, by=.1))
#'
confine(x, min.level, max.level) %::% numeric : numeric : numeric : numeric
confine(x, min.level=-1, max.level=1) %when% {
  min.level < max.level
} %as% {
  y <- ifelse(x < min.level, min.level, x)
  ifelse(y > max.level, max.level, y)
}


#' Partition a sequence into adjacent windows and apply a metric function 
#' to each window
#'
#' This function transforms a sequence into a rolling set of adjacent
#' windows separated by a pivot point. Each window is passed to a 
#' metric function that yields a scalar value. The result is effectively
#' a coordinate pair that represents the two adjacent windows.
#'
#' @section Usage:
#' partition(x, metric, radius) %::% . : Function : numeric : matrix
#' partition(x, metric=median, radius=10)
#' 
#' @section Details:
#' Many analysis approaches explore ways to reduce the dimensionality
#' of a data set to make it easier to model. The opposite situation is 
#' when there is not enough information in the data structure as is.
#' This circumstance requires a technique that can add dimensionality
#' to a data structure, which is what this function does.
#'
#' The idea is that a sequence can yield additional information by
#' comparing the neighborhoods around a given point. For this function,
#' a point is an index of the sequence. In the 1D case, given 
#' an index k and a radius r, the left neighbohood is defined by
#' [k-r+1, k] and the right neighborhood is defined by
#' [k+1, k+r]. The values associated with each neighborhood
#' are then applied to a metric function m: A^r -> R.
#' This output becomes the coordinate pair (left, right).
#'
#' At the edges of the sequence the above formalism is not completely
#' accurate. This is because at the edge, the neighborhood will be 
#' smaller than the radius, with a minimum size of 1.
#' Hence the first iteration on a sequence will yield a left neighborhood
#' of 1, while the right neighborhood will be [2, 1+r]. Whether this is
#' acceptable is case-specific.
#'
#' In the future, a wrap parameter might be included that would emulate
#' a loop instead of a sequence. This would be useful if a sequence
#' represented a stationary time series.
#'
#' @name partition 
#' @param x A sequence
#' @param metric A function that maps a vector to a real-valued scalar
#' @param radius The extent of the neighborhood about the index point
#' 
#' @return A length(x)-1 by 2 matrix where each row represents
#' the value of the metric applied to left and right neighborhoods
#' about an index point.
#'  
#' @author Brian Lee Yung Rowe
#'
#' @examples
#' partition(1:10, mean, radius=2)
#'
partition(x, metric, radius) %::% . : Function : numeric : matrix
partition(x, metric=median, radius=10) %when% {
  is.null(dim(x))
} %as% {
  f <- function(x,i) {
    c(left=metric(x[max(1,i-radius+1):i]),
      right=metric(x[(i+1):min(length(x),i+radius)]))
  }
  t(sapply(1:(length(x)-1), function(i) f(x,i)))
}


#' Segment a sequence into shifted versions of itself
#'
#' Create a shifted version of a sequence to make it easier to do
#' certain types of analysis.
#'
#' @section Usage:
#' segment(x, do.pad=FALSE)
#'
#' @section Details:
#' Segmenting sequences into offset versions of itself is useful for
#' detecting patterns in sequences. This approach is compatible with
#' a functional programming style as each row can then be passed to
#' a map-vectorized function for further processing.
#'  
#' The advantage over an iterative approach is that the map-vectorized
#' function can focus on a row-specific model independent of data
#' management mechanics like maintaining proper indices while iterating
#' over the sequence, as this is handled by segment.
#'
#' @name segment
#' @param x A vector
#' @param do.pad Whether the vector should be padded to contain 
#' the edges of the sequence
#' @return The return value is a data.frame with dimensions 
#' length(x) - 1 by 2 or length(x) + 1 by 2 if do.pad == TRUE.
#' A data.frame is used to support arbitrary types. For example,
#' using a Date vector will result in a numeric output, which is
#' inconvenient.
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{partition}} \code{\link{maprange}}
#' @note The segment function is a convenience and can be implemented
#' using the general functions partition and also maprange. If you want
#' more than two columns, use maprange.
#'
#' @examples
#' segment(1:10)
#'
#' # Notice how the ends of the sequence are given their own rows
#' segment(1:10, TRUE)
#'
#' # Emulate segment using partition
#' partition(1:10, function(x) x, 1)
#'
#' # Emulate segment using maprange
#' t(maprange(1:10, 2, function(x) x))
#'
#' # Create four shifted copies instead of two
#' maprange(1:10, 4, function(x) x)
segment(x, do.pad=FALSE) %when% {
  is.null(dim(x))
} %as% {
  x <- onlyif(do.pad, function(y) pad(y,1,1), x)
  data.frame(a=x[1:(length(x)-1)], b=x[2:length(x)])
}


#' Find contiguous ranges of a given value within a sequence
#'
#' Identify the index ranges for a given value in a sequence and return
#' the minimum and maximum values of the ranges.
#'
#' @section Usage:
#' range_for(target, x)
#'
#' @name range_for
#' @param target A value to find in x
#' @param x A vector
#' @return A data.frame where each row specifies the end points
#' of a contiguous range that contains the target value
#'
#' @author Brian Lee Yung Rowe
#'
#' @examples
#' # Find all contiguous ranges containing 2
#' x <- sample(c(1,2,2,2,3,4), 20, replace=TRUE)
#' range_for(2,x) 
#'
range_for(target, x) %when% {
  is.null(dim(x))
} %as% {
  y <- segment(x, TRUE)
  a <- y[,'a']
  b <- y[,'b']
  idx <- 1:nrow(y)
  idx.inf <- (is.na(a) | a != target) & (!is.na(b) & b == target)
  idx.sup <- (!is.na(a) & a == target) & (is.na(b) | b != target)
  data.frame(min=idx[idx.inf], max=idx[idx.sup]-1)
}


