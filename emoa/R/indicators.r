##
## pareto_utilities.r - Operators relating to pareto optimality
##
## Author:
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##

##' Scale point cloud
##'
##' Rescale all points to lie in the box bounded by \code{minval}
##' and \code{maxval}.
##'
##' @param points Matrix containing points, one per column.
##' @param minval Optional lower limits for the new bounding box.
##' @param maxval Optional upper limits for the new bounding box.
##' @return Scaled points.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' @export
normalize_points <- function(points, minval, maxval) {
  if (missing(minval))
    minval <- apply(points, 1, min)
  if (missing(maxval))
    maxval <- apply(points, 1, max)
  ## FIXME: This is ugly!
  (points - minval)/(maxval - minval)
}

##' Binary quality indicators
##'
##' Calculates the quality indicator value of the set of points given in
##' \code{x} with respect to the set given in \code{o}. As with all
##' functions in \code{emoa} that deal with sets of objective values
##' these are stored by column.
##'
##' @param points Matrix of points for which to calculate the indicator
##'   value stored one per column.
##' @param o Matrix of points of the reference set.
##' @param ref Reference point, if omitted, the nadir of the point sets
##'   is used.
##' @param ideal Ideal point of true Pareto front. If omited the ideal
##'   of both point sets is used.
##' @param nadir Nadir of the true Pareto front. If ommited the nadir
##'   of both point sets is used.
##' @param lambda Number of weight vectors to use in estimating the
##'   utility.
##' @param utility Name of utility function.
##' @return  Value of the quality indicator.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##'
##' @references
##'   Zitzler, E., Thiele, L., Laumanns, M., Fonseca, C., and
##'   Grunert da Fonseca, V (2003): Performance Assessment of
##'   Multiobjective Optimizers: An Analysis and Review. IEEE
##'   Transactions on Evolutionary Computation, 7(2), 117-132.
##'
##' @export
##' @rdname binary_indicator
hypervolume_indicator <- function(points, o, ref) {
  if (missing(ref))
    ref <- pmax(apply(points, 1, max), apply(o, 1, max))

  hvx <- dominated_hypervolume(points, ref)
  hvo <- dominated_hypervolume(o, ref)
  return(hvo - hvx)
}


##' @export
##' @rdname binary_indicator
epsilon_indicator <- function(points, o) {
  stopifnot(is.matrix(points), is.numeric(points),
            is.matrix(o), is.numeric(o))
  if (any(points < 0) || any(o < 0))
    stop("The epsilon indicator is only defined for strictly positive objective values.")
  
  .Call(do_eps_ind, points, o)
}

##
## R indicators:
##
r_indicator <- function(points, o, ideal, nadir, lambda, utility, summary) {
  ## (OME): Order of utility functions is important. It translates
  ## into the method number in the C code!
  utility.functions <- c("weighted sum", "Tchebycheff", "Augmented Tchebycheff")
  utility <- match.arg(utility, utility.functions)
  method <- which(utility == utility.functions)
  
  if (missing(ideal)) 
    ideal <- pmin(apply(points, 1, min), apply(o, 1, min))
  if (missing(nadir))
    nadir <- pmax(apply(points, 1, max), apply(o, 1, max))

  dim <- nrow(points)
  if (missing(lambda)) {
    lambda <- if (dim == 2) { 500 }
    else if (dim == 3) { 30  }
    else if (dim == 4) { 12  }
    else if (dim == 5) { 8   }
    else               { 3   }
  }
  
  ix <- .Call(do_r_ind, points, ideal, nadir,
              as.integer(lambda), as.integer(method))
  io <- .Call(do_r_ind, o, ideal, nadir,
              as.integer(lambda), as.integer(method))

  return(summary(ix, io))
}

##' @export
##' @rdname binary_indicator
r1_indicator <- function(points, o, ideal, nadir, lambda, utility="Tchebycheff")
  r_indicator(points, o, ideal, nadir, lambda, utility,
              function(ua, ur) mean(ua > ur) + mean(ua == ur)/2)

##' @export
##' @rdname binary_indicator
r2_indicator <- function(points, o, ideal, nadir, lambda, utility="Tchebycheff") 
    r_indicator(points, o, ideal, nadir, lambda, utility,
                function(ua, ur) mean(ur - ua))

##' @export
##' @rdname binary_indicator
r3_indicator <- function(points, o, ideal, nadir, lambda, utility="Tchebycheff") 
    r_indicator(points, o, ideal, nadir, lambda, utility,
                function(ua, ur) mean((ur - ua)/ur))

##' Unary R2 indicator
##' 
##' @param points Matrix of points for which to calculate the indicator
##'   value stored one per column.
##' @param weights Matrix of weight vectors stored one per column.
##' @param ideal Ideal point of true Pareto front. If omited the ideal
##'   of \code{points} is used.
##' @return Value of unary R2 indicator.
##' 
##' @export
##' @author Olaf Mersmann \email{olafm@@p-value.net}
unary_r2_indicator <- function(points, weights, ideal) {
  if (missing(ideal)) 
    ideal <- apply(points, 1, min)

  .Call(do_unary_r2_ind, points, weights, ideal)
}
