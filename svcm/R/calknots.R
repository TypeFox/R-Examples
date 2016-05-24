## NAME:
##
##   calknots
##
## PURPOSE:
##
##   calculate extended knot partition
##
## CATEGORY:
##
##   spline
##
## USAGE:
##
##   k = calknots(coord, degree, nrknots, delta=FALSE)
##
## INPUTS:
##
##   vnumber :  number of voxels
##   vsize   :  size of voxels
##   degree  :  integer value indicating the degree of the spline.
##   nrknots :  an integer number of inner knots of the desired partition.
##   delta   :  logical, indicating whether the equidistance 'delta' between
##              adjacent knots is to be returned or the vector of knot
##              positions.
##
## VALUE:
##
##   either the value of the equidistance between adjacent knots or the numeric
##   vector of equidistant knot positions of the extended partition (= inner
##   knots plus 'degree' knots).
##
## NOTE:
##
##   The range [0, vnumber*vsize] determines the uniform spreading of the inner
##   knots. The final partition, that will be returned, consists of
##   nrknots + 2*degree knots.
##   If then passed to splineDesign(), package 'splines', a design matrix with
##   (nrknots + 2*degree) - (degree + 1) = (nrknots + degree - 1) basis
##   functions is generated.
## 
## EXAMPLES:
##
##   k <- calknots(vnumber = 15, vsize = 1.739, degree = 3, nrknots = 5)
##   k <- calknots(15, 1.739, 3, 5, delta=TRUE)
##   library(splines)
##   B <- splineDesign(k, 1:15, ord = 4)
##
## HISTORY:
##
##   July 11, 2005    S.Heim: written;
##   Oct  31, 2005    S.Heim: replaced equidistant knot setting on the ad-hoc
##                            and somewhat arbitrary extended coordinate range
##                            by 1% (Eilers and Marx) with equidistant knot
##                            setting on the full range [0, max(dim)*vsize].
##                            This became necessary for proper increment of the
##                            resolution.
##--
calknots <- function(vnumber, vsize, degree, nrknots, delta=FALSE) {
  
  lower <- 0
  upper <- vnumber * vsize
  dist <- (upper - lower)/(nrknots - 1)
  knotvec <- seq(lower - degree * dist, upper + degree * dist, by = dist)

  if (delta == FALSE) {
    invisible(knotvec)
  } else {
    invisible(dist)
  }
  
}
