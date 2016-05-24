#' Triplet intensity function
#' 
#' Summarise the number of r-close triangles in a stationary and isotropic point pattern (2d,3d).
#'   
#' @param x Point pattern
#' @param r Vector of distances to estimate the function
#' @param ... ignored.
#' 
#' @details
#' Border correction is done via minus sampling.
#' 
#' See 
#' 
#' \strong{Schladitz, Baddeley: A Third order point process characteristic, SJS, vol 27, 657-671, 2000.}
#' 
#' for details.
#' 
#' @return 
#' \code{\link{fv}}-object.
#' 
#' @useDynLib SGCS
#' @import spatstat
#' @export

Tfun <- function(x, r, ...) {
  ### prepare data
  x <- internalise_pp(x)
  r <- default_r(x, r)
  ### border correction
  x$edgeDistances <- edge_distance(x)
  ### pairwise distances for speed
  x$pairwise_distances <- pairwise_distances(x)
  ### Compute:
  res <- .External("SGCS_Tfun_c",
                   x,
                   r,
                   PACKAGE="SGCS"
  )
  # scale:
  lambda <- x$n/x$area
  res <- res/lambda^2
  # theoretical
  theo <- if(x$dim == 2) 0.5*pi*(pi-3/4 * sqrt(3))*r^4 else 5/12 * pi^2 * r^6
  # make fv suitable
  C.final<-fv( data.frame(T=res, r=r, theo=theo),
               argu = "r",
               alim = range(r),
               ylab = substitute(K(r),NULL),
               desc = c("Triplet intensity", "Poisson", "range"),
               valu = "T",
               fmla = ".~r",
               fname="T"
  )
  
  C.final  
}
