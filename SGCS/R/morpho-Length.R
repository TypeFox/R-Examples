#' Morphologicals: Relative boundary length of diluted pattern
#' 
#' Dilute pattern with b(o,r) and compute boundary length.
#' @param x Point pattern
#' @param r Vector of distances to estimate the function
#' @param ... Ignored.
#' 
#' @return
#' Only for 2D.
#' 
#' The default plotted curve, "u(r)" is U(r)/(2 * r * lambda * pi)
#' with U(r) the boundary length fraction.
#' 
#' Reduced sample border correction.
#' 
#' @export

morphoLength <- function(x, r, ...){
  ### prepare data
  x <- internalise_pp(x)
  if(x$dim==3) stop("Area fraction function only for 2d.")
  ### range
  r <- default_r(x, r)
  ### Distances for speed
  x$pairwise_distances <- pairwise_distances(x)
  ### Border distance for correction
  x$edgeDistances <- edge_distance(x)
  ### compute
  res <- .External("SGCS_morphoLength_c",
                   x,
                   r,
                   PACKAGE="SGCS"
  )
  ### Use spatstat:
  pp <- internal_to_ppp(x)
  windows <- lapply(r, erosion, w=pp$window)
  lambdas <- sapply(windows, function(w) intensity(pp[w]) )
  areas <- sapply(windows, area)
  U <- res/areas
  rU <- U/(2 * pi * r * lambdas)
  
  U[r==0] <- 0
  rU[r==0] <- 1
  # theoretical for Poisson
  l <- pi* intensity(pp) * r^2
  theo <- exp(-l)
  # make fv suitable
  A.final<-fv( data.frame(r=r, theo=theo, u=rU, U=U),
               argu = "r",
               alim = range(r),
               ylab = substitute(u(r), NULL),
               desc = c("distance argument r", "Poisson", "Relative border length fraction", "Border length fraction"),
               valu = "u",
               fmla = "cbind(u, theo)~r",
               fname="u"
  )
  
  A.final
}