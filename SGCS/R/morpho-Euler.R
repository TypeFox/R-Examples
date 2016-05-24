#' Morphologicals: Euler number of dilated pattern
#' 
#' Dilate pattern with b(o,r) and compute the Euler number.
#' 
#' @param x Point pattern
#' @param r Vector of distances to estimate the function
#' @param ... Ignored.
#' 
#' @return
#' Only for 2D.
#' 
#' Reduced sample border correction.
#' 

morphoEuler <- function(x, r, ...){
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
  res <- .External("SGCS_morphoEuler_c",
                   x,
                   r,
                   PACKAGE="SGCS"
  )
  ### Use spatstat:
  pp <- internal_to_ppp(x)
  windows <- lapply(r, erosion, w=pp$window)
  ns <- sapply(windows, function(w) pp[w]$n )
  lambda <- intensity(pp)
  E <- res/ns  / lambda
  
  #E[r==0] <- 0
  # theoretical for Poisson
  l <- pi* lambda * r^2
  theo <- (1-l)*exp(-l)
  # make fv suitable
  A.final<-fv( data.frame(r=r, theo=theo, E=E),
               argu = "r",
               alim = range(r),
               ylab = substitute(E(r), NULL),
               desc = c("distance argument r", "Poisson", 
                        "Euler number"),
               valu = "E",
               fmla = ".~r",
               fname="E"
  )
  
  A.final
}
