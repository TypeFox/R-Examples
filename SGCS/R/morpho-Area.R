#' Morphologicals: Relative area fraction of diluted pattern
#' 
#' Dilute pattern with b(o,r) and compute area.
#' @param x Point pattern
#' @param r Vector of distances to estimate the function
#' @param ... Ignored.
#' 
#' @return
#' The default plotted curve, "af", is relative to lambda * pi * r^2. 
#' The component "AF" holds the area fraction.
#' 
#' @export

morphoArea <- function(x, r, ...){
  ### prepare data
  x <- internalise_pp(x)
  if(x$dim==3) stop("Area fraction function only for 2d.")
  ### range
  r <- default_r(x, r)
  ## use spatstat for now
  #   ### Distances for speed
  #   x$pairwise_distances <- pairwise_distances(x)
  #   ### Border distance for correction
  #   x$edgeDistances <- edge_distance(x)
  #   ### compute
  #   res <- .External("SGCS_morphoArea_c",
  #                    x,
  #                    r,
  #                    PACKAGE="SGCS"
  #   )
  #   
  ### Use spatstat:
  pp <- internal_to_ppp(x)
  res <- Fest(pp, r=r, correction="border")$rs
  windows <- lapply(r, erosion, w=pp$window)
  lambdas <- sapply(windows, function(w) intensity(pp[w]) )
  areas <- sapply(windows, area)
  AF <- res
  rAF <- AF/(lambdas * pi * r^2)
  
  # theoretical for Poisson
  l <- pi*x$n/x$area * r^2
  theo <- (1-exp(-l))/l
  # make fv suitable
  A.final<-fv( data.frame(r=r, theo=theo, rAF=rAF, AF=AF),
               argu = "r",
               alim = range(r),
               ylab = substitute(rAF(r), NULL),
               desc = c("distance argument r", "Poisson", "Relative area fraction", "Area fraction"),
               valu = "rAF",
               fmla = "cbind(rAF, theo)~r",
               fname="rAF"
  )
  
  A.final
}