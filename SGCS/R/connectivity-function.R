#' Connectivity function and cumulative connectivity function.
#' 
#' Amount of network connected pairs as a function of distance.
#' 
#' @param x Point pattern
#' @param r Vector of distances to estimate the function
#' @param R The radius for generating the network, as geometric graph.
#' @param h Smoothing parameter. h=0 and h>0 mean different things, see Details.
#' @param adjust Adjust h by this factor (default=1).
#' @param preGraph Precomputed network/graph, as a spatgraph-object. Alternative to R.
#' @param ... ignored.
#' 
#' @details
#' If h=0 we compute the cumulative version of the connecitivity function, 
#' corresponding to Ripley's K-function under the condition that the points 
#' in each pair must belong to the same component in an underlying network. 
#' 
#' The underlying network can be given, or it will be computed as a geometric graph with 
#' parameter 'R'. If given as 'preGraph', it must be a \pkg{spatgraphs}-object, with same
#' dimensions as the point pattern. 
#' 
#' If h>0: Compute the probability of a pair being in the same component given their distance is ~ r. 
#' Uses kernel smoothing with bandwidth h. 
#' 
#' Sensible defaults are computed for h and R if not given.
#' 
#' Border correction is done via translation correction. The bias is unknown as the network censoring 
#' is quite complex.
#' 
#' Theoretical values are unknown due to the graph conditioning.
#' 
#' @return 
#' fv-object, see \link{spatstat} for more. Theoretical values unknown.
#' 
#' @examples 
#' \dontrun{
#' x <- rMatClust(10, 0.1, 10)
#' plot(Cx<-confun(x,h=0, R=0.1))
#' 
#' # fit wrong model
#' ftho <- thomas.estpcf(x)
#' yf <- function()rThomas(ftho$par[1], ftho$par[2], x$n/ftho$par[1])
#' CC <- envelope(x, fun=confun, h=0, sim=yf, R=0.1)
#' C <- envelope(x, fun=confun, sim=yf, R=0.1)
#'
#' plot(CC)
#' plot(C)
#' }
#' 
#' @useDynLib SGCS
#' @import spatstat
#' @export

confun <- function(x, r, R, h, adjust=1, preGraph=NULL, ...) {
  ### prepare data
  x <- internalise_pp(x)
  ### default range for generating the component network
  if(missing(R)) R <- 1/(x$n/x$area)^(1/x$dim)
  ### default smoothing
  if(missing(h)) h <- 0.15*R
  h <- adjust * h
  ### range
  r <- default_r(x, r)
  
  ### Translation weights for correction
  x$weights <- translation_weights(x)
  ### Distances for speed
  x$pairwise_distances <- pairwise_distances(x)
  
  ### check preGraph
  if(!is.null(preGraph)) if(class(preGraph) != "sg") stop("preGraph is not of class 'sg' (see package 'spatgraphs')")
  ### Compute:
  res <- .External("SGCS_confun_c",
                   x,
                   r,
                   c(R,h), # function parameters
                   preGraph,
                   PACKAGE="SGCS"
                   )
  
  # scale if needed
  cumu <- NULL
  if(h==0){ 
    res <- res / (x$n/x$area)^2
    cumu <- "Cumulative"
  }
  # paramaeters
  pars <- paste0("(", ifelse(h==0, "", paste0("h=",h, ", ")), "R=", R, ")")
  
  
  # make fv suitable
  C.final<-fv( data.frame(r=r, C=res),
               argu = "r",
               alim = range(r),
               ylab = substitute(C(r), NULL),
               desc = c("distance argument r", paste(cumu, "Connectivity Function", pars)),
               valu = "C",
               fmla = ".~r",
               fname="C"
  )
  
  C.final
}

