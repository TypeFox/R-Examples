#' Clustering function versio 2
#' 
#' Ratio of T and K
#' 
#' @param x Point pattern
#' @param r Vector of distances
#' @param correction Border correction. "none" or "border" (reduced window) supported.
#' @param scaled Scale with theoretical value?
#' @param ... ignored.
#' 
#' @details
#' 
#' Reduced border correction available.
#'
#' @return    
#' \code{\link{fv}}-object, see \link{spatstat} for more.
#' 
#' @examples
#' \dontrun{
#' en <- envelope(rcell(nx=15), fun=Rfun)
#' plot(en)
#' }
#' @useDynLib SGCS
#' @import spatstat
#' @export

Rfun <- function(x, r, correction="border", scaled=FALSE, ...) {
  ### prepare data
  x <- internalise_pp(x)
  ### range
  r <- default_r(x, r)
  
  ### Distances for speed
  x$pairwise_distances <- pairwise_distances(x)
  
  ### Border distances for correction
  correction_i <- correction %in% c("border","best")
  x$edgeDistances <- if(correction_i) edge_distance(x) else rep(max(r), x$n)
  
  
  ### Compute:
  res <- .External("SGCS_Rfun_c",
                   x,
                   r,
                   PACKAGE="SGCS"
  )
  
  Te <- if(x$dim==2) 0.5*pi*(pi-3*sqrt(3)/4)*r^4 else 5*pi^2*r^6/12
  Ke <- if(x$dim==2) pi*r^2 else pi*r^3 * 4/3
  # div 0
  i <- which(Ke==0)
  Ke[i] <- 1e-9
  # 
  lam <- x$n/x$area 
  lpr <- lam * Ke
  p  <- 1 - exp(-lpr) * (lpr+1) # p(K>0)
  theo <- (Te / Ke^2) * p 
  # if we scale away the Poisson value
  if(scaled){
    res <- res/theo
    res[r==0] <- Inf
    theo <- theo/theo
    theo[r==0] <- Inf
  }
  # make fv suitable
  c.final<-fv( data.frame(r=r, theo=theo, R=res),
               argu = "r",
               alim = range(r),
               ylab = substitute(R(r), NULL),
               desc = c("distance argument r", "Theoretical values unknown", "Ratio Function"),
               valu = "R",
               fmla = ".~r",
               fname="R"  
  )
  
  c.final
}