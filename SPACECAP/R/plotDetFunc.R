
# Produces a plot of the detection function.

# Not exported, no sanity checks, no documentation

# Calculates and displays the minimum distance between a trap and the edge of the
#  state space and the probability of detection at that distance.

# Requires: the e2dist function.

# Arguments:
# lam0 and sigma are point estimates of the detection parameters.
# habGrid is a 2-column matrix with the coordinates of all the pixels
#  making up the state space. This should include bad habitat as well
#  as good.
# trapGrid is a 2-column matrix with the trap locations.
# detFunc specifies the detection function to plot, either half-normal (HN)
#  or negative exponential (NE).

# Returns: a vector with the minimum distance and the detection probability 
#  at that distance.

# Note: This generates a "Note" in R CMD check, 
#   "no visible binding for global variable 'x'",
#   which is due to non-standard use of 'x' by the graphics::curve
#   function. It is NOT a problem.
utils::globalVariables("x") ## This sorts the issue of the global variable [AMG]

plotDetFunc <- function(lam0, sigma, habGrid=NULL, trapGrid=NULL,
    detFunc=c("HN", "NE")) {

  detFunc <- match.arg(detFunc)
  if (detFunc == "HN") {
    DF <- function(x, lam0, sigma) {lam0 * exp(-x^2/(2*sigma^2))}
    main <- "Half normal detection function"
  } else {
    DF <- function(x, lam0, sigma) {lam0 * exp(-x/sigma)}
    main <- "Negative exponential detection function"
  }

  ### Find minimum distance between traps and edge of state space
  minDistToEdge <- NULL
  if (!is.null(habGrid) && !is.null(trapGrid)) {
    try({      # interPix is huge, sometimes not enough memory.
      # Find the minimum distance between pixels:
      interPix <- e2dist(habGrid, habGrid)
      diag(interPix) <- NA
      pixWidth <- min(interPix, na.rm=TRUE)
      # Find the number of neighbours within this distance (+20%) for each pixel
      numnb <- colSums(interPix < (pixWidth * 1.2), na.rm=TRUE)
      if (sum(numnb > 4) == 0) {      # should not exceed 4
        edge <- habGrid[numnb < 4, ]  # Coordinates of the pixels on the edge
        minDistToEdge <- min(e2dist(trapGrid, edge))
      }
    }, silent=TRUE)
  }

  ### Do the basic plot
  curve(DF(x, lam0, sigma), 0, max(minDistToEdge*1.1, sigma*4),
    ylim=c(0, lam0*1.05), main=main, lwd=2, xaxs='i', yaxs='i', bty='l',
    xlab="Distance between trap and home range center",
    ylab="Probability of detection") # or "Expected number of captures" ?
    
  ### Add the minimum distance to edge information
  lamMin <- NULL
  if(!is.null(minDistToEdge)) {
    lamMin <- DF(minDistToEdge, lam0, sigma)
    message <- paste("Minimum distance between trap\nand edge of state space:\n",
      signif(minDistToEdge, 3),
      "\n\nProbability of detection here:\n", signif(lamMin, 3))
    arrows(minDistToEdge, lam0, minDistToEdge, lamMin, col='red', length=0.1, lwd=2)
    text(minDistToEdge, lam0*0.5, message, pos=2, col='red')  
  }
  
  return(c(minDistToEdge=minDistToEdge, lamMin=lamMin))
}

