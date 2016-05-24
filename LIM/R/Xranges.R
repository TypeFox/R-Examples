
##==============================================================================
##  Calculates ranges of inverse components (linp problem)
##==============================================================================

Xranges <- function(lim,...) {

  res   <- xranges (lim$A, lim$B,
                    lim$G, lim$H, ...)

  rownames(res) <- lim$Unknowns

  return(res)

}
