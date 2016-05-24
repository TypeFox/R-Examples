addtoedges <- function(profileEdges=blackbox.getOption("profileEdges"), locedge) {
  if(is.null(locedge)) return(profileEdges)
  oldnr <- nrow(profileEdges)
  if (is.null(oldnr)) oldnr <- 0
  if ( oldnr>1000 ) { message.redef("profileEdges growing dangerously...")}
  profileEdges <- rbind(profileEdges, locedge)
  nr <- nrow(profileEdges)
  if (is.null(nr)) nr <- 0
  if ( oldnr<1001 & nr>1000 ) {
    cat("(Reducing profileEdges)", "\n")
    profileEdges <- resetCHull(profileEdges, formats="vertices", redundant.mode="double")$vertices
  }
  blackbox.options(profileEdges=profileEdges)
  return(profileEdges) ## but the two return may not be used
}
