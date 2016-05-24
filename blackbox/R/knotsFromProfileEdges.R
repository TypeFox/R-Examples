## profileEdges has accumulated edges encountered during all profile computations and may be very redundant
## but they may all be in on plane and convhulln will fail (); this is where redundant is still useful
# alternatively, use lower level linear programming assuggested by http://stackoverflow.com/questions/25399417/algorithm-intersecting-polytope-and-half-line
knotsFromProfileEdges <- function(lrthreshold) {
  profedges <- NULL
  profileEdges <- blackbox.getOption("profileEdges")
  if(! is.null(nrow(profileEdges))) {
    liks <- apply(profileEdges, 1, purefn, testhull=F) ## possibly many low likelihood values here
    profedges <- profileEdges[which(liks>lrthreshold), , drop=FALSE] ## retains those with high enough lik
    if(nrow(profedges)>0) {
      profedges <- q2d(redundant(d2q(cbind(0, 1, profedges)), representation="V")$output[, -c(1:2), drop=FALSE])
      colnames(profedges) <- colnames(profileEdges)
    }
  }
  return(profedges)
}

