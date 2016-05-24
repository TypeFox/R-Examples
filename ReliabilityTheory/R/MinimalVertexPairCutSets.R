# Find collection of minimal (s,t) vertex cut sets
minimalVertexPairCutSets <- function(graph, start, terminate) {
  # Get all separators
  stCutSets <- minimal.st.separators(graph)
  # Double check separator (there was a bug in igraph early 0.6 releases)
  stCutSets <- subset(stCutSets, vapply(stCutSets, is.minimal.separator, TRUE, graph=graph))
  # Obviously, if s or t are in the cut set is it not of interest
  stCutSets <- stCutSets[vapply(lapply(stCutSets, match, match(c(start, terminate), V(graph)$name), 0), sum, 1)==0]
  # Check that the cut set *is* an s,t separator, since igraph gets all separators
  cutGraphs <- lapply(stCutSets, delete.vertices, graph=graph)
  reachable <- lapply(cutGraphs, subcomponent, start, "out")
  isSTsep <- rep(TRUE, length(stCutSets))
  for(j in 1:length(stCutSets)) { # Can't do nice vectorised stuff, because when deleting vertices igraph changes all the ids, so each check different
    if(match(terminate, V(cutGraphs[[j]])$name) %in% reachable[[j]]) {
      isSTsep[j] <- FALSE
    }
  }
  if(sum(isSTsep) < 1) return(NULL);
  stCutSets <- stCutSets[isSTsep]
}
