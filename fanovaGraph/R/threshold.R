threshold <- function(graphlist, delta, scaled = TRUE, robust = FALSE) {
  d <- graphlist$d
  if (robust & ncol(graphlist$tii)!=4){
    warning("confidence intervals missing in argument 'graphlist', robust = TRUE ignored")
    robust <- FALSE
  }
  if (robust){
    valuesToCut <- graphlist$tii[,4]
  } else{
    valuesToCut <- graphlist$tii[,1]
  }
  if (scaled){
    valuesToCut <- valuesToCut/graphlist$V
  } 
  graphlistNew <- graphlist
  graphlistNew$tii[which(valuesToCut <= delta),] <- 0
  graphlistNew$tii.scaled[which(valuesToCut <= delta),] <- 0
  # estimate cliques
    E <- t(combn(d,2)[,graphlistNew$tii[,1] > 0])
    graphlistNew$cliques <- maximal.cliques(graph(as.vector(t(E)), d, FALSE))
  return(graphlistNew)
} 