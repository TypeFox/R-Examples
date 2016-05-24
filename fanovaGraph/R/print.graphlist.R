print.graphlist <- function(x,...){
  cat("Overall variance\n")
  print(x$V)
  cat("\nCliques\n")
  cl <- x$cliques
  for (i in 1:length(cl)) cat(paste(paste(cl[[i]], collapse=""),"\n"))
  cat("\nScaled total interaction indices\n")
  print(x$tii.scaled)
  }