summary.rsm <-
function(object,...){
  Y = object
  if(class(Y) != "rsm") stop("Entry must be from class rsm")

  vec <- rep(-Inf, max(Y$Klist))
  
  for(K in Y$Klist) vec[K] <- Y$output[[K]]$lower
  vec <- vec[-c(1:(Y$Klist[1]-1))]
  
  cat(paste("Initial settings:\n", Y$N, "vertices \n", Y$R, "subgraphs\n", Y$C, "relations types\n"))
  
  cat("\n The optimal number of cluster is K = ", Y$K_star, "; \n with the lower bound equal: ", max(vec), "\n")
  
}
