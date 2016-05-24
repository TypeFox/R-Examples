coherentNetworksOfOrder <- function(n) {
  # Create the adjacency matrix, with two additional vertices for s/t
  adjacency <- matrix(0,ncol=n+2,nrow=n+2)
  # How many matrices do we try to be exhaustive?
  tot <- 2^(((n+2)*(n+2)-n-2)/2)-1
  res <- list()
  graphs <- list()
  i <- 1
  progress <- 0.0
  for(m in 1:tot) {
    if(round(m*100/tot, 1) != progress) { progress <- round(m*100/tot, 1); cat("\r", progress, "%   "); }
    # Make the graph
    adjacency[upper.tri(adjacency)] <- c(digitsBase(m, 2, ((n+2)*(n+2)-n-2)/2))
    # Exclude possibility of direct s,t connection
    if(adjacency[1,n+2] == 1) next;
    g <- graph.adjacency(adjacency, mode="upper")
    # Unless there is only one connected component we're not interested
    if(length(subcomponent(g, 1)) != n+2) next;
    # Ok, so we can name the vertices
    V(g)$name <- c("s",1:n,"t")
    # Check if the graph represents a coherent system <=> union of minimal cuts is all edges
    stsep <- minimalEdgeCutSets(g, "s", "t")
    # Now, check the union of all minimal cuts contains all edges
    if(length(unique(unlist(stsep))) == sum(c(digitsBase(m, 2, ((n+2)*(n+2)-n-2)/2)))) {
#      iso <- lapply(graphs, graph.count.isomorphisms.vf2, graph2=g, vertex.color1=c(1,rep(2, n),1), vertex.color2=c(1,rep(2, n),1))
#      if(sum(unlist(iso))>0) next;
      if(i > 1) {
        skipG <- FALSE
        for(j in (i-1):1) { # Do backwards -- we are most likely to be equivalent to a recently generated system
          if(graph.count.isomorphisms.vf2(graphs[[j]], g, vertex.color1=c(1,rep(2, n),1), vertex.color2=c(1,rep(2, n),1))) {
            skipG <- TRUE
            break;
          }
        }
        if(skipG) {
          next;
        }
      }
      # Ok, it's genuinely a new one
      res[[i]] <- list(graph=g, cutsets=stsep, signature=computeNetworkSignature(g, stsep))
      graphs[[i]] <- g
      i <- i+1
    }
  }
  res
}

# cnO2 <- coherentNetworksOfOrder(2)
# o <- order(unlist(lapply(cnO2, function(x) { expectedSignatureLifetimeExp(x$signature)})))
# cnO2 <- cnO2[o]
# save(cnO2, file="cnO2.RData", compress=TRUE)
# plot(cnO2[[1]]$graph, vertex.color=c(2,3,3,2))

# cnO3 <- coherentNetworksOfOrder(3)
# o <- order(unlist(lapply(cnO3, function(x) { expectedSignatureLifetimeExp(x$signature)})))
# cnO3 <- cnO3[o]
# save(cnO3, file="cnO3.RData", compress=TRUE)
# plot(cnO3[[1]]$graph, vertex.color=c(2,3,3,3,2))

# cnO4 <- coherentNetworksOfOrder(4)
# o <- order(unlist(lapply(cnO4, function(x) { expectedSignatureLifetimeExp(x$signature)})))
# cnO4 <- cnO4[o]
# save(cnO4, file="cnO4.RData", compress=TRUE)
