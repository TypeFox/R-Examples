isomorphicCutSets <- function(x, y, n) {
  if(length(x)!=length(y[[1]])) { return(FALSE); } # Different number of cutsets!
  if(prod(sort(unlist(lapply(x, length))) == sort(unlist(lapply(y[[1]], length)))) == 0) { return(FALSE); } # Different sizes in the cutsets!
  len <- length(y[[1]])
  for(i in 1:n) {
    l <- length(intersect(x, y[[i]]))
    if(len == l) {
      return(TRUE)
    }
  }
  return(FALSE)
}

coherentSystemsOfOrder <- function(n) {
  facN <- factorial(n)
  # Create the adjacency matrix, with two additional vertices for s/t
  adjacency <- matrix(0,ncol=n+2,nrow=n+2)
  # How many matrices do we try to be exhaustive?
  tot <- 2^(((n+2)*(n+2)-n-2)/2)-1
  res <- list()
  stseps <- list()
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
    # Check if the graph represents a coherent system <=> union of minimal separators is all nodes
    stsep <- minimalVertexPairCutSets(g, "s", "t")
      # Now, check the union of all minimal separators contains all nodes
    if(length(unique(unlist(stsep))) == n) {
      stsepPerms <- cutAndPathSetPerms(stsep, n)
      # Check for equality of cut sets -- note that we must check all permutations of the vertex labels too
#       eq <- lapply(stseps, function(x, y, n) {
#                               if(length(x)!=length(y[[1]])) { return(0); }
#                               len <- length(y[[1]])
#                               for(i in 1:n) {
#                                 l <- length(intersect(x, y[[i]]))
#                                 if(len == l) {
#                                   return(1)
#                                 }
#                               }
#                               return(0)
#                             }, stsepPerms, factorial(n))
      if(i > 1) {
        skipG <- FALSE
        for(j in (i-1):1) { # Do backwards -- we are most likely to be equivalent to a recently generated system
          if(isomorphicCutSets(stseps[[j]], stsepPerms, facN)) {
            skipG <- TRUE
            break;
          }
        }
        if(skipG) {
          next;
        }
      }
      # Ok, it's genuinely a new one
      res[[i]] <- list(graph=g, cutsets=stsep, signature=computeSystemSignature(g, stsep))
      stseps[[i]] <- stsep
      i <- i+1
    }
  }
  res
}

# sccsO2 <- coherentSystemsOfOrder(2)
# o <- order(unlist(lapply(sccsO2, function(x) { expectedSignatureLifetimeExp(x$signature)})))
# sccsO2 <- sccsO2[o]
# save(sccsO2, file="sccsO2.RData", compress=TRUE)

# sccsO3 <- coherentSystemsOfOrder(3)
# o <- order(unlist(lapply(sccsO3, function(x) { expectedSignatureLifetimeExp(x$signature)})))
# sccsO3 <- sccsO3[o]
# save(sccsO3, file="sccsO3.RData", compress=TRUE)

# sccsO4 <- coherentSystemsOfOrder(4)
# o <- order(unlist(lapply(sccsO4, function(x) { expectedSignatureLifetimeExp(x$signature)})))
# sccsO4 <- sccsO4[o]
# save(sccsO4, file="sccsO4.RData", compress=TRUE)

# sccsO5 <- coherentSystemsOfOrder(5)
# o <- order(unlist(lapply(sccsO5, function(x) { expectedSignatureLifetimeExp(x$signature)})))
# sccsO5 <- sccsO5[o]
# save(sccsO5, file="sccsO5.RData", compress=TRUE)

# sccsO6 <- coherentSystemsOfOrder(6)
# o <- order(unlist(lapply(sccsO6, function(x) { expectedSignatureLifetimeExp(x$signature)})))
# sccsO6 <- sccsO6[o]
# save(sccsO6, file="sccsO6.RData", compress=TRUE)
