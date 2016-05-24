# Function to compute the signature of an arbitrary system
# Define a graph of the system using the igraph package, with additional nodes named ``s'' and ``t'' for either end of the system diagram (see example below)
# Pass this to the graph object to the function to compute the signature
computeSystemSignature <- function(graph, cutsets=NULL, frac=FALSE) {
  # Function to check which order statistic caused failure
  cutNum <- function(x) {
    # Cycle through each order statistic ...
    for(i in 1:length(x)) {
      # ... and compare subset upto that order against all cut sets
      for(cut in cutsets) {
        if(sum(cut %in% x[1:i])==length(cut)) {
          # order statistic i caused failure here
          return(i)
        }
      }
    }
  }
  
  if(is.null(cutsets)) {
    # Find the minimal s,t separating vertices
    cutsets <- minimalVertexPairCutSets(graph, "s", "t")
  }
  # Translate igraph vertex IDs to naming scheme passed to graph.formula()
  cutsets <- lapply(cutsets, function(x) { as.numeric(V(graph)$name[x]) })
  
  # Compute signature
  s <- tabulate(as.vector(permn(1:(length(V(graph)$name)-2), cutNum), mode="numeric"), nbins=length(V(graph)$name)-2)
  if(frac) {
    denom <- factorial(length(V(graph)$name)-2)
    cat("s = ( ")
    for(val in s) {
      cd <- gcd(val, denom)
      cat(val/cd, "/", denom/cd, ", ", sep="")
    }
    cat(" )\n")
  } else {
    return(s/factorial(length(V(graph)$name)-2))
  }
}

## \textbf{Example}
## This is the bridged system of Figure \ref{fig:HardStructure}
#g <- graph.formula(s -- 1:3 -- 5 -- 2:4 -- t, 1 -- 2, 3 -- 4)
#computeSystemSignature(g)
## Output is:
#  2   3   4 
#0.2 0.6 0.2
## \textbf{NB:} $s_i=0$ if value for $i$ is not listed 

# Function to compute the signature of an arbitrary coherent network
# Define a graph of the system using the igraph package, with additional nodes named ``s'' and ``t'' for either end of the system diagram (see example below)
# Pass this to the graph object to the function to compute the signature
computeNetworkSignature <- function(graph, cutsets=NULL, frac=FALSE) {
  # Function to check which order statistic caused failure
  cutNum <- function(x) {
    # Cycle through each order statistic ...
    for(i in 1:length(x)) {
      # ... and compare subset upto that order against all cut sets
      for(cut in cutsets) {
        if(sum(cut %in% x[1:i])==length(cut)) {
          # order statistic i caused failure here
          return(i)
        }
      }
    }
  }
  
  if(is.null(cutsets)) {
    # Find the minimal s,t separating vertices
    cutsets <- minimalEdgeCutSets(graph, "s", "t")
  }
  # Translate igraph vertex IDs to naming scheme passed to graph.formula()
  #cutsets <- lapply(cutsets, function(x) { as.numeric(V(graph)$name[x]) })
  
  # Compute signature
  s <- tabulate(as.vector(permn(1:(length(as.numeric(E(graph)))), cutNum), mode="numeric"), nbins=length(as.numeric(E(graph))))
  if(frac) {
    denom <- factorial(length(as.numeric(E(graph))))
    cat("s = ( ")
    for(val in s) {
      cd <- gcd(val, denom)
      cat(val/cd, "/", denom/cd, ", ", sep="")
    }
    cat(" )\n")
  } else {
    return(s/factorial(length(as.numeric(E(graph)))))
  }
  
}
