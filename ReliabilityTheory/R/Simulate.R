# Simulate n system failure times from the given system
# system may be defined by a graph, cutsets or signature
simulateSystem <- function(system, n, rdens, ...) {
  if(class(system) == "igraph") {
    sig <- computeSystemSignature(system)
  } else if(class(system) == "list") {
    stop("Cutset system definition not supported yet.")
  } else if(class(system) == "numeric") {
    if(min(system) >= 0 && sum(system) == 1) {
      sig <- system
    } else {
      stop("Vector passed is not a signature.")
    }
  }
  i <- sample.int(length(sig), n, TRUE, sig)
  samps <- lapply(split(rdens(n*length(sig), ...), rep(1:n, each=length(sig))), sort)
  as.numeric(unlist(samps)[(1:n-1)*length(sig)+i])
}

