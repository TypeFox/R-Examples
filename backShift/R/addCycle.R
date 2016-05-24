addCycle <- function(A, p, minCoef, maxCoef){
  # add cycles
  sample.func <- function(i) 
    sample(c(runif(i, min = minCoef, max = maxCoef), 
             runif(i, min = -maxCoef, max = -minCoef)), i)
  
  while(!hasCycles(A)){
    G <- graph.adjacency(A,mode="directed",weighted="a")
    nz <- which(abs(A) > 0, arr.ind = TRUE)
    sample.node <- sample(as.numeric(nz[,1]), 1)
    dfs.result <- graph.dfs(G, root = sample.node, unreachable = FALSE, dist = TRUE) 
    furthest.dis.from.sample.node <- max(dfs.result$dist, na.rm = TRUE) + 1
    furthest.from.sample.node <- which.max(dfs.result$dist)
    A[furthest.from.sample.node, sample.node] <- sample.func(1)
    diag(A) <- 0
  }
  list(A = A, sizeCycle = furthest.dis.from.sample.node)
}

