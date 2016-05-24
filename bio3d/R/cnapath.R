# Correlation network suboptimal path analysis
#
# Reference 
# Yen, J.Y. (1971) Finding the K Shortest Loopless Paths in a Network.
# Management Science. 17(11):712-716.

cnapath <- function(cna, from, to, k = 10, ncore = NULL, ...) {
  
  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) 
     stop("igraph package missing: Please install, see: ?install.packages")
  
  if(!inherits(cna, "cna")) 
     stop("Input cna is not a 'cna' object")

  ncore = setup.ncore(ncore)

  graph = cna$network

  # which path from the list is the shortest?
  select.shortest.path <- function(variants){
     return( which.min( unlist( lapply( variants, function(x){x$dist} ) ) ) )
  }

  # does a list contain this path?
  contains.path <- function(variants, variant){
     return( any( unlist( lapply( variants, function(x){ isTRUE(all.equal(x$path, variant)) } ) ) ) )
  }

  # first shortest path
  k0 <- igraph::get.shortest.paths(graph, from, to, output='both', ...)
  
  # if no shortest path found, network contains isolated parts.
  if(length(k0$vpath[[1]]) == 0) {
     cat("  No path found.\n", 
         "  Please check if the network contains isolated parts!\n\n", sep="")
     return(NULL)
  }

  # number of currently found shortest paths
  kk <- 1

  # All shortest paths are stored in container A in order
  dist = sum(igraph::E(graph)$weight[k0$epath[[1]]])
  A <- list(list(path=as.integer(k0$vpath[[1]]), epath=as.integer(k0$epath[[1]]), dist=dist))

  # All candidates are stored in container B
  B <- list()

  # For progress bar
  pb <- txtProgressBar(min=0, max=k, style=3)

  # until k shortest paths are found
  while(kk < k){
    # take last found shortest path
    last.path <- A[[length(A)]]
    
    tmpB <- mclapply(1:(length(last.path$path)-1), function(i) {    
       spurNode <- last.path$path[i]
       rootPath <- last.path$path[1:i]
       if(i==1) rootePath = NULL
       else rootePath = last.path$epath[1:(i-1)]

       # Remove edges that coincide with the next step from the spur node on
       # those shortest paths stored in A that share the same root path here
       g <- graph
       for(j in 1:length(A)) {
          if(length(A[[j]]$path) > i && isTRUE(all.equal(rootPath, A[[j]]$path[1:i]))) {
             nn = A[[j]]$path[i+1]
             ee = igraph::E(g)[igraph::'%--%'(spurNode, nn)]
             if(length(ee)>0) g <- igraph::delete.edges(g, ee)
          }
       }
       # Remove all edges that link to nodes on the root path (excluding the spur node)
       if(i > 1) {
          for(j in rootPath[-(length(rootPath))]) {
             ee = igraph::E(g)[from(j)] 
             if(length(ee)>0) g <- igraph::delete.edges(g, ee)
          }
       }

       # Suppress warnings because some nodes are intentionally isolated 
       spurPath <- suppressWarnings(igraph::get.shortest.paths(g, spurNode, to, output='both'), ...)

       if(length(spurPath$vpath[[1]]) > 0 ) {
          vpath = c(rootPath, as.integer(spurPath$vpath[[1]][-1]))
          if(!contains.path(B, vpath)) {
             spurPath$epath <- as.integer(igraph::E(graph, path=as.integer(spurPath$vpath[[1]])))
             epath = c(rootePath, spurPath$epath)
             return (list(path=vpath, epath = epath, dist = sum(igraph::E(graph)$weight[epath])) )
          }
       }
       NULL
    } ) 
    tmpB <- tmpB[ !sapply(tmpB, is.null) ]
    B <- c(B, tmpB)
    if(length(B) == 0) break
    
    # find shortest candidate
    sp <- select.shortest.path(B)

    # add to A, increase kk, remove shortest path from list of B
    A <- c(A, B[sp])
    kk <- kk + 1
    B <- B[-sp]

    setTxtProgressBar(pb, kk)
  }

  # stopped before reaching k paths
  if(kk < k) {
    setTxtProgressBar(pb, k)
    warning("Reaching maximal number of possible paths (", kk, ")")
  }
  close(pb)

  out <- list(path=lapply(A, "[[", "path"),  
              epath = lapply(A, "[[", "epath"), 
              dist = sapply(A, "[[", "dist"))
  class(out) <- c("cnapath", "list")
  return(out)
}

