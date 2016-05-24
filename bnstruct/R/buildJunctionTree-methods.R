#' @rdname build.junction.tree
#' @aliases build.junction.tree,InferenceEngine
#' @import igraph Matrix bitops methods
setMethod("build.junction.tree",
          c("InferenceEngine"),
          function(object, ...) {
            # Calculate junction tree of a given graph.
            # Input parameter is the adjacency matrix of a directed graph.
            
            dgraph <- dag(bn(object))
            
            graph <- moralization(dgraph)    # adj. matrix
            graph <- directed.to.undirected.graph(graph)   # adj. matrix
            graph <- triangulation(graph)   # adj. matrix
            ctout <- clique.tree(graph)     # (adj. matrix, list of lists of nodes)
            
            ctree <- ctout$clique.tree
            cs    <- ctout$cliques
            
            junction.tree(object)      <- ctree
            #num.nodes(object)          <- length(cs)
            jt.cliques(object)         <- cs
            # triangulated.graph(object) <- graph
            
            validObject(object)
            object
          }
)


directed.to.undirected.graph <- function(dg)
{
  # Take a NxN adjacency matrix representing a directed graph, remove edge directionality,
  # and return a NxN matrix representing the resulting undirected graph.
  # There is probably a better, native way for this...
  
  N <- nrow(dg)
  
  for ( i in 1:N )
  {
    for ( j in 1:N ) # cannot assume triangular matrices
    {
      if (dg[i,j] == 1 && dg[j,i] == 0)
      {
        #print(paste(i,j,"add bidirectional edge"))
        dg[j,i] <- 1
      }
    }
  }
  
  return(dg)
}

moralization <- function(graph)
{
  # Compute moral graph of the given NxN directed graph
  # by connecting non-adjacent nodes that share a child.
  
  # get number of nodes, create empty moral graph
  N <- nrow(graph)
  moral <- graph #matrix(0, N, N)
  
  # Iterate through nodes.
  # For each node, if it has 2+ parents, check if they are adjacent.
  # If not, add an edge in the moral graph.
  # Works assuming that graph[i,j] contains directed arc i->j.
  for (i in 1:N)
  {
    found <- which(graph[,i] == 1)
    
    if (length(found) > 1) # got the parents, check 'moral edge'
    {
      for (j in 1:(length(found)-1))
      {
        for (k in (j+1):length(found))
        {
          if (graph[found[j],found[k]] == 0 && graph[found[j],found[k]] == 0)
          {
            #print(paste(found[j],found[k],"share a child:",i))
            moral[found[j],found[k]] <- 1
            moral[found[k],found[j]] <- 1
          }
        }
      }
    } # end if
    
  } # for i
  
  # print("moral graph")
  # print(moral)
  
  #   # Finally, remove edge directionality in the original graph.
  #   # If a directed edge appears in the given graph, then add
  #   # the corresponding undirected edge to the moral graph.
  #   for (i in 1:N)
  #   {
  #     for (j in 1:N) # cannot assume triangular matrices
  #     {
  #       if (graph[i,j] == 1)
  #       {
  #         moral[i,j] <- 1
  #         moral[j,i] <- 1
  #       }
  #     }
  #   }
  
  return(moral)
}

triangulation <- function(graph)
{
  # Get an undirected graph, return a triangulated (chordal) graph,
  # i.e. a graph which has no cycles of length >3 without an edge (chord)
  # connecting two nodes that are not adjacent in the cycle.
  
  # create igraph object from adjacency matrix
  ig <- graph.adjacency(graph, "undirected", weighted=NULL, diag=TRUE,
                        add.colnames=NULL, add.rownames=NA)
  
  # All of the following could be replaced (I hope) with
  # > is.chordal(ig, newgraph=TRUE)$newgraph
  # but in the current version of igraph there is a fatal bug when ran in linux
  # Under windows it is reported to work, but I haven't tested it yet
  
  # for what I understand, the result is a list of the form a b c d ...
  # meaning that the fill-in edges are (a,b), (c,d), ...
  fill.in.edges <- is.chordal(ig, fillin=TRUE)$fillin
  
  #print("fill in edges")
  #print(fill.in.edges)
  
  # manually add the edges to my original graph
  i <- 1
  while(i <= length(fill.in.edges))
  {
    #print(paste("adding edge",fill.in.edges[i],fill.in.edges[i+1]))
    graph[fill.in.edges[i],fill.in.edges[i+1]] <- 1
    graph[fill.in.edges[i+1],fill.in.edges[i]] <- 1
    i <- i+2
  }
  
  return(graph)
}

clique.tree <- function(graph)
{
  # Compute clique tree of the given graph, i.e. a graph whose nodes
  # are the maximal cliques in the given graph, and the edges form the 
  # maximum spanning tree for the clique tree.
  # Returns the clique tree and the cliques.
  
  # create igraph object from adjacency list
  ig <- graph.adjacency(graph, "undirected", weighted=NULL, diag=TRUE,
                        add.colnames=NULL, add.rownames=NA)  
  
  # get list containing all the maximal cliques in the graph
  igraph::cliques(ig, min=NULL, max=NULL)
  cs <- maximal.cliques(ig)
  # print("---")
  # print(cs)
  num.cliques <- length(cs)
  # print(num.cliques)
  
  # for (i in 1:num.cliques)
  # {
  #   print("....")
  #   print(cs[i])
  #   cs[i] <- list(sort(unlist(cs[i])))
  # }
  
  # construct distance matrix for the clique tree:
  # - nodes are the cliques
  # - edges between cliques are associated with sepsets (the nodes they share)
  #     ctree[i,j] <- |sepset(C_i, C_j)|
  #
  # Also, a complementary graph is computed, in order to use the igraph
  # minimum.spanning.tree method to compute a maximum spanning tree
  ctree <- matrix(0, num.cliques, num.cliques)
  compl <- matrix(0, num.cliques, num.cliques)
  for (i in 1:num.cliques)
    for (j in 1:num.cliques)
    {
      # why the previous sorting (cs[i] <- list(sort(unlist(cs[i])))) has no effect?
      # print(paste(i,j,cs[i],cs[j],intersect(sort(unlist(cs[i])),sort(unlist(cs[j])))))
      # cat(unlist(cs[i])," ",unlist(cs[j])," ", length(intersect(sort(unlist(cs[i])),sort(unlist(cs[j])))), "\n")
      ctree[i,j] <- length(intersect(sort(unlist(cs[i])),sort(unlist(cs[j]))))
      compl[i,j] <- -ctree[i,j]#min(length(unlist(cs[i])),length(unlist(cs[j]))) - ctree[i,j]
    }
  #print(ctree)
  #print(compl)
  ig <- graph.adjacency(compl, "undirected", weighted=TRUE, diag=FALSE,
                        add.colnames=NULL, add.rownames=NA)  
  # compute max-spanning tree using the complementary graph, and get adjacency matrix
  ig.mst <- minimum.spanning.tree(ig)
  mst <- get.adjacency(ig.mst, type="both", attr=NULL, edges=FALSE)#, names=TRUE,)
  #sparse=getIgraphOpt("sparsematrices"))  
  # 'edges=TRUE' + get.edgelist may help with order...
  # print(mst)
  
  
  # convert mst adjacency matrix to a format we can use elsewhere
  mst <- matrix(data = mst, nrow = num.cliques, ncol = num.cliques, byrow = TRUE)
  
  # now manually insert the correct edge weights
  for (i in 1:num.cliques)
    for (j in 1:num.cliques)
      if (mst[i,j] == 1)
        mst[i,j] <- ctree[i,j]
  
  # print(mst)
  
  return(list("clique.tree" = mst, "cliques" = cs))
}
