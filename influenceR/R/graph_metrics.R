# graph_metrics.R: R code for vertex importance metrics.
# AUTHOR: Simon Jacobs <sdjacobs@uchicago.edu>
# LICENSE: GPLv2


#' Convert a CSV file to an igraph graph object.
#'
#' The first column should be sources, the second should be targets.
#'
#' @param fname A filename
#' @return An igraph graph object built from the filename.
#' 
#' @examples 
#' \dontrun{ig.csv <- csv.to.igraph("edgelist.csv") }
#' 
#' @export
csv.to.igraph <- function(fname) {
    x <- utils::read.csv(fname) # this may be dangerous because of users' settings.
                         # See: http://r-pkgs.had.co.nz/r.html
    el <- as.matrix(x[c(1,2)])
    if(!is.character(el))
      el <- apply(el, 2, as.character)
    
    igraph::graph.edgelist(el, directed=F)
}

#' Vertex Betweenness centrality measure.
#'
#' The Betweenness centrality score of a node u is the sum over all pairs s,t of the
#' proportion of shortest paths between s and t that pass through u. This 
#' function allows the use of either the SNAP betweenness implementation (default), or 
#' the igraph betweenness function. The SNAP version makes use of OpenMP for 
#' parallelization, and may be faster in some circumstances.
#'
#' @references \url{http://snap-graph.sourceforge.net/}
#'
#' @param g The igraph object to analyze
#' @param snap True to use the SNAP betweenness code, False to use igraph::betweenness
#' @return A numeric vector with the betweenness centrality score for each vertex
#'
#' @examples
#' ig.ex <- igraph::erdos.renyi.game(100, p.or.m=0.3) # generate an undirected 'igraph' object
#' betweenness(ig.ex) # betweenness scores for each node in the graph
#'
#' @export
betweenness <- function(g, snap=T) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
# 1/2 the values of our betweenness code, which is because this is UNDIRECTED for real
  if (!snap)
    return(igraph::betweenness(g))
  
  el <- igraph::get.edgelist(g, names=F)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el))
  m <- as.integer(length(el)/2) # TODO: for directed too?
  
  vals <- .Call("snap_betweenness_R", el_i, n, m, PACKAGE="influenceR")
  vals[vals<2^-128] <- 0
  vals[is.nan(vals)] <- 0
  names(vals) <- igraph::V(g)$name
  vals
}

#' Compute a KPP-Pos set for a given graph.
#'
#' @description 
#' The "Key Player" family of node importance algorithms (Borgatti 2006) involves the selection
#' of a metric of node importance and a combinatorial optimization strategy to
#' choose the set S of vertices of size k that maximize that metric.
#' 
#' @details 
#' This function implements KPP-Pos, a metric intended to identify k nodes which
#' optimize resource diffusion through the network. We sum over all vertices
#' not in S the reciprocal of the shortest distance to a vertex in S.
#' 
#' According to Borgatti, a number of off-the-shelf optimization algorithms may
#' be suitable to find S, such as tabu-search, K-L, simulated annealing, or
#' genetic algorithms. He presents a simple greedy algorithm, which we excerpt
#' here:
#' 
#'  \enumerate{
#'    \item Select k nodes at random to populate set S
#'    \item Set F = fit using appropriate key player metric.
#'    \item For each node u in S and each node v not in S:
#'      \itemize{\item DELTAF = improvement in fit if u and v were swapped}
#'    \item Select pair with largest DELTAF
#'      \itemize{
#'        \item If DELTAF <= [tolerance] then terminate
#'        \item Else, swap pair with greatest improvement in fit and set F = F + DELTAF
#'      }
#'    \item Go to step 3.
#' }
#' 
#' Our implementation uses a different optimization method which we call
#' stochastic gradient descent. In tests on real world data, we found that
#' our method discovered sets S with larger fits in less computation time.
#' The algorithm is as follows:
#' 
#' \enumerate{
#'  \item Select k nodes at random to populate set S
#'  \item Set F = fit using appropriate key player metric (KPP-Pos in our case)
#'  \item Get a new state:
#'  \itemize{
#'    \item Pick a random u in S and v not in S.
#'    \item F' = fit if u and v were swapped
#'    \item If F' > F, swap u and v in S. Else, repeat step 3. (Alternatively, if a positive value is given for the `prob' parameter, a swap will be accepted with a small probability regardless of whether it improves the fit).
#'  }
#'  \item If F' - F < tolerance or our maximum computation time is exceeded, return S. Else, go to step 3.
#' }
#'     
#' This implementation uses OpenMP (if available on the host system) so that
#' multiple workers can explore the solution space in parallel. After a given
#' of time, the workers synchronize their sets S to the one which maximizes
#' the metric.
#'
#' @references \url{http://www.bebr.ufl.edu/sites/default/files/Borgatti\%20-\%202006\%20-\%20Identifying\%20sets\%20of\%20key\%20players\%20in\%20a\%20social\%20networ.pdf}
#'
#' @param g The igraph object to analyze.
#' @param k The size of the KP-set
#' @param prob probability of accepting a state with a lower value
#' @param tol tolerance within which to stop the optimization and accept the current value
#' @param maxsec The total computation budget for the optimization, in seconds
#' @param roundsec Number of seconds in between synchronizing workers' answer
#' @return a vector with the vertex number of each vertex in the selected set S.
#'
#' @examples
#' ig.ex <- igraph::erdos.renyi.game(100, p.or.m=0.3) # generate an undirected 'igraph' object
#' keyplayer(ig.ex, k=10) # key-player set consisting of 10 actors
#'
#' @export
keyplayer <- function(g, k, prob = 0.0, tol = 0.0001, maxsec = 120, roundsec = 30) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  
  if (roundsec > maxsec)
    roundsec <- maxsec
  
  el <- igraph::get.edgelist(g, names=F)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el))
  m <- as.integer(length(el)/2)

  converges <- vector("integer", 1) # just allocate space for an integer

  s <- .Call("snap_keyplayer_R", el_i, n, m, as.integer(k), prob, tol, as.integer(maxsec), as.integer(roundsec), converges, PACKAGE="influenceR")
  
  if (converges == 1)
    warning("Maximum computation time (arg 'maxsec') exceeded!")
  
  igraph::V(g)[which(s>0)]
}

#' Valente's Bridging vertex measure.
#'
#' Edges that reduces distances in a network are important structural bridges. Here we implement Valente and Fujimoto's metric,
#' where a node's bridging score is the average decrease in cohesiveness if each of
#' its edges were removed from the graph.
#' 
#' @references \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2889704/}
#'
#' @param g The igraph object to analyze.
#' @return A numeric vector with the bridging score for each vertex
#'
#' @examples
#' ig.ex <- igraph::erdos.renyi.game(100, p.or.m=0.3) # generate an undirected 'igraph' object
#' bridging(ig.ex) # bridging scores for each node in the graph
#' 
#' @export
bridging <- function(g) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  el <- igraph::get.edgelist(g, names = F)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el_i))
  m <- as.integer(length(el_i) / 2)
  
  x <- .Call("snap_bridging_R", el_i, n, m, as.integer(FALSE), as.integer(0), PACKAGE = "influenceR")
  names(x) <- igraph::V(g)$name
  x
}

#' Burt's Effective Network Size and Constraint index.
#' The next two functions below provide ways to measure the actors' access to structural holes in a network. Structural holes 
#' "provide opportunities to broker connections between people" (Burt 2008).

#' @param g The igraph object to analyze.
#' @return A numeric vector with the effective network size for each vertex
#'
#' @examples
#' ig.ex <- igraph::erdos.renyi.game(100, p.or.m=0.3) # generate an undirected 'igraph' object
#' ens(ig.ex) # Effective Network Size scores for each node in the graph
#'
#' @references \url{http://faculty.chicagobooth.edu/ronald.burt/research/files/NNappB.pdf}
#' @export
ens <- function(g) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  A <- igraph::get.adjacency(g)   # This will be sparse, which is great.
  S <- Matrix::crossprod(A)       # S[i,j] = # of shared neighbors between i,j
  Q <- A * S              # Q[i,j] = # of shared neighbors if i and j are neighbors, 0 else
  qsum <- Matrix::rowSums(Q)
  deg <- Matrix::rowSums(A)
  ens <- deg - (qsum / deg)
  ens[is.nan(ens)] <- 0 # If a vertex has no neighbors, make its ENS 0
  names(ens) <- igraph::V(g)$name
  ens
}

#' Burt's Constraint Index.
#'
#' The igraph package provides an implementation of Constraint; this is an alternate implementation.
#'
#' @param g The igraph object to analyze.
#' @param v vertices over which to compute constraint (default to all)
#' @return A numeric vector with the constraint score for each vertex in v
#' 
#' @examples
#' ig.ex <- igraph::erdos.renyi.game(100, p.or.m=0.3) # generate an undirected 'igraph' object
#' constraint(ig.ex) # constraint scores for each node in the graph
#'
#' @export
constraint <- function(g, v=igraph::V(g)) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  
  process_sparse <- function(A, Ai, deg) {
    M <- methods::as(A, 'TsparseMatrix')
    x <- .Call("process_sparse_R", M@i, M@j, M@x, Ai, deg, Matrix::nnzero(M), PACKAGE = "influenceR")
    M@x <- x
    M
  }
  
  A <- igraph::get.adjacency(g, sparse=T)
  n <- dim(A)[1]
  deg <- Matrix::rowSums(A)

  constraint_i <- function(i) {
    # process sparse does this: jq <- drop0(t(A*A[,i]) * A[,i]); jqd <- drop0(jq * deg)
    jqd <- process_sparse(A, A[i, ], deg)
    
    jqd <- Matrix::drop0(jqd)
    jqd@x <- (1 / jqd@x) * (1 / deg[i])
         
    Sj <- Matrix::colSums(jqd)
  
    idx <- as.numeric(igraph::neighbors(g, i))
    Sj[idx] <- Sj[idx] + (1 / deg[i])
    
    Sj2 <- Sj * Sj
    sum(Sj2)
  }
  
  vals <- sapply(v, constraint_i)
  names(vals) <- v$name
  vals
}

