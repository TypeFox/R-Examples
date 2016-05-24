# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

##################################################################################
### Node-level measure
### mreach.degree
#' Compute the M-reach Degree Centrality Score in a Netwrok
#'
#' \code{mreach.degree} computes the size of reachable nodes from a particular node within M steps.
#' M-reach degree centrality generalizes the \code{\link[sna]{degree}} centrality
#' by delimiting specific neighborhoods.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param node Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If not specified, scores for all nodes will be reported.
#'
#' @param M Number indicating the maximum geodistance between two nodes,
#' above which the two nodes are considered disconnected.
#' M hence defines the reachable set. The default is \code{Inf}.
#'
#' @param binary Logical scalar. If \code{TRUE}, the adjacency matrix is binarized,
#' and thus M essentially means steps. If \code{FALSE}, the edge values are
#' considered.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures and is the default.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'
#' @details The interprtation of the measure in binary and weighted adjacency matrix
#' are slightly different. In binary networks, the reachable set of nodes is defined by nodes that
#' are reachable within M steps. In weighted networks, the reachable set
#' is defined by nodes that are reachable within geodistance M.
#'
#' @return A vector indicating the outdegree, indegree, or total-degree
#' mreach.degree score of the chosen node; or a data frame containing all
#' the above information.
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua and Yu-Hsin Liu (2016). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' Working Paper, Indiana Univeristy.\cr
#'
#' Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
#' version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr
#'
#' Csardi, G and Nepusz, T (2006). "The igraph software package for complex network research."
#' InterJournal, Complex Systems 1695. \url{http://igraph.org} \cr
#'
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix,
#' # where edge values represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#'
#' # List the 2-reach degree scores for every node where W is binarized
#' mreach.degree(W,M=2,cmode="all",large=FALSE)
#'
#' @seealso
#' \code{\link[sna]{geodist}};
#' \code{\link[igraph]{shortest.paths}};
#' \code{\link{mreach.closeness}};
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#' @export

mreach.degree=function(adj.matrix, node, M = Inf, binary = TRUE, cmode = "all", large = TRUE, geodist.precomp = NULL){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  rownames(d) <- c(colnames(d))

  # create geodistance matrix
  if (is.null(geodist.precomp)){
    if (isTRUE(large)){
      if (isTRUE(binary)){
        g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = NULL)
      }else{
        g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = TRUE)
      }
      distances=igraph::shortest.paths(g, mode = "out")
    }else{
      distances=sna::geodist(d,ignore.eval = binary)$gdist
    }
  }else{
    distances=geodist.precomp
  }

  colnames(distances) <- c(colnames(d))
  rownames(distances) <- c(colnames(d))

  # set M constraint using distances matrix
  if(missing(M)){
    distances = distances
  }else{
    for(i in 1:ncol(distances)){
      for(j in 1:nrow(distances)){
        if(distances[i,j]>M) distances[i,j]=Inf
      }
    }
  }
  # calculate the score
  diag(distances)=Inf
  weights=1/distances
  weights[weights>0]=1
  outdegree=apply(weights,1,sum) # column vector, outdegree
  indegree=apply(weights,2,sum) # row vector, indegree
  total= outdegree + indegree

  if(missing(node)){
    kppos=cbind(outdegree, indegree, total)
  }else{
    kppos=c(outdegree[node], indegree[node], total[node])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
  }
  # reports
  if (cmode=="outdegree"){
    outdegree[node]
  }else if (cmode=="indegree"){
    indegree[node]
  }else if (cmode=="total"){
    total[node]
  }else{
    kppos
  }
}

###############################################################################
### mreach.closeness
#' Compute the M-reach Closeness Centrality Score in a Netwrok
#'
#' \code{mreach.closeness} refines the \code{\link{mreach.degree}} centrality by
#' using the (inverse) geodistance as weights.
#' The edge values should be properly interpreted as distances.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param node Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If not specified, scores for all nodes will be reported.
#'
#' @param M Number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected.
#' M hence defines the reachable set. The default is \code{Inf}.
#'
#' @param binary Logical scalar. If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures and is the default.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'
#' @details \code{mreach.closeness} refines the \code{\link{mreach.degree}} centrality
#' by using the (inverse) geodistance as weights, just as \code{\link[sna]{closeness}}
#' centrality refines \code{\link[sna]{degree}} centrality.
#' It captures the degree centrality when M is properly set (e.g. M=1 in a binarized network).
#' It captures the Gil-Schmidt power index (Gil and Schmidt, 1996)
#' and the cohesion centrality (Borgatti, 2006) when M is sufficiently large
#' (unconstrained). The normalization factor takes care of non-binary
#' edge values. Also note that the geodistance matrix does
#' not necessarily to be symmetric.
#'
#' @return A vector indicating the outdegree, indegree, or total-degree
#' cohesion score of the chosen player; or a data frame containing all
#' the above information. Note that the outdegree and indegree scores are normalized
#' to [0,1]. This means that the total-degree score is between [0,2].
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua and Yu-Hsin Liu (2016). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' Working Paper, Indiana Univeristy.\cr
#'
#' Borgatti, Stephen P. (2006). "Identifying Sets of Key Players in a Network."
#' \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr
#'
#' Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
#' version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr
#'
#' Csardi, G and Nepusz, T (2006). "The igraph software package for complex network research."
#' InterJournal, Complex Systems 1695. \url{http://igraph.org} \cr
#'
#' Gil, J and Schmidt, S (1996). "The Origin of the Mexican Network of Power."
#' Proceedings of the International Social Network Conference, Charleston, SC, 22-25.\cr
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix, where edge values
#' # represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # Transform the edge value to distance interpretaion
#' A <- W
#' A[W!=0] <- 1/W[W!=0]
#'
#' # List all types of 2-reach closeness scores for every node
#' mreach.closeness(A,M=2,cmode="all",large=FALSE)
#'
#' @seealso
#' \code{\link[sna]{geodist}};
#' \code{\link[igraph]{shortest.paths}};
#' \code{\link{mreach.degree}};
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#'
#' @export

mreach.closeness=function(adj.matrix, node, M = Inf, binary = FALSE, cmode = "all", large = TRUE, geodist.precomp = NULL){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  rownames(d) <- c(colnames(d))

  # create matrix of geodesic distances
  if (is.null(geodist.precomp)){
    if (isTRUE(large)){
      if (isTRUE(binary)){
        g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = NULL)
      }else{
        g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = TRUE)
      }
      distances=igraph::shortest.paths(g, mode = "out")
    }else{
      distances=sna::geodist(d,ignore.eval = binary)$gdist
    }
  }else{
    distances=geodist.precomp
  }

  colnames(distances) <- c(colnames(d))
  rownames(distances) <- c(colnames(d))

  # set threshold constraint using distances matrix
  if(missing(M)){
    distances = distances
  }else{
    for(i in 1:ncol(distances)){
      for(j in 1:nrow(distances)){
        if(distances[i,j]>M) distances[i,j]=Inf
      }
    }
  }

  # calculate the score
  diag(distances)=Inf
  weights=1/distances
  sum.out=apply(weights,1,sum) # column vector, outdegree
  sum.in=apply(weights,2,sum) #  row vector, indegree
  outdegree = sum.out/((ncol(distances)-1)*max(weights)) # normalize to [0,1]
  indegree = sum.in/((ncol(distances)-1)*max(weights)) # normalize to [0,1]
  total=outdegree + indegree
  if(missing(node)){
    kppos=cbind(outdegree, indegree, total)
  }else{
    kppos=c(outdegree[node], indegree[node],total[node])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
  }
  # reports
  if (cmode=="outdegree"){
    outdegree[node]
  }else if (cmode=="indegree"){
    indegree[node]
  }else if (cmode=="total"){
    total[node]
  }else{
    kppos
  }
}

#####################################################################################
### fragment

#' Compute the Fragmentation Centrality Score in a Netwrok
#'
#' \code{fragment} measures the extent of fragmentation of a network after a
#' set of nodes is removed from the network. The more fragmented the residual network is, the more central a node is.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If there are multiple players,
#' use \code{c(index1,index2,...)}.
#' If not specified, scores for all nodes will be reported.
#'
#' @param M Number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected.
#' M hence defines the reachable set. The default is \code{Inf}.
#'
#' @param binary Logical scalar. If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. The default is \code{FALSE}.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'
#' @details A natural way to apply the fragmentation centrality is in the
#' context of counter-terrorism, as shown in Borgatti (2006).
#' The measure uses geodistances to compute the fragmentation level of the
#' residual network, and thus edge values should be properly adjusted to
#' distance interpretation. The fragmentation centrality is not directional
#' as edge values are counted aggregately at the network level.
#
#'
#' @return Vector indicating fragment score(s) of the chosen player(s).
#' Score is normalized to [0,1].
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua and Yu-Hsin Liu (2016). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' Working Paper, Indiana Univeristy.\cr
#'
#' Borgatti, Stephen P. 2006. "Identifying Sets of Key Players in a Network."
#' \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr
#'
#' Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
#' version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr
#'
#' Csardi, G and Nepusz, T (2006). "The igraph software package for complex network research."
#' InterJournal, Complex Systems 1695. \url{http://igraph.org} \cr
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix, where edge values
#' # represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # Transform the edge value to distance interpretaion
#' A <- W
#' A[W!=0] <- 1/W[W!=0]
#'
#' # List the fragmentation centrality scores for every node
#' fragment(A)
#'
#' @seealso
#' \code{\link[sna]{geodist}};
#' \code{\link[igraph]{shortest.paths}};
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#'
#' @export

fragment=function(adj.matrix, nodes, M=Inf, binary=FALSE, large=TRUE, geodist.precomp=NULL){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  # geodesic distances
  if (is.null(geodist.precomp)){
    if (isTRUE(large)){
      if (isTRUE(binary)){
        g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = NULL)
      }else{
        g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = TRUE)
      }
      distances=igraph::shortest.paths(g, mode = "out")
    }else{
      distances=sna::geodist(d,ignore.eval = binary)$gdist
    }
  }else{
    distances=geodist.precomp
  }

  diag(distances)=Inf # replace the diagonal with Infs
  weights=1/distances # take the reciprocal of distances
  m=max(weights)

  if(missing(nodes)){
    s <- c(1:ncol(adj.matrix))
    for(k in 1:ncol(adj.matrix)){
      d <- data.matrix(adj.matrix, rownames.force = NA)
      d <- d[-k,-k]
      # create matrix of geodesic distances
      distances=sna::geodist(d,ignore.eval = binary)$gdist
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
      # cauculate the score
      diag(distances)=Inf
      weights=1/distances
      sum=sum(weights) # sum for the matrix elements, sum is a scalar
      s[k]=1-sum/(ncol(distances)*(ncol(distances)-1)*m)
    }
    s<-t(s)
    rownames(s)<-"fragment"
    t(s)
  }else{
    # remove objective nodes
    ColToDelete <- nodes
    d <- d[,-ColToDelete]
    RowToDelete <- t(nodes)
    d <- d[-RowToDelete,]

    # create matrix of geodesic distances
    distances=sna::geodist(d,ignore.eval = binary)$gdist

    # set threshold using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }
    # cauculate the score
    diag(distances)=Inf
    weights=1/distances
    sum=sum(weights) # sum for the matrix elements, sum is a scalar
    1-sum/(ncol(distances)*(ncol(distances)-1)*m)
  }
}


################################################################################
### diffusion
#' Compute the Diffusion Centrality Score in a Network
#'
#' \code{diffusion} measures player's ability to disseminate information through all the
#' possible paths. For each path from i to j there is a reaching probability
#' P_{ij}, which is specified in the inputted adjacency matrix.
#'
#'
#' @param adj.matrix Matrix indicating the probability matrix.
#'
#' @param node Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If not specified, scores for all nodes will be reported.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. In the first iteration, the adjacency matrix
#' is as the input. In the nth iteration, the adjacency matrix becomes
#' the input adjacency matrix to the power of n. By default, T is the network size.
#'
#' @return A vector indicating the defusion centrality score(s) of
#' the chosen player(s).
#'
#' @details The diffusion centrality measures the expected number of information receivers from a particular node (Banerjee et.al. 2013). The measure can approximate the degree, Katz-Bonacich, or
#' eigenvector centrality when proper parameters are chosen. See Banerjee et.al. (2014) for details and proofs.
#'
#' In its original parametrization (Banerjee et.al. 2013), P=q*g, where q is a measure of the information passing probability and g the adjacency matrix. For simplication and consistency with other centrality measures, the current packages asks users to input the probability matrix P directly.
#' With information on q and the adjacency matrix, the probability matrix P can easily be calculated by their product.
#'
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua and Yu-Hsin Liu (2016). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' Working Paper, Indiana Univeristy.\cr
#'
#' Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2013):
#' "Diffusion of Microfinance," \emph{Science}, Vol. 341. p.363\cr
#'
#' Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2014):
#' "Gossip: Identifying Central Individuals in a Social Network,"
#' Working Paper.\cr
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix, where edge values
#' # represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # Transform the edge value to probability interpretaion
#' P <- W *0.2
#'
#' # List the diffusion centrality score for every node
#' diffusion(P, T = 2)
#'
#' @seealso
#' \code{\link[matpow]{matpow}};
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#'
#'
#' @export


diffusion=function(adj.matrix, node, T=ncol(adj.matrix)){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  s <- 0
  for (t in 1:T){
    s=s+matpow::matpow(d,t)$prod1
  }
  # sum for each node (row), sum is a column vector, outdegree
  diffusion=apply(s,1,sum)
  if(missing(node)){
    score <- t(diffusion)
    rownames(score)<-"diffusion"
    t(score)
  }else{
    diffusion[node]
  }
}


##################################################################################
##################################################################################
### Group-level measure (also able to measure the node level)
### contract function
#' Group the Chosen Players in a Network
#'
#' \code{contract} combines selected nodes into one large pseudo-node and provides a reduced network.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#' The inputted adjacency matrix for the diffusion centrality should
#' be properly transfomred to the probability interpretation.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix.
#'
#' @param method Indication of which grouping criterion should be used.\cr
#' \code{method="min"} indicates the "minimum" criterion (edge values as distances).\cr
#' \code{method="max"} indicates the "maximum" criterion (edge values as non-cummulative strengths).\cr
#' \code{method="add"} indicates the "addition" criterion (edge values as cummulative strengths).\cr
#' \code{method="union"} indicates the "union" criterion (edge values as probability).\cr
#' The default is "min". See details for examples.
#'
#' @details Minimum Criterion: the edge value between a group and an outside node
#' is measured as the minimal value among all the (nonzero) edge values between
#' any node in the group and the outside node. Suggested if edge values are
#' interpreted as distances.\cr
#' \emph{Example: suppose node A to C has distance 2 and B to C has distance 1,
#' then according to the minimum criterion, the distance between C and
#' the merged set AB is 1. Note that if B and C are not connected,
#' the algorithm takes the distance between A and C to describe
#' the distance between AB and C.}
#'
#' Maximun Criterion: the edge value between a group and an outside node
#' is measured as the maximal value among all the (nonzero) edge values between
#' any node in the group and the outside node. Suggested if edge values are
#' interpreted as non-cummulative strengths. \cr
#' \emph{Example: we keep using the above example, but the figure now indicates
#' the strength of tie. According to the maximum criterion, the strength of tie
#' between AB and C is 2.}
#'
#' Addition Criterion: the edge value between a group and an outside node
#' is measured as the sum of all the edge values between any node in the group
#' and the outside node. Suggested if edge values are as cummulative strengths. \cr
#' \emph{Example: according to the addition criterion, the strength of tie between
#' AB and C is 3}
#'
#' Union Criterion: the edge value between a group and an outside node is measured
#' as the probability that there is at least one path connecting the group with
#' the outside node. Suggested if edge values are as probability. \cr
#' \emph{Example: suppose A has probability 0.2 to reach C and B has probability
#' 0.5 to reach C, then C can be reached from merged AB with probability
#' 1-(1-0.2)*(1-0.5)=0.6 according to the union criterion.}
#'
#'
#' @return A new adjacency matrix after contracting the chosen nodes (named
#' \code{set}).
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua and Yu-Hsin Liu (2016). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' Working Paper, Indiana Univeristy.\cr
#'
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix, where edge values
#' # represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # If the strength is believed to be non-accumulative for a group of nodes,
#' # it is proper to use the "maximum" criterion to contract node 2 and 3
#' contract(W,c(2,3),"max")
#'
#' # Transform the edge value to probability interpretaion
#' P <- W *0.2
#'
#' # Contract node 2 and 3 using the "union" criterion as it is proper for
#' # probability matrix input
#' contract(P,c(2,3),"union")
#'
#' @seealso
#' \code{\link{kpcent}};
#' \code{\link{kpset}}
#'
#' @export

contract=function(adj.matrix, nodes, method=c("min","max","union","add")){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  colnames(d) <- 1:ncol(d)
  rownames(d) <- 1:nrow(d)

  if (missing(method)){
    # contract the objective nodes; row sum for outdegree; column sum for indegree
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      max <- max(d)+1
      d <- replace(d,d==0,max)
      set <- d[,ncol(d)]
      for(i in 1:nrow(d)){
        set[i]=min(d[i,nodes])
      }
      set.r <- d[ncol(d),]
      for(j in 1:ncol(d)){
        set.r[j]=min(d[nodes,j])
      }
      d <- replace(d,d==max,0)
      set <- replace(set,set==max,0)
      set.r <- replace(set.r,set.r==max,0)
    }
  }else if (method=="add"){
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      set <- rowSums(d[,nodes])
      set.r <- colSums(d[nodes,])
    }
  }else if (method=="union"){
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      set <- rowSums(d[,nodes])
      for(i in 1:nrow(d)){
        set[i]=1-prod(1-d[i,nodes])
      }
      set.r <- colSums(d[nodes,])
      for(j in 1:ncol(d)){
        set.r[j]=1-prod(1-d[nodes,j])
      }
    }
  }else if(method=="max"){
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      set <- d[,ncol(d)]
      for(i in 1:nrow(d)){
        set[i]=max(d[i,nodes])
      }

      set.r <- d[ncol(d),]
      for(j in 1:ncol(d)){
        set.r[j]=max(d[nodes,j])
      }
    }
  }else{
    if (length(nodes)==1){
      set <- d[,nodes]
      set.r <- d[nodes,]
    }else{
      max <- max(d)+1
      d <- replace(d,d==0,max)
      set <- d[,ncol(d)]
      for(i in 1:nrow(d)){
        set[i]=min(d[i,nodes])
      }
      set.r <- d[ncol(d),]
      for(j in 1:ncol(d)){
        set.r[j]=min(d[nodes,j])
      }
      d <- replace(d,d==max,0)
      set <- replace(set,set==max,0)
      set.r <- replace(set.r,set.r==max,0)
    }
  }

  # update the adj.matrix
  set.r <- t(set.r)
  set.r <- cbind(set.r,0)
  d <- cbind(d,set)
  d <- rbind(d,set.r)
  rownames(d)=c(colnames(d))
  ColToDelete <- nodes
  d <- d[,-ColToDelete] # if v is names, use d <- d[,!(names(x) %in% ColToDelete)]
  RowToDelete <- t(nodes)
  d <- d[-RowToDelete,] # if v is names, use d <- d[!(names(x) %in% RowToDelete),]
  d
}


##################################################################################
### group.mreach.degree
#' Compute the group-level mreach.degree Centrality Score in a Netwrok
#'
#' \code{group.mreach.degree} computes the size of the reachable nodes from a particular node within M steps.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If there are multiple players,
#' use \code{c(index1,index2,...)}
#'
#' @param M Number indicating the maximum distance between two nodes,
#' above which the two nodes are considered disconnected.
#' m hence defines a reachable set. The default is \code{Inf}.
#'
#' @param method Indication of which grouping criterion should be used.
#' \code{method="min"} indicates the "minimum" criterion (edge values as distances).
#' \code{method="max"} indicates the "maximum" criterion (edge values as non-cummulative strengths).
#' \code{method="add"} indicates the "addition" criterion (edge values as cummulative strengths).
#' \code{method="union"} indicates the "union" criterion (edge values as probability).
#' The default is the "minimum" criterion for mreach.degree centrality.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The default is to report the total degree.
#' \code{"outdegree"} and \code{"indegree"} refer to indegree and outdegree
#' respectively. If \code{"all"}, all the three types are reported.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'
#' @return A vector indicating the outdegree, indegree, or total-degree
#' mreach.degree score of the chosen player(s); or a data frame containing all
#' the above information.
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @keywords internal

group.mreach.degree=function(adj.matrix, nodes, M = Inf, method = "min", binary = TRUE, cmode = "total", large = TRUE, geodist.precomp = NULL){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  rownames(d) <- c(colnames(d))

  if(NCOL(d[,nodes])==1){
    # create matrix of geodesic distances
    if (is.null(geodist.precomp)){
      if (isTRUE(large)){
        if (isTRUE(binary)){
          g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = NULL)
        }else{
          g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = TRUE)
        }
        distances=igraph::shortest.paths(g, mode = "out")
      }else{
        distances=sna::geodist(d,ignore.eval = binary)$gdist
      }
    }else{
      distances=geodist.precomp
    }
    colnames(distances) <- c(colnames(d))
    rownames(distances) <- c(colnames(d))

    # set step constraint using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }

    # calculate the score
    diag(distances)=Inf
    weights=1/distances
    weights[weights>0]=1
    kppos.out=apply(weights,1,sum) # column vector, outdegree
    kppos.in=apply(weights,2,sum) # row vector, indegree
    kppos.total=kppos.out + kppos.in
    kppos=c(kppos.out[nodes], kppos.in[nodes],kppos.total[nodes])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
    # reports
    if (missing(cmode)){
      kppos.total[nodes]
    }else if (cmode=="outdegree"){
      kppos.out[nodes]
    }else if (cmode=="indegree"){
      kppos.in[nodes]
    }else if (cmode=="all"){
      kppos
    }else{
      kppos.total[nodes]
    }
  }else{
    # geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    diag(distances)=Inf
    weights=1/distances
    d=contract(d,nodes,method)

    # create matrix of geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    colnames(distances) <- c(colnames(d))
    rownames(distances) <- c(colnames(d))

    # set step using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }

    # calculate the score
    diag(distances)=Inf
    weights=1/distances
    weights[weights>0]=1
    kppos.out=apply(weights,1,sum) # column vector, outdegree
    kppos.in=apply(weights,2,sum) # row vector, indegree
    kppos.total=kppos.out + kppos.in
    kppos=c(kppos.out["set"], kppos.in["set"], kppos.total["set"])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
    if (missing(cmode)){
      kppos.total["set"]
    }else if (cmode=="outdegree"){
      kppos.out["set"]
    }else if (cmode=="indegree"){
      kppos.in["set"]
    }else if (cmode=="all"){
      kppos
    }else{
      kppos.total["set"]
    }
  }
}

##################################################################################
### group.mreach.closeness
#' Compute the Group-level mreach.closeness Centrality Score in a Netwrok
#'
#' \code{mreach.closeness} refines the \code{\link{mreach.degree}} centrality by
#' using the (inverse) geodistance as weights.
#' The edge values should be properly interpreted as distances.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If there are multiple players,
#' use \code{c(index1,index2,...)}
#'
#' @param M Number indicating the maximum distance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}.
#'
#' @param method Indication of which grouping criterion should be used.
#' \code{method="min"} indicates the "minimum" criterion (edge values as distances).
#' \code{method="max"} indicates the "maximum" criterion (edge values as non-cummulative strengths).
#' \code{method="add"} indicates the "addition" criterion (edge values as cummulative strengths).
#' \code{method="union"} indicates the "union" criterion (edge values as probability).
#' The default is the "minimum" criterion for mrach.closeness centrality.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The default is to report the total degree.
#' \code{"outdegree"} and \code{"indegree"} refer to indegree and outdegree
#' respectively. If \code{"all"}, all the three types are reported.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'
#' @return A vector indicating the outdegree, indegree, or total-degree
#' cohesion score of the chosen players; or a data frame containing all
#' the above information. Note that the outdegree and indegree scores are normalized
#' to [0,1]. This means that the total-degree score is between [0,2].
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @keywords internal


group.mreach.closeness=function(adj.matrix, nodes, M = Inf, method = "min", binary = FALSE, cmode = "total", large = TRUE, geodist.precomp = NULL){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  rownames(d) <- c(colnames(d))

  if(NCOL(d[,nodes])==1){
    # create matrix of geodesic distances
    if (is.null(geodist.precomp)){
      if (isTRUE(large)){
        if (isTRUE(binary)){
          g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = NULL)
        }else{
          g <- igraph::graph.adjacency(adj.matrix, mode = "directed", weighted = TRUE)
        }
        distances=igraph::shortest.paths(g, mode = "out")
      }else{
        distances=sna::geodist(d,ignore.eval = binary)$gdist
      }
    }else{
      distances=geodist.precomp
    }
    colnames(distances) <- c(colnames(d))
    rownames(distances) <- c(colnames(d))

    # set threshold constraint using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }

    # calculate the score
    diag(distances)=Inf
    weights=1/distances
    sum.out=apply(weights,1,sum) # column vector, outdegree
    sum.in=apply(weights,2,sum) # row vector, indegree
    kppos.out = sum.out/((ncol(distances)-1)*max(weights)) # normalize to [0,1]
    kppos.in = sum.in/((ncol(distances)-1)*max(weights)) # normalize to [0,1]
    kppos.total=kppos.out + kppos.in
    kppos=c(kppos.out[nodes], kppos.in[nodes],kppos.total[nodes])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
    # reports
    if (missing(cmode)){
      kppos.total[nodes]
    }else if (cmode=="outdegree"){
      kppos.out[nodes]
    }else if (cmode=="indegree"){
      kppos.in[nodes]
    }else if (cmode=="all"){
      kppos
    }else{
      kppos.total[nodes]
    }
  }else{
    # geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    diag(distances)=Inf
    weights=1/distances
    m=max(weights)

    d=contract(d,nodes,method)

    # create matrix of geodesic distances
    distances=sna::geodist(d,ignore.eval=binary)$gdist
    colnames(distances) <- c(colnames(d))
    rownames(distances) <- c(colnames(d))

    # set threshold using distances matrix
    if(missing(M)){
      distances = distances
    }else{
      for(i in 1:ncol(distances)){
        for(j in 1:nrow(distances)){
          if(distances[i,j]>M) distances[i,j]=Inf
        }
      }
    }

    # calculate the score
    diag(distances)=Inf
    weights=1/distances
    sum.out=apply(weights,1,sum) # column vector, outdegree
    sum.in=apply(weights,2,sum) # row vector, indegree
    kppos.out = sum.out/((ncol(distances)-1)*m) # normalize to [0,1]
    kppos.in = sum.in/((ncol(distances)-1)*m) # normalize to [0,1]
    kppos.total=kppos.out + kppos.in
    kppos=c(kppos.out["set"], kppos.in["set"], kppos.total["set"])
    names(kppos) <- c("outdegree","indegree","total")
    kppos=as.data.frame(t(kppos))
    rownames(kppos) <- "score"
    if (missing(cmode)){
      kppos.total["set"]
    }else if (cmode=="outdegree"){
      kppos.out["set"]
    }else if (cmode=="indegree"){
      kppos.in["set"]
    }else if (cmode=="all"){
      kppos
    }else{
      kppos.total["set"]
    }
  }
}


###################################################################################
### kpcent

#' Compute Group Centraltiy in a Network
#'
#' \code{kpcent} reports the group-level centrality scores.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network or the probability matrix in the case of calculating diffusion centrality.
#'
#' @param nodes Integer indicating the column index of the chosen player
#' in the adjacenncy matrix. If there are multiple players,
#' use \code{c(index1,index2,...)}
#'
#' @param type
#' \code{type="betweenness"} for \code{\link[sna]{betweenness}} centrality. \cr
#' \code{type="closeness"} for \code{\link[sna]{closeness}} centrality. \cr
#' \code{type="degree"} for \code{\link[sna]{degree}} centraslity. \cr
#' \code{type="diffusion"} for \code{\link{diffusion}} centrality. \cr
#' \code{type="evcent"} for \code{\link[sna]{evcent}} (eigenvector) centrality. \cr
#' \code{type="fragment"} for \code{\link{fragment}} centrality. \cr
#' \code{type="mreach.degree"} for \code{\link{mreach.degree}} centrality. \cr
#' \code{type="mreach.closeness"} for \code{\link{mreach.closeness}} centrality. \cr
#'
#'
#' @param method Indication of which grouping criterion should be used. \cr
#' \code{"min"} indicates the "minimum" criterion and is the default for
#' betweenness, closeness, fragmentation, and M-reach centralities. \cr
#' \code{"max"} indicates the "maximum" criterion and is the default for
#' degree and eigenvector centralities.\cr
#' \code{"add"} indicates the "addition" criterion.\cr
#' \code{"union"} indicates the "union" criterion and is the default for
#' diffusion centrality.\cr
#' See Details section for explanations on grouping method.
#'
#' @param M Positive number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}. The option is applicable to mreach.degree, mreach.closeness,
#' and fragmentation centralities.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. For diffusion centrality only. By default, T is the network size.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. By default, \code{binary=FALSE}
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The option is applicable to degree and M-reach centralities.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures. The default is to report the total degree.
#' Note for closeness centrality, we use the Gil-Schmidt power index when \code{large=FALSE}.
#' See \code{\link[sna]{closeness}} for explanation. When
#' large=\code{TRUE}, the function reports the standard closeness score.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'
#' @return A vector indicating the centrality score of a group.
#'
#' @details The basic idea of measuring the group-level centrality is to treat
#' a group of nodes as a large pseudo-node. We propose several methods to measure
#' the tie status between this pseudo node and other nodes, responding to several
#' common edge value interpretations (An and Liu, 2015).
#'
#' Minimum Criterion: the edge value between a group and an outside node
#' is measured as the minimal value among all the (nonzero) edge values between
#' any node in the group and the outside node. Suggested if edge values are
#' interpreted as distances.\cr
#' \emph{Example: suppose node A to C has distance 2 and B to C has distance 1,
#' then according to the minimum criterion, the distance between C and
#' the merged set AB is 1. Note that if B and C are not connected,
#' the algorithm takes the distance between A and C to describe
#' the distance between AB and C.}
#'
#' Maximun Criterion: the edge value between a group and an outside node
#' is measured as the maximal value among all the (nonzero) edge values between
#' any node in the group and the outside node. Suggested if edge values are
#' interpreted as non-cummulative strengths. \cr
#' \emph{Example: we keep using the above example, but the figure now indicates
#' the strength of tie. According to the maximum criterion, the strength of tie
#' between AB and C is 2.}
#'
#' Addition Criterion: the edge value between a group and an outside node
#' is measured as the sum of all the edge values between any node in the group
#' and the outside node. Suggested if edge values are as cummulative strengths. \cr
#' \emph{Example: according to the addition criterion, the strength of tie between
#' AB and C is 3}
#'
#' Union Criterion: the edge value between a group and an outside node is measured
#' as the probability that there is at least one path connecting the group with
#' the outside node. Suggested if edge values are as probability. \cr
#' \emph{Example: suppose A has probability 0.2 to reach C and B has probability
#' 0.5 to reach C, then C can be reached from merged AB with probability
#' 1-(1-0.2)*(1-0.5)=0.6 according to the union criterion.}
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua. (2015). "Multilevel Meta Network Analysis with Application to Studying Network Dynamics of Network Interventions." \emph{Social Networks} 43: 48-56.\cr
#'
#' An, Weihua and Yu-Hsin Liu (2016). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' Working Paper, Indiana Univeristy.\cr
#'
#' Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2013):
#' "Diffusion of Microfinance," \emph{Science}, Vol. 341. p.363\cr
#'
#' Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2014):
#' "Gossip: Identifying Central Individuals in a Social Network,"
#' Working Paper.\cr
#'
#' Borgatti, Stephen P. (2006). "Identifying Sets of Key Players in a Network."
#' \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr
#'
#' Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
#' version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr
#'
#' Csardi, G and Nepusz, T (2006). "The igraph software package for complex network research."
#' InterJournal, Complex Systems 1695. \url{http://igraph.org} \cr
#'
#' @seealso
#' \code{\link{contract}}
#' \code{\link{kpset}}
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix,
#' # where edge values represent the strength of tie
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # List the degree centrality for group of node 2 and 3
#' kpcent(W,c(2,3),type="degree")
#'
#' # Transform the edge value to distance interpretaion
#' # Compute the fragmentation centrality for node 2
#' A <- W
#' A[W!=0] <- 1/W[W!=0]
#' kpcent(A,2,type="fragment")
#'
#' # Replicate the group-level degree centrality (normalized) when the weights
#' # are given by the inverse distances and report the outgoing score only
#' kpcent(A,c(2,3),type="mreach.closeness",binary=TRUE,M=1,cmode="outdegree")
#'
#' # Transform the edge value to probability interpretation
#' # Compute the diffusion centrality with number of iteration 20
#' P <- 0.1*W
#' kpcent(P,c(2,3),type="diffusion",T=20)
#'
#' @export

kpcent=function(adj.matrix, nodes, type, M=Inf, T=ncol(adj.matrix), method, binary=FALSE, cmode, large = TRUE, geodist.precomp = NULL){
  d <- data.matrix(adj.matrix, rownames.force = NA)

  if(type == "betweenness"){
    if (length(nodes)==1){
      if(isTRUE(large)){
        if (isTRUE(binary)){
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = NULL)
        }else{
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = TRUE)
        }
        betweenness = igraph::betweenness(g, v = nodes)
      }else{
        betweenness = sna::betweenness(d, nodes = nodes, ignore.eval = binary, rescale = FALSE, geodist.precomp = geodist.precomp)
      }
    }else{
      if(missing(method)){
        d <- contract(d, nodes, method="min")
      }else{
        d <- contract(d, nodes, method)
      }
      s = ncol(d)
      if(isTRUE(large)){
        if (isTRUE(binary)){
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = NULL)
        }else{
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = TRUE)
        }
        betweenness = igraph::betweenness(g)[s]
      }else{
       betweenness <- sna::betweenness(d, nodes = s, ignore.eval = binary, rescale = FALSE)
      }
    }
    score <- betweenness

  }else if(type == "closeness"){
    if(isTRUE(large)){
      if(length(nodes)==1){
        if(isTRUE(binary)){
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = NULL)
        }else{
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = TRUE)
        }
        closeness = igraph::closeness(g, v = nodes, mode = "out", normalized = TRUE)
      }else{
        if(missing(method)){
          d <- contract(d, nodes, method="min")
        }else{
          d <- contract(d, nodes, method)
        }
        s = ncol(d)

        if(isTRUE(binary)){
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = NULL)
        }else{
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = TRUE)
        }
        closeness = igraph::closeness(g, v = s, mode = "out", normalized = TRUE)
      }
    }else{
        if (length(nodes)==1){
          closeness = sna::closeness(d, nodes = nodes, cmode = "suminvdir", ignore.eval = binary)
        }else{
          if(missing(method)){
            d <- contract(d, nodes, method="min")
          }else{
            d <- contract(d, nodes, method)
          }
          s = ncol(d)
          closeness <- sna::closeness(d, nodes = s, cmode = "suminvdir", ignore.eval = binary)
        }
    }
    score <- closeness

  }else if(type == "degree"){
    if(length(nodes)==1){
      indegree = sna::degree(d, nodes = nodes, cmode = "indegree", ignore.eval = binary)
      outdegree = sna::degree(d, nodes = nodes, cmode = "outdegree", ignore.eval = binary)
      total = sna::degree(d, nodes = nodes, cmode = "freeman", ignore.eval = binary)
    }else{
      if (missing(method)){
        d <- contract(d,nodes, method="max")
      }else{
        d <- contract(d, nodes, method)
      }
      s = ncol(d)
      indegree = sna::degree(d, nodes = s, cmode = "indegree", ignore.eval = binary)
      outdegree = sna::degree(d, nodes = s, cmode = "outdegree", ignore.eval = binary)
      total = sna::degree(d, nodes = s, cmode = "freeman", ignore.eval = binary)
    }

    set.degree=c(outdegree, indegree, total)
    names(set.degree) <- c("outdegree","indegree","total")
    set.degree=as.data.frame(t(set.degree))
    rownames(set.degree) <- "score"

    if (missing(cmode)){
      score <- total
    }else if (cmode=="outdegree"){
      score <- outdegree
    }else if (cmode=="indegree"){
      score <- indegree
    }else if (cmode=="all"){
      score <- set.degree
    }else{
      score <- total
    }
  }else if(type == "diffusion"){
    s <- 0
    if(NCOL(d[,nodes])==1){
      for (t in 1:T){
        s=s+matpow::matpow(d,t)$prod1
      }
      diffusion=apply(s,1,sum)
      score <- diffusion[nodes]
    }else{
      if(missing(method)){
        d <- contract(d, nodes, method="union")
      }else{
        d <- contract(d, nodes, method)
      }
      s=0
      for (t in 1:T){
        s=s+matpow::matpow(d,t)$prod1
      }
      diffusion=apply(s,1,sum)
      score <- diffusion["set"]
    }
  }else if(type == "evcent"){
    if(length(nodes)==1){
      if(isTRUE(large)){
        if (isTRUE(binary)){
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = NULL)
        }else{
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = TRUE)
        }
        evcent = igraph::evcent(g, directed = TRUE, scale = FALSE)$vector[nodes]
      }else{
        evcent = sna::evcent(d, nodes = nodes, ignore.eval = binary)
      }
    }else{
      if(missing(method)){
        d <- contract(d, nodes, method="max")
      }else{
        d <- contract(d, nodes, method)
      }
      s = ncol(d)
      if(isTRUE(large)){
        if (isTRUE(binary)){
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = NULL)
        }else{
          g <- igraph::graph.adjacency(d, mode = "directed", weighted = TRUE)
        }
        evcent = igraph::evcent(g, directed = TRUE, scale = FALSE)$vector[s]
      }else{
        evcent = sna::evcent(d, nodes = s, ignore.eval = binary)
      }
    }
    score <- evcent
  }else if(type == "fragment"){
    score <- fragment(adj.matrix, nodes, M, binary, large, geodist.precomp)
  }else if(type == "mreach.degree"){
    if(missing(method)){
      score <- group.mreach.degree(adj.matrix, nodes, M, method="min", binary, cmode, large, geodist.precomp)
    }else{
      score <- group.mreach.degree(adj.matrix, nodes, M, method, binary, cmode, large, geodist.precomp)
    }
  }else if(type == "mreach.closeness"){
    if(missing(method)){
      score <- group.mreach.closeness(adj.matrix, nodes, M, method="min", binary, cmode, large, geodist.precomp)
    }else{
      score <- group.mreach.closeness(adj.matrix, nodes, M, method, binary, cmode, large, geodist.precomp)
    }
  }
  names(score) <- NULL
  return(score)
}

####################################################################################
### kpset
### Step 1. create delta.score to compute the change from switching i to j

#' The change on objective function for greedy search implimentation
#'
#' \code{delta.score} calculates the change in group centrality score.
#' @param adj.matrix The adjacency matrix of the network.
#' @param candidate An initial set of players selected.
#' @param residual The remaining players in the network without members in the initial set.
#' @param i The specific member in the candidate set to be replaced
#' @param j The specific member in the residual set to replace \code{i} in the candidate set
#' @param type Choose
#' \code{type="betweenness"} for betweenness centrality,
#' \code{type="closeness"} for closeness centrality,
#' \code{type="degree"} for degree centraslity,
#' \code{type="diffusion"} for diffusion centrality.
#' \code{type="evcent"} for eigenvector centrality,
#' \code{type="fragment"} for fragmentation centrality,
#' \code{type="mreach.degree"} for mreach.degree centrality, and
#' \code{type="mreach.closeness"} for mreach.closeness centrality.
#'
#' @param M Positive number indicating the maximum distance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}.
#' The option is applicable to mreach.degree, mreach.closeness, fragmentation,
#' and diffusion centralities.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. For diffusion centrality only.
#' In the first iteration, the adjacency matrix
#' is as the input. In the nth iteration, the adjacency matrix becomes
#' the input adjacency matrix to the power of n. By default, T is the network size.
#'
#' @param method Indication of which grouping criterion should be used.
#' \code{method="min"} indicates the "minimum" criterion (edge values as distances).
#' \code{method="max"} indicates the "maximum" criterion (edge values as non-cummulative strengths).
#' \code{method="add"} indicates the "addition" criterion (edge values as cummulative strengths).
#' \code{method="union"} indicates the "union" criterion (edge values as probability).
#' By default, the minimun criterion is used for betweenness, closeness, fragmentation,
#' mreach.degree, and mreach.closeness centralities.
#' The maximun criterion is used for degree and eigenvector centralities.
#' The union criterion is used for diffusion centrality.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. By default, \code{binary=FALSE}
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The option is applicable to degree, mreach.degree, and mreach.closeness centralities.
#' The default is to report the total degree.
#' \code{cmode="outdegree"} and \code{cmode="indegree"} refer to indegree and outdegree
#' respectively. If \code{cmode="all"}, all the three types are reported.
#' The option can also applicable to closeness centrality.
#' See \code{\link[sna]{closeness}} Details section.
#' The default is to use the Gil-Schmidt power index.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'
#' @return The change in group centrality score
#' @keywords internal

delta.score=function(adj.matrix,candidate,residual,i,j,type,M=Inf,T=ncol(adj.matrix),method,binary=FALSE,cmode,large=TRUE,geodist.precomp=NULL){

  s <- kpcent(adj.matrix,candidate,type,M,T,method,binary,cmode,large,geodist.precomp=NULL)
  candidate[i] <- residual[j]
  delta.s = kpcent(adj.matrix,candidate,type,M,T,method,binary,cmode,large,geodist.precomp=NULL)-s
  delta.s
}

################################################################################
#' Return the most central player in sequentially reduced networks
#'
#' \code{kpnode} returns the node with the highest centrality score.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param type
#' \code{type="betweenness"} for \code{\link[sna]{betweenness}} centrality. \cr
#' \code{type="closeness"} for \code{\link[sna]{closeness}} centrality. \cr
#' \code{type="degree"} for \code{\link[sna]{degree}} centraslity. \cr
#' \code{type="diffusion"} for \code{\link{diffusion}} centrality. \cr
#' \code{type="evcent"} for \code{\link[sna]{evcent}} (eigenvector) centrality. \cr
#' \code{type="fragment"} for \code{\link{fragment}} centrality. \cr
#' \code{type="mreach.degree"} for \code{\link{mreach.degree}} centrality. \cr
#' \code{type="mreach.closeness"} for \code{\link{mreach.closeness}} centrality. \cr
#'
#'
#' @param M Positive number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}. The option is applicable to mreach.degree, mreach.closeness,
#' and fragmentation centralities.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. For diffusion centrality only.
#' In the first iteration, the adjacency matrix
#' is as the input. In the nth iteration, the adjacency matrix becomes
#' the input adjacency matrix to the power of n. By default, T is the network size.
#'
#' @param method Indication of which grouping criterion should be used. \cr
#' \code{"min"} indicates the "minimum" criterion and is suggested for
#' betweenness, closeness, fragmentation, and M-reach centralities. \cr
#' \code{"max"} indicates the "maximum" criterion and is suggested for
#' degree and eigenvector centralities.\cr
#' \code{"add"} indicates the "addition" criterion and is suggested for
#' degree and eigenvector centralities as an altenative of "max".\cr
#' \code{"union"} indicates the "union" criterion and is suggested for
#' diffusion centrality.\cr
#' By default, kpset uses "min".
#' See \code{\link{kpcent}} Details section for explanations on grouping method.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. By default, \code{binary=FALSE}
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The option is applicable to degree and M-reach centralities.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures. The default is to report the total degree.
#' The option also applies to closeness centrality, but with different options.
#' The default is to use the Gil-Schmidt power index as the closeness measure.
#' See \code{\link[sna]{closeness}} for complete options.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'
#'
#' @return The most central player and its centrality score
#' @keywords internal

kpnode=function(adj.matrix, type, M=Inf, T=ncol(adj.matrix),method, binary=FALSE, cmode, large=TRUE, geodist.precomp=NULL){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  N <- c(1:NCOL(d))
  C <- sample(N,1)
  R <- N[!N %in% C]
  for (j in 1:NROW(R)){
    if(delta.score(adj.matrix,C,R,1,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) C[1]<-R[j]
  }
  call <- match.call()
  col1 <- C
  col1 <- sort(col1)
  col2 <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
  names(col2) <- NULL
  result <- list(keyplayers=col1, centrality=col2)
  result
}

################################################################################
#' Orders players by individual centrality from high to low
#'
#' \code{topnode} Orders players by individual centrality from high to low
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network.
#'
#' @param type
#' \code{type="betweenness"} for \code{\link[sna]{betweenness}} centrality. \cr
#' \code{type="closeness"} for \code{\link[sna]{closeness}} centrality. \cr
#' \code{type="degree"} for \code{\link[sna]{degree}} centraslity. \cr
#' \code{type="diffusion"} for \code{\link{diffusion}} centrality. \cr
#' \code{type="evcent"} for \code{\link[sna]{evcent}} (eigenvector) centrality. \cr
#' \code{type="fragment"} for \code{\link{fragment}} centrality. \cr
#' \code{type="mreach.degree"} for \code{\link{mreach.degree}} centrality. \cr
#' \code{type="mreach.closeness"} for \code{\link{mreach.closeness}} centrality. \cr
#'
#'
#' @param M Positive number indicating the maximum geodistance between two nodes,
#' above witch the two nodes are considered disconnected. The default is
#' \code{Inf}. The option is applicable to mreach.degree, mreach.closeness,
#' and fragmentation centralities.
#'
#' @param T Integer indicating the maximum number of iterations
#' of communication process. For diffusion centrality only.
#' In the first iteration, the adjacency matrix
#' is as the input. In the nth iteration, the adjacency matrix becomes
#' the input adjacency matrix to the power of n. By default, T is the network size.
#'
#' @param method Indication of which grouping criterion should be used. \cr
#' \code{"min"} indicates the "minimum" criterion and is suggested for
#' betweenness, closeness, fragmentation, and M-reach centralities. \cr
#' \code{"max"} indicates the "maximum" criterion and is suggested for
#' degree and eigenvector centralities.\cr
#' \code{"add"} indicates the "addition" criterion and is suggested for
#' degree and eigenvector centralities as an altenative of "max".\cr
#' \code{"union"} indicates the "union" criterion and is suggested for
#' diffusion centrality.\cr
#' By default, kpset uses "min".
#' See \code{\link{kpcent}} Details section for explanations on grouping method.
#'
#' @param binary If \code{TRUE}, the adjacency matrix is binarized.
#' If \code{FALSE}, the edge values are considered. By default, \code{binary=FALSE}
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The option is applicable to degree and M-reach centralities.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and (total) degree respectively. \code{"all"} reports
#' all the above measures. The default is to report the total degree.
#' The option also applies to closeness centrality, but with different options.
#' The default is to use the Gil-Schmidt power index as the closeness measure.
#' See \code{\link[sna]{closeness}} for complete options.
#'
#' @param large Logical scalar, whether the computation method for large network is
#' implemented. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used; otherwise the method implemented in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the graph to be
#' analyzed (optional).
#'#'
#' @return A vector of indices by individual centrality from high to low
#' @keywords internal

topnode=function(adj.matrix, type, M=Inf, T=ncol(adj.matrix),method, binary=FALSE, cmode, large=TRUE, geodist.precomp=NULL){
  d <- data.matrix(adj.matrix, rownames.force = NA)
  N <- NCOL(d)
  id <- seq(1, N, 1)
  topd <- NULL
  for (i in 1:N) {
    topd[i] <- kpcent(adj.matrix,i,type,M,T,method,binary,cmode,large,geodist.precomp)
  }
  ID <- id[order(topd, decreasing = TRUE)]
  ID
}

################################################################################
### Step 2. keyplayers greedy search algorithm

#' Selecting the Most Central Group of Players in a Network
#'
#' \code{kpset} helps identify the most central group of players in a social network given a sepcified centraliy measure and a target group size.
#'
#' @param adj.matrix Matrix indicating the adjacency matrix of the network or in the case of diffusion centrality a probability matrix.
#'
#' @param size Integer indicating the target size of players.
#'
#' @param type A string indicating the type of centrality measure to be used. Should be one of \code{"degree"} for degree centrality,
#' \code{"closeness"} for closeness centrality, \code{"betweenness"} for betweenness centrality, \code{"evcent"} for eigenvector centrality,
#' \code{"mreach.degree"} for M-reach degree centrality, \code{"mreach.closeness"} for M-reach closeness centrality,
#' \code{"fragment"} for fragment centrality, and \code{"diffusion"} for diffusion centrality.
#'
#' @param method Indication of which grouping criterion should be used. \cr
#' \code{"min"} indicates the "minimum" criterion and is suggested for
#' betweenness, closeness, fragmentation, and M-reach centralities. \cr
#' \code{"max"} indicates the "maximum" criterion and is suggested for
#' degree and eigenvector centralities.\cr
#' \code{"add"} indicates the "addition" criterion and is suggested for
#' degree and eigenvector centralities as an altenative of "max".\cr
#' \code{"union"} indicates the "union" criterion and is suggested for
#' diffusion centrality.\cr
#' The default is "min".
#' See \code{\link{kpcent}} Details section for explanations on grouping method.
#'
#' @param M Positive number indicating the maximum geodistance between two nodes,
#' above which the two nodes are considered disconnected. The default is \code{Inf}.
#' The option is applicable to M-reach degree, M-reach closeness, and fragmentation centralities..
#'
#' @param T Integer indicating the maximum number of iterations
#' in the communication process. By default, T is the network size.
#'
#' @param binary If \code{TRUE}, the input matrix is binarized.
#' If \code{FALSE}, the edge values are considered. The default is \code{FALSE}.
#'
#' @param cmode String indicating the type of centrality being evaluated.
#' The option is applicable to degree and M-reach centralities.
#' \code{"outdegree"}, \code{"indegree"}, and \code{"total"} refer to
#' indegree, outdegree, and total degree, respectively. \code{"all"} reports
#' all the above measures. The default is to report the total degree.
#' Note for closeness centrality, we use the Gil-Schmidt power index when \code{large=FALSE}.
#' See \code{\link[sna]{closeness}} for explanation. When
#' large=\code{TRUE}, the function reports the standard closeness score.
#'
#' @param large Logical scalar. If \code{TRUE} (the default), the method implmented in \pkg{igraph} is
#' used for computing geodistance and related centrality measures; otherwise the method in \pkg{sna} is used.
#'
#' @param geodist.precomp Geodistance precomputed for the network to be
#' analyzed (optional).
#'
#' @param seed String indicating the seeding method or a vector of the seeds specified by user.
#' If \code{"top"}, players with the high individual centrality
#' are used as the seeds. If \code{"random"}, seeds are randomly sampled.
#' The default is \code{"top"} for efficiency.
#'
#' @param parallel Logical scalar. IF \code{TRUE}, the parallel computation is
#' implement. The default is \code{FALSE}.
#'
#' @param cluster Integer indicating the number of CPU cores to be used for parallel computation.
#'
#' @param round Integer indicating the "length" of search,
#' namely, the number of loops over the nodes in the candidate set.
#'
#' @param iteration Integer indicating the "width" of search in each round,
#' namely, the number of loops over the nodes in the residual set.
#'
#' @return \code{kpset} returns the column indices of the players who form
#' the most central set and its centrality score.
#'
#' @details
#' The most central group of players in a network is not necessarily the set of players who are the most
#' central as individuals because there may be redundancy in their connections.
#' Currenlty a greedy search algorithm is implemented in this package to identify the most central group of key players. The basic steps are shown as follows.
#' \enumerate{
#' \item Select an initial candidate set \emph{C}. The residual set is denoted as \emph{R}.
#' \item Update the candidate set \emph{C}.
#' \itemize{
#'  \item Start with the first node in \emph{C}. Try to swap it with nodes in \emph{R} sequentially (loop 1).
#'  Make the swap if it improves the centrality score of the resulting \emph{C}.
#'  The number of loop 1 is defined as the number of iterations (over the nodes in the residual set).
#'  \item Repeat step 1 for each node in \emph{C} sequentially (loop 2). The number of loop 2 is
#'   defined as the number of rounds (over the nodes in the candidate set).
#'   \item Stop if (a) the change in \emph{C}'s centrality score is negligible
#'         (i.e. it is smaller than a pre-specified threshold determined by both
#'         the network size and edge values.)
#'         or (b) the process reaches a specified number
#'         of rounds.
#'   }
#' \item Return the final candidate set and the centrality score.
#' }
#'
#' It is recommended to run \code{kpset} several times with different seeds so that the algorithm will not be trapped in a local optimum.
#' To facilitate the search in large networks, users may specify a reasonable number of iterations
#' or rounds and/or utilize parallel computation. During parallel computation, for each cluster and each iteration
#' the algorithm randomly picks a node from the candidate set and the residual set, respectively,
#' and swaps the two if it improves the centrality score of the candidate set. It repeats this process until exhausting the specified iterations and rounds
#' and then compare and combine the results from the clusters.
#'
#' @author Weihua An \email{weihuaan@@indiana.edu}; Yu-Hsin Liu \email{yuhsliu@@indiana.edu}
#'
#' @references
#' An, Weihua. (2015). "Multilevel Meta Network Analysis with Application to Studying Network Dynamics of Network Interventions." \emph{Social Networks} 43: 48-56.\cr
#'
#' An, Weihua and Yu-Hsin Liu (2016). "keyplayer: An R Package for Locating Key Players in Social Networks."
#' Working Paper, Indiana Univeristy.\cr
#'
#' Borgatti, Stephen P. (2006). "Identifying Sets of Key Players in a Network." \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr
#'
#' Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
#' version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr
#'
#' Csardi, G and Nepusz, T (2006). "The igraph software package for complex network research."
#' InterJournal, Complex Systems 1695. \url{http://igraph.org} \cr
#'
#' @seealso
#' \code{\link{kpcent}}
#'
#' @examples
#' # Create a 5x5 weighted and directed adjacency matrix
#' W <- matrix(
#'   c(0,1,3,0,0,
#'     0,0,0,4,0,
#'     1,1,0,2,0,
#'     0,0,0,0,3,
#'     0,2,0,0,0),
#'     nrow=5, ncol=5, byrow = TRUE)
#'
#' # Find the most central player set sized 2 in terms of the degree centrality
#' kpset(W,size=2,type="degree")
#'
#' # Find two most central players in terms of indegree
#' # via parallel computation using 5 cpu cores
#' kpset(W,size=2,type="degree", cmode="indegree", parallel = TRUE, cluster = 2)
#'
#' @export

kpset=function(adj.matrix, size, type = "degree", M=Inf, T=ncol(adj.matrix), method="min", binary=FALSE, cmode="total", large=TRUE, geodist.precomp=NULL, seed="top", parallel=FALSE, cluster=2, round=10, iteration=ncol(adj.matrix)){

  # general set-up
  d <- data.matrix(adj.matrix, rownames.force = NA)
  N <- c(1:NCOL(d))
  d1 = d
  d1[d1==0] = Inf
  m <- max(max(d),1/min(d1))
  storage <- c(0:round)

  if (is.numeric(seed)){
    if(isTRUE(parallel)){
      # designated seed
      C <- seed
      R <- N[!N %in% C]
      # round: parallel computations in each round, choose the best of each cluster
      for (r in 2:(round+1)){
        par <- function(...){
          for (t in 1:iteration){
            i <- sample.int(NROW(C),1)
            j <- sample.int(NROW(R),1)
            if(delta.score(adj.matrix,C,R,i,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) {C[i]<-R[j]}
            R <- N[!N %in% C]
          }
          col1 <- C
          col1 <- sort(col1)
          col2 <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
          names(col2) <- NULL
          result1 <- list(keyplayers=col1,centrality=col2)
          result1
        }
        df <- do.call(rbind, parallel::mclapply(seq_len(cluster), par))
        rank <- rank(as.numeric(df[,2]),ties.method = "first")
        df <- cbind(df,rank)
        for(i in 1:NROW(df)){
          if (df[i,3]==cluster) C<-df[i,1]
        }
        C <- C[[1]]
        R <- N[!N %in% C]
        storage[r] <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
      }
      if (storage[round+1]-storage[round]>(1/(2*(NCOL(d)^2)*(m^3)))) warning('The specified number of rounds have been reached. The solution may not converge yet.')
    }else{
      # designated seed
      C <- seed
      R <- N[!N %in% C]
      for (r in 2:(round+1)){
        for (i in 1:NROW(C)){
          for (j in 1:min(NROW(R),iteration)){
            if(delta.score(adj.matrix,C,R,i,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) {C[i] <- R[j]}
          }
          R <- N[!N %in% C]
        }
        storage[r] <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
      }
      if (storage[round+1]-storage[round]>(1/(2*(NCOL(d)^2)*(m^3)))) warning('The specified number of round has been reached. The solution may not converge yet.')
    }
  }else if (seed=="top"){
    C <- sample(N,size)
    R <- N[!N %in% C]
    if(NROW(C)==1){
      for (j in 1:NROW(R)){
        if(delta.score(adj.matrix,C,R,1,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) C[1]<-R[j]
      }
    }else{
      # "good" seed
      # Conventional method uses the players with the highest individual centrality as seeds. Likely a local trap.
      C <- topnode(adj.matrix,type,M,T,method,binary,cmode,large,geodist.precomp)[1:size]

#       # Similar to the KI method, but does not remove the friends of the selected.
#       R.temp <- d
#       for (i in 1:NROW(C)){
#         C[i] <- kpnode(R.temp,type,M,T,method,binary,cmode,large,geodist.precomp)$keyplayers
#         R.temp <- d[-C[1:i],-C[1:i]]
#       }

      R <- N[!N %in% C]
      if(isTRUE(parallel)){
        # round: parallel computations in each round, choose the best of each cluster
        for (r in 2:(round+1)){
          par <- function(...){
            for (t in 1:iteration){
              i <- sample.int(NROW(C),1)
              j <- sample.int(NROW(R),1)
              if(delta.score(adj.matrix,C,R,i,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) C[i]<-R[j]
              R <- N[!N %in% C]
            }
            col1 <- C
            col1 <- sort(col1)
            col2 <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
            names(col2) <- NULL
            result1 <- list(keyplayers=col1,centrality=col2)
            result1
          }
          df <- do.call(rbind, parallel::mclapply(seq_len(cluster), par))
          rank <- rank(as.numeric(df[,2]),ties.method = "first")
          df <- cbind(df,rank)
          for(i in 1:NROW(df)){
            if (df[i,3]==cluster) C<-df[i,1]
          }
          C <- C[[1]]
          R <- N[!N %in% C]
          storage[r] <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
        }
        if (storage[round+1]-storage[round]>(1/(2*(NCOL(d)^2)*(m^3)))) warning('The specified number of round has been reached. The solution may not converge yet.')
      }else{
        for (r in 2:(round+1)){
          for (i in 1:NROW(C)){
            for (j in 1:min(NROW(R),iteration)){
              if(delta.score(adj.matrix,C,R,i,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) C[i]<-R[j]
            }
            R <- N[!N %in% C]
          }
          storage[r] <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
        }
        if (storage[round+1]-storage[round]>(1/(2*(NCOL(d)^2)*(m^3)))) warning('The specified number of round has been reached. The solution may not converge yet.')
      }
    }
  }else if (seed=="random"){
    C <- sample(N,size)
    R <- N[!N %in% C]
    if(NROW(C)==1){
      for (j in 1:NROW(R)){
        if(delta.score(adj.matrix,C,R,1,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) C[1]<-R[j]
      }
    }else{
      if(isTRUE(parallel)){
        # random seed parallel
        par <- function(...){
          C <- sample(N,size)
          R <- N[!N %in% C]
          storage <- c(0:round)
          for (r in 2:(round+1)){
            for (i in 1:NROW(C)){
              R.sub <- sample(R,min(NROW(R),iteration)) # allow random order and subset
              for (j in 1:NROW(R.sub)){
                if(delta.score(adj.matrix,C,R,i,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) C[i]<-R.sub[j]
              }
              R <- N[!N %in% C]
            }
            storage[r] <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
          }
          if (storage[round+1]-storage[round]>(1/(2*(NCOL(d)^2)*(m^3)))) warning('The specified number of round has been reached. The solution may not converge yet.')

          col1 <- C
          col1 <- sort(col1)
          col2 <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
          names(col2) <- NULL
          result1 <- list(keyplayers=col1,centrality=col2)
          result1
        }
        df <- do.call(rbind, parallel::mclapply(seq_len(cluster), par))
        rank <- rank(as.numeric(df[,2]),ties.method = "first")
        df <- cbind(df,rank)
        for(i in 1:NROW(df)){
          if (df[i,3]==cluster) C<-df[i,1]
        }
        C <- C[[1]]

      }else{
        for (r in 2:(round+1)){
          for (i in 1:NROW(C)){
            for (j in 1:min(NROW(R),iteration)){
              if(delta.score(adj.matrix,C,R,i,j,type,M,T,method,binary,cmode,large,geodist.precomp)>0) C[i]<-R[j]
            }
            R <- N[!N %in% C]
          }
          storage[r] <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)

        }
        if (storage[round+1]-storage[round]>(1/(2*(NCOL(d)^2)*(m^3)))) warning('The specified number of round has been reached. The solution may not converge yet.')
      }
    }
  }
  call <- match.call()
  col1 <- C
  col1 <- sort(col1)
  col2 <- kpcent(adj.matrix,C,type,M,T,method,binary,cmode,large,geodist.precomp)
  names(col2) <- NULL
  result <- list(keyplayers=col1, centrality=col2)
  result
}


### Data: Friends
#' The friendship network of 21 managers in a high-tech company
#'
#' @usage data("Friends")
#'
#' @format A network object for the friendship network of 21 managers in a high-tech company (Krackhardt, 1987).
#' It is a directed network including 21 nodes and 60 edges. There are two node atrributes. "Dept" shows the department affiliations (1-4).
#' "Level" shows the rank of the managers in the company (1-3).
#'
#'  @source D. Krackhardt. "Cognitive social structures." Social Networks, 9:109-134, 1987.
"Friends"



