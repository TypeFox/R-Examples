rNetwork <- function(p, pi, alpha = c(1), directed = FALSE, name = "a network", signed=TRUE) {
  
  ## SOME INITIALIZATIONS
  A  <- matrix(0,p,p)
  Q  <- length(alpha)
  pi <- as.matrix(pi)
  fixed.edges <- FALSE
  
  ## CHECK if Pi represents either the connectivity or the number of edges by class
  if (max(pi) > 1) {
    if (sum(pi) %% floor(sum(pi)) == 0) {
      fixed.edges <- TRUE
    } else {      
      stop("\nPi does not sum to one nor to an integer...")
    }
  }
  if (sum(alpha) < 1) {
    cat("\nNormalizing the vector of prior proportions which does not sum up to one\n")
    alpha <- alpha/sum(alpha)
  }
  if (!directed & !isSymmetric(pi)) {
    stop("\nPi should be symmetric for undirected graph...")
  }
  
  ## CHECK PARAMETERS CONSISTENCY
  if (ncol(pi) != Q) {
    stop("\nPi and alpha do not have consistent dimensions...\n") 
  }
  if (fixed.edges) {
    if (!directed) {
      pi.bounds <- floor(1/2*(p*alpha) %*% t(p*alpha) - diag(p*alpha/2, nrow=Q))
    } else {
      pi.bounds <- floor((p*alpha) %*% t(p*alpha))
    }
    if (any(pi > pi.bounds)) {
      cat("\nMatrix Pi is too dense: here is the adjusted version\n")
      pi[pi > pi.bounds] <- pi.bounds[pi > pi.bounds]
      print(pi)
    }    
  }
  
  ## CLASSES' LEVELS AND VECTOR OF CONNECTIVITY ACCORDING TO THE STRUCTURE
  cl.levels <- factor(1:ncol(pi))
  if (is.null(rownames(pi))) {
    colnames(pi) <- rownames(pi) <- 1:ncol(pi)
  }
  
  ## GENERATE NODES CLASS BELONGING
  names(alpha) <- cl.levels
  if (!fixed.edges) {
    clusters <- factor(sample(cl.levels,p,prob=alpha,replace=TRUE),levels=cl.levels)
  } else {
    cl.prop  <- round(p*alpha)
    clusters <- factor(sample(rep(cl.levels,cl.prop)),levels=cl.levels)
  }
  ## get all couple of edges with class annotations
  couple <- expand.grid(clusters,clusters)
  
  ## ADJACENCY MATRIX GENERATION
  if (directed) {
    ## couples possibly connected for directed graph (non dust)
    possible.edges <- rep(TRUE,nrow(couple))
  } else {
    ## get upper tri indices
    up <- upper.tri(A,diag=FALSE)
    ## couples possibly connected for undirected graph (upper tri non dust)
    possible.edges <- rep(TRUE,nrow(couple)) & up
  }

  for (q in 1:Q) {
    for (l in 1:Q) {
      ## indices of couples (q,l)
      ql <- which(couple[,1] == q & couple[,2] == l & possible.edges)      
      if (!fixed.edges) {
        ## sample the edges according to the probability of connection of the classes
        A[ql] <- rbinom(length(ql),1,pi[q,l])
      } else {
        ## distribute the edges among the K available according to the
        ## density of the classes
        pi[q,l] <- min(pi[q,l],length(ql))
        if (pi[q,l]> 0){A[sample(ql,pi[q,l])] <- 1}
      }
    }
  } 
  ## UNDIRECTED GRAPH SYMMETRIZATION
  if (!directed) {
    A <- (A | t(A)) * 1    
  }

  ## names of A according to the nodes' class
  dimnames(A) <- list(clusters,clusters)
  
  ## GENERATING THE MATRIX OF PARAMETERS

  ## first, the sign pattern (if required) 
  if (signed) {
    switchsign <- matrix(sample(c(1,-1),p*p,replace=TRUE),p,p)
  } else {
    switchsign <- 1
  }

  ## THE THETA MATRIX
  ## either the concentration matrix for undirected network (signed laplacian)
  if (!directed) {
    ## The concentration matrix is the signed Laplacian of A
    L     <- laplacian(A)
    Theta <- - switchsign * t(switchsign) * L
    diag(Theta) <- rep(1,p)
  } else {
    ## The VAR1 matrix is the is normalized by the largest eigen value to
    ## keep the system stable
    Theta <- A
    corr <- runif(p*p)
    Theta[Theta != 0] <- corr[rank(corr)>sum(A==0)] * ((switchsign * A)[A != 0])
    if (!(max(Mod(eigen(Theta)$values))<1))  {
      Theta <- Theta/(max(Mod(eigen(Theta)$values))+1e-9)
    }
  }
  nodes          <- as.character(paste("g",1:p,sep=""))
  dimnames(Theta) <- list(nodes,nodes)

  return(structure(list(A        = A,
                        Theta    = Theta,
                        directed = directed,
                        clusters = clusters,
                        name     = name), class="simone.network"))
}

laplacian <- function(M){
  ## INPUT
  ##  M : adjacency matrix
  ## OUTPUT
  ##  L:  normalized laplacian
  n <- nrow(M)
  diag(M) <- 1
  D <- colSums(M)
  L <- diag(rep(1,n)) -  diag(D^(-1/2))%*% M %*% diag(D^(-1/2))
  L[is.nan(L)] <- 0
  return(L)
}

rTranscriptData <- function (n, graph, ..., mu=rep(0,p), sigma=0.1) {
  
  ## Collect all the graphs in the same list
  graphs  <- list(graph, ...)

  for (g in graphs) {
    if (!is.simone.network(g)) {
      stop("Not a (simone) network object")
    }
  }
  
  ## Check that all graphs are of the same size
  graphsize  <- sapply(graphs,function(g) return(ncol(g$A)))
  if (all(graphsize == graphsize[1])) {
    p  <- graphsize[1]
  } else {
    stop("\nAll graphs must share the same number of nodes.")
  }

  ngraphs   <- length(graphs)
  ndirected <- sapply(graphs,function(g) return(g$directed))

  ## Some initializations
  if (length(n) == 1) {
    ns <- rep(n,ngraphs)
  } else
  if (length(n) == ngraphs) {
    ns <- n
  } else {
    ns <- rep(n[1],ngraphs)
  }
  data <- c()
  
  ## UNDIRECTED GRAPHS -> STATIC DATA
  if (sum(ndirected) == 0) {
    k <- 1
    for (g in graphs) {

      ## Sigma is the pseudo inverse of Theta for undirected graphs
      eig <- eigen(g$Theta)
      D <- eig$values
      P <- eig$vectors
      
      ## Removing the null eigen value
      Dminus <- rep(0,p)
      Dminus[D>1e-12] <- 1/D[D>1e-12]
      Dminus <- diag(Dminus)
      Sigma <- P %*% Dminus %*% t(P)

      ## simulate data matrix and stock the results
      Y <- mu %*% t(rep(1,ns[k])) +
        t(chol(Sigma)) %*% matrix(rnorm(ns[k]*p),p,ns[k]) + rnorm(p,0,sigma)
      data <- rbind(data,t(Y))
      
      k <- k+1
    }
    tasks <- factor(rep(1:ngraphs,ns))    
  } else
  
  ## DIRECTED GRAPHS -> DYNAMIC DATA
  if (sum(ndirected) == ngraphs) {
    k <- 1
    for (g in graphs) {

      ## Xt / VAR1 generation
      X <- matrix(0,ns[k],p)
      X[1,] <- rnorm(p) %*% g$Theta + mu + rnorm(p,0,sigma)
      for (t in 2:ns[k]){
        X[t,] = X[(t-1),] %*% g$Theta + mu + rnorm(p,0,sigma)
      }
      
      data <- rbind(data,X)
      k <- k+1
    }
    tasks <- factor(rep(1:ngraphs,ns))    
  } else {
    stop("\nThe specified graphs do not share the same orientation.\n")
  }

  colnames(data) <- colnames(graphs[[1]]$Theta)
  rownames(data) <- 1:sum(ns)
  
  return(list(X = data, tasks = tasks))
}

coNetwork <- function(graph, delta, name="a co-network") {

  if (!is.simone.network(graph)) {
    stop("Not a (simone) network object")
  }
  
  if (delta > sum(graph$A != 0) ) {
    delta <- sum(graph$A != 0)
    cat("\n Can't remove more edges than available: delta is reduced accordinly (=",delta,").\n",sep="")
  }
  
  p <- ncol(graph$A)
  ## Select edges and non edges
  if (graph$directed) {
    ## all over the matrix (directed graph)
    edges <- which(graph$A != 0)
    zeros <- which(graph$A == 0)    
  } else {
    ## in the upper triangle (symmetric matrix / undirected graph)
    edges <- which(graph$A[upper.tri(graph$A)] != 0)
    zeros <- which(graph$A[upper.tri(graph$A)] == 0)
  }
  nb.edges <- length(edges)
  nb.zeros <- length(zeros)
  
  co.g <- graph

  ## Get some random edges to remove and to add
  edges.to.remove <- sample(edges,delta)
  edges.to.add    <- sample(zeros,delta)

  if (graph$directed) {
    co.g$A[edges.to.remove] <- 0
    co.g$A[edges.to.add]    <- sample(1,delta,replace=TRUE)
    signs <- sign(co.g$Theta)
    signs[edges.to.add] <- sample(c(-1,1),delta,replace=TRUE)

  } else {
    co.g$A[upper.tri(co.g$A)][edges.to.remove] <- 0
    co.g$A[upper.tri(co.g$A)][edges.to.add]    <- sample(1,delta,replace=TRUE)
    tmp <- t(co.g$A)
    co.g$A[lower.tri(co.g$A)] <- tmp[lower.tri(tmp)]
    signs <- sign(co.g$Theta)
    signs[upper.tri(signs)][edges.to.add] <- sample(c(-1,1),delta,replace=TRUE)
    tmp <- t(signs)
    signs[lower.tri(signs)] <- tmp[lower.tri(tmp)]
  }

  
  ## THE THETA MATRIX
  ## either the concentration matrix for undirected network (signed laplacian)
  if (!co.g$directed) {
    ## The concentration matrix is the signed Laplacian of A
    L     <- laplacian(co.g$A)
    Theta <- - signs * L
    diag(Theta) <- rep(1,p)
  } else {
    ## The VAR1 matrix is the is normalized by the largest eigen value to
    ## keep the system stable
    Theta <- co.g$A
    corr <- runif(p*p)
    Theta[Theta != 0] <- corr[rank(corr)>sum(co.g$A==0)]*(signs[co.g$A != 0])
    if (!(max(Mod(eigen(Theta)$values))<1))  {
      Theta <- Theta/(max(Mod(eigen(Theta)$values))+1e-9)
    }
  }
  nodes           <- as.character(paste("g",1:p,sep=""))
  dimnames(Theta) <- list(nodes,nodes)
  
  A <- sign(Theta)
  if (!graph$directed)
    diag(A) <- 0

  return(structure(list(A        = A,
                        Theta    = Theta,
                        directed = graph$directed,
                        clusters = graph$clusters,
                        name     = name), class="simone.network"))
}
