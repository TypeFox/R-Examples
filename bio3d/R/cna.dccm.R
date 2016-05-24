cna.dccm <-  function(cij, cutoff.cij=0.4, cm=NULL,  vnames=colnames(cij),
                  cluster.method="btwn", collapse.method="max", 
                  cols=vmd.colors(), minus.log=TRUE, ...){

    
  ## Check for presence of igraph package
 oops <- requireNamespace("igraph", quietly = TRUE)
 if (!oops) {
    stop("igraph package missing: Please install, see: ?install.packages")
 }

  if (dim(cij)[1] != dim(cij)[2]) {
    stop("Input 'cij' should be a square matrix as obtained from the 'dccm()' function")
  }

  ## Check vnames/colnames if present. These are used to name nodes
  if( is.null( vnames ) ) {
    vnames <- 1:ncol(cij)
  }
  if( length(vnames) != ncol(cij) ) {
    stop("Length of input 'vnames' and number of cols in input 'cij' do not match")
  }
  colnames(cij) <- vnames
  
  ## Check 'cm' contact map if present.
  if(!is.null(cm)){
    if (dim(cm)[1] != dim(cm)[2]) {
      stop("Input 'cm' should be a square contact matrix as obtained from the 'cmap()' function")
    }
    if (any(range(cm, na.rm=T) != c(0,1))) {
      stop("Input 'cm' should be a binary contact matrix as obtained from the 'cmap()' function")
    }
    if (dim(cm)[1] != dim(cij)[1]) {
      stop("Inputs 'cij' and 'cm' should have the same dimensions")
    }
    # convert NAs to 0
    cm[is.na(cm)] = 0
  }

  
  ##-- Functions for later
  cluster.network <- function(network, cluster.method="btwn"){
    
    ## Function to define community clusters from network,
    ##  cluster methods can be one of of 'cluster.options'
    cluster.options=c("btwn", "walk", "greed")
    cluster.method <- match.arg(tolower(cluster.method), cluster.options)
    comms <- switch( cluster.method,
              btwn = igraph::edge.betweenness.community(network, directed=FALSE),
              walk = igraph::walktrap.community(network),
              greed = igraph::fastgreedy.community(network) )
    
    names(comms$membership) <- igraph::V(network)$name
    return(comms)
  }

  contract.matrix <- function(cij.network, membership,## membership=comms$membership,
                              collapse.method="max", minus.log=minus.log){  ## Changed from minus.log=TRUE
    
    ## Function to collapse a NxN matrix to an mxm matrix
    ##  where m is the communities of N. The collapse method
    ##  can be one of the 'collapse.options' below

    ## convert to the original cij values if "-log" was used
    
    if(minus.log){
      cij.network[cij.network>0] <- exp(-cij.network[cij.network>0])
    }
    
    collapse.options=c("max", "median", "mean", "trimmed")
    collapse.method <- match.arg(tolower(collapse.method), collapse.options)

    ## Fill a 'collapse.cij' nxn community by community matrix
    node.num <- max(membership)
    if(node.num > 1){
      collapse.cij <- matrix(0, nrow=node.num, ncol=node.num)
      inds <- pairwise(node.num)

      for(i in 1:nrow(inds)) {
        comms.1.inds <- which(membership==inds[i,1])
        comms.2.inds <- which(membership==inds[i,2])
        submatrix <- cij.network[comms.1.inds, comms.2.inds]

        ## Use specified "collapse.method" to define community couplings
        collapse.cij[ inds[i,1], inds[i,2] ] = switch(collapse.method,
                    max = max(submatrix),
                    median = median(submatrix),
                    mean = mean(submatrix),
                    trimmed = mean(submatrix, trim = 0.1))
      }
      
      if(minus.log){
        collapse.cij[collapse.cij>0] <- -log(collapse.cij[collapse.cij>0])
      }
      
      ## Copy values to lower triangle of matrix and set colnames
      collapse.cij[ inds[,c(2,1)] ] = collapse.cij[ inds ]
      colnames(collapse.cij) <- 1:ncol(collapse.cij)
    }
    else{
      warning("There is only one community in the $communities object.
               $community.cij object will be set to 0 in the 
               contract.matrix() function.")

      collapse.cij <- 0
    }

    class(collapse.cij) <- c("dccm", "matrix")
    return(collapse.cij)
  }
  
  ## Store the command used to submit the calculation
  cl <- match.call()
  
  ##- Take absolute value of 'cij'
  cij.abs <- abs(cij)

  ## Filter: set to 0 all values below the cutoff
  cij.abs[cij.abs < cutoff.cij] = 0

  if(minus.log){
    ##-- Calculate the -log of cij
    ## change cij >= 0.9999 to 0.9999 to avoid numerical problems 
    ##   (-log is too close to zero)

    cij.network <- cij.abs
    cij.network[cij.network >= 0.9999] = 0.9999
    cij.network <- -log(cij.network)
    ## remove infinite values
    cij.network[is.infinite(cij.network)] = 0
  }
  else{
    cij.network <- cij.abs
  }

  if(!is.null(cm)){  
    ##-- Filter cij by contact map
    cij.network <- cij.network * cm
  }

  ##  cij.network contains either the -log(abs.cij) or just abs.cij. 
  ##   (the default is minus.log=TRUE)
  
  ##-- Make an igraph network object
  network <- igraph::graph.adjacency(cij.network,
                             mode="undirected",
                             weighted=TRUE,
                             diag=FALSE)
  
  ##-- Calculate the first set of communities
  communities <- cluster.network(network, cluster.method)

  ##-- Coarse grain the cij matrix to a new cluster/community matrix
  community.cij <- contract.matrix(cij.network, communities$membership, 
                                   collapse.method, minus.log)

  ##-- Generate a coarse grained network --##
  if(sum(community.cij)>0){
    community.network <-  igraph::graph.adjacency(community.cij,
                               mode="undirected",
                               weighted=TRUE,
                               diag=FALSE)

    ##-- Cluster the community network to obtain super-communities -- OLD VERSION
    ## clustered.communities <- cluster.network(community.network, cluster.method)

    ##-- Annotate the two networks with community information
    ## Check for duplicated colors
    if(max(communities$membership) > length(unique(cols)) ) {
      warning("The number of communities is larger than the number of unique 
              'colors' provided as input. Colors will be recycled")
    }
  
    ## Set node colors
    igraph::V(network)$color <- cols[communities$membership]
    igraph::V(community.network)$color <- cols[ 1:max(communities$membership)]
  
    ## Set node sizes
    igraph::V(network)$size <- 1
    igraph::V(community.network)$size <- table(communities$membership)

  } else{
    warning("The $communities structure does not allow a second clustering 
            (i.e. the collapsed community.cij matrix contains only 0). 
            'community.network' object will be set to NA")
      
    community.network <- NA
    clustered.communities <- NA
      
    if(max(communities$membership) > length(unique(cols)) ) {
      warning("The number of communities is larger than the number of unique 
              'colors' provided as input. Colors will be recycled")
    }
  
    ## Set node colors
    igraph::V(network)$color <- cols[communities$membership]

    ## Set node sizes
    igraph::V(network)$size <- 1
  }
  
  ## Output
  output <- list("network"=network,
                 "communities"=communities,
                 "community.network"=community.network,
                 "community.cij"=community.cij,
                 "cij"=cij.network,
                 call = cl)

  class(output)="cna"

  return(output)
}
  
