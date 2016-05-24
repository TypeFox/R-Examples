#' Linear MRCA function
#'
#' Function to make the most recent common ancestor (MRCA) matrix of a tree, where entry (i,j) gives the MRCA of tips i and j.
#' The function is linear, exploiting the fact that the tree is rooted.
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree an object of the class \code{phylo} which should be rooted.
#' @param k (optional) number of tips in tree, for faster computation
#'
#' @importFrom phangorn Descendants
#' @importFrom phangorn Children
#' @importFrom combinat combn
#' @importFrom compiler cmpfun
#'
#' @examples
#'
#' ## generate a random tree
#' x <- rtree(6)
#'
#' ## create matrix of MRCAs: entry (i,j) is the node number of the MRCA of tips i and j
#' linearMrca(x,6)
#'
#'
#' @export
linearMrca <- function(tree,k=0) { # k is number of tips, which can be passed to the function to save on computation
  if(!is.rooted(tree)){stop("This function requires the tree to be rooted")}
  if (k==0) {k <- length(tree$tip.label)}
  M <- matrix(0, nrow=k, ncol=k); # initialise matrix
  T <- tree$Nnode # total number of internal nodes
  # traverse internal nodes from root down
  for (tmp in (k+1):(k+T)){
    # find the children of tmp. Then tmp is the MRCA of all pairs of tips descending from different children
    tmp.desc <- Children(tree,tmp)
    Desc <- sapply(1:length(tmp.desc), function(x) Descendants(tree,tmp.desc[[x]],type="tips"))
      if (length(tmp.desc)==2) {  # tmp is the MRCA of tips descending from child one and tips from child two
        I <- Desc[[1]]; J <- Desc[[2]]
          for (i in I)  {
           for (j in J)  {
            M[i,j] <- M[j,i] <- tmp
           }
          }
      }
      else { # for each pair of children of tmp, tmp is the MRCA of their descendant tips
        pairs <- combn(length(Desc),2)
        for (p in 1:length(pairs[1,])) {
          for (i in Desc[[pairs[1,p]]])  {
           for (j in Desc[[pairs[2,p]]])  {
            M[i,j] <- M[j,i] <- tmp
           }
          }
        }
      }
  }
  diag(M) <- 1:k # we define the diagonal elements of M to be the tips themselves
  return(M)
}
linearMrca <- compiler::cmpfun(linearMrca) # compile


#' Tree vector function
#'
#' Function which takes an object of class phylo and outputs the vector for the Kendall Colijn metric.
#' The elements of the vector are numeric if \code{return.lambda.function=FALSE} (default),
#' and otherwise they are functions of lambda.
#'
#' @author Jacob Almagro-Garcia \email{nativecoder@@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree an object of the class \code{phylo}
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{return.lambda.function=TRUE}.
#' @param return.lambda.function If true, a function that can be invoked with different lambda values is returned. This function returns the vector of metric values for the given lambda.
#' @param emphasise.tips an optional list of tips whose entries in the tree vector should be emphasised. Defaults to \code{NULL}.
#' @param emphasise.weight applicable only if a list is supplied to \code{emphasise.tips}, this value (default 2) is the number by which vector entries corresponding to those tips are emphasised.
#'
#' @return The vector of values according to the metric, or a function that produces the vector given a value of lambda.
#'
#' @import ape
#' @importFrom Rcpp evalCpp
#' @importFrom combinat combn2
#' @useDynLib treescape
#'
#' @examples
#'
#' ## generate a random tree
#' tree <- rtree(6)
#' ## topological vector of mrca distances from root:
#' treeVec(tree)
#' ## vector of mrca distances from root when lambda=0.5:
#' treeVec(tree,0.5)
#' ## vector of mrca distances as a function of lambda:
#' vecAsFunction <- treeVec(tree,return.lambda.function=TRUE)
#' ## evaluate the vector at lambda=0.5:
#' vecAsFunction(0.5)
#' 
#'
#' @export
treeVec <- function(tree, lambda=0, return.lambda.function=FALSE, emphasise.tips=NULL, emphasise.weight=2) {
  if(lambda<0 || lambda>1) stop("Pick lambda in [0,1]")
  if(class(tree)!="phylo") stop("Tree should be of class phylo")
  if(is.rooted(tree)!=TRUE) stop("Metric is for rooted trees only")
  if(is.null(tree$edge.length)) {
    warning("Tree edge lengths are not defined, setting edges to have length 1")
	tree$edge.length <- rep(1,length(tree$edge))
	}

  num_leaves <- length(tree$tip.label)
  num_edges <- nrow(tree$edge)

  # We work with ordered labels, using this vector to transform indices.
  tip_order <- match(1:num_leaves, order(tree$tip.label))

  # Ordering the edges by first column places the root at the bottom.
  # Descendants will be placed always before parents.
  edge_order <- ape::reorder.phylo(tree, "postorder", index.only=T)
  edges <- tree$edge[edge_order,]
  edge_lengths <- tree$edge.length[edge_order]
  
  # To emphasise the position of certain tips, if needed:
  if (is.null(emphasise.tips)==FALSE){
    # translate important tip label names into order
    emphasise.tips.order <- tip_order[which(tree$tip.label%in%emphasise.tips)]
    # find the positions where these tips appear in the k choose 2 elements of the final vector
    pairs <- combn2(1:num_leaves)
    emphasise.vector.elements <- union(which(pairs[,1]%in%emphasise.tips.order ),which(pairs[,2]%in%emphasise.tips.order ))
    # make a vector to multiply positions concerning important tips by chosen weighting
    tip.weighting <- c(sapply(1:length(pairs[,1]), function(x) if(x%in%emphasise.vector.elements){emphasise.weight} else{1}),
                       sapply(1:num_leaves, function(y) if(y%in%emphasise.tips.order){emphasise.weight} else{1}))
  }
  else{tip.weighting <- rep(1,0.5*num_leaves*(num_leaves+1))}
  
  
    # We annotated the nodes of the tree in this list. In two passes we are going to
  # compute the partition each node induces in the tips (bottom-up pass) and the distance
  # (in branch length and number of branches) from the root to each node (top-down pass).
  annotated_nodes <- list()

  # Bottom up (compute partitions, we store the branch lengths to compute distances
  # to the root on the way down).
  for(i in 1:num_edges) {

    parent <- edges[i,1]
    child <- edges[i,2]

    # Initialization (leaves).
    if(child <= num_leaves) {
      # We translate the index for the sorted labels.
      child <- tip_order[child]
      # Leaves have as children themselves.
      annotated_nodes[[child]] <- list(root_distance=NULL, edges_to_root=1, partitions=list(child))
    }

    # Aggregate the children partitions (only if we are not visiting a leaf).
    aggregated_partitions <- annotated_nodes[[child]]$partitions[[1]]
    if((child > num_leaves)) {
      for(p in 2:length(annotated_nodes[[child]]$partitions))
        aggregated_partitions <- c(aggregated_partitions, annotated_nodes[[child]]$partitions[[p]])
    }

    # Update the branch length on the child.
    annotated_nodes[[child]]$root_distance <- edge_lengths[i]

    # We have not visited this internal node before.
    if(parent > length(annotated_nodes) || is.null(annotated_nodes[[parent]])) {
      # Assume the first time we get the left child partition.
      annotated_nodes[[parent]] <- list(root_distance=NULL, edges_to_root=1, partitions=list(aggregated_partitions))
    }
    # This is not the first time we have visited the node.
    else {
      # We store the next partition of leaves.
      annotated_nodes[[parent]]$partitions[[length(annotated_nodes[[parent]]$partitions)+1]] <- aggregated_partitions
    }
  }

  # Update the distance to the root at the root (i.e. 0)
  # And the number of edges to the root (i.e. 0).
  annotated_nodes[[num_leaves+1]]$root_distance <- 0
  annotated_nodes[[num_leaves+1]]$edges_to_root <- 0

  # Top down, compute distances to the root for each node.
  for(i in num_edges:1) {
    parent <- edges[i,1]
    child <- edges[i,2]

    # If the child is a leaf we translate the index for the sorted labels.
    if(child <= num_leaves)
      child <- tip_order[child]

    annotated_nodes[[child]]$root_distance <- annotated_nodes[[child]]$root_distance + annotated_nodes[[parent]]$root_distance
    annotated_nodes[[child]]$edges_to_root <- annotated_nodes[[child]]$edges_to_root + annotated_nodes[[parent]]$edges_to_root
  }

  # Distance vectors
  vector_length <- (num_leaves*(num_leaves-1)/2) + num_leaves
  length_root_distances <- double(vector_length)
  topological_root_distances <- integer(vector_length)

  # Fill-in the leaves (notice the involved index translation for leaves).
  topological_root_distances[(vector_length-num_leaves+1):vector_length] <- 1
  length_root_distances[(vector_length-num_leaves+1):vector_length] <- edge_lengths[match(1:num_leaves, edges[,2])][order(tree$tip.label)]

  # Instead of computing the lexicographic order for the combination pairs assume we
  # are filling in a symmetric distance matrix (using only the triangular upper part).
  # We just need to "roll" the matrix indices into the vector indices.
  # Examples for (k=5)
  # The combination c(1,4) would be located at position 3 on the vector.
  # The combination c(2,1) would be located at position 1 on the vector because d(2,1) = d(1,2).
  # The combination c(2,3) would be located at position 5 on the vector.

  index_offsets <- c(0, cumsum((num_leaves-1):1))

  # This is the slow part, we compute both vectors as gain would be marginal.
  sapply(annotated_nodes, function(node) {

    # We skip leaves and the root (if the MRCA for M groups of leaves is at the root
    # all combinations of leaves -among different groups- have 0 as distance to the root).
    # For large trees this can spare us of computing a lot of combinations.
    # Example: In a perfectly balanced binary tree (N/2 leaves at each side of the root),
    # at the root we'd save (N/2) * (N/2) combinations to update. Worst case scenario is completely
    # unbalanced tree (N-1,1), we'd save in that case only N-1 combinations.

    # Make sure we are not visiting a leaf or the root.
    if(length(node$partitions) > 1 && node$root_distance > 0) {

      # Update all combinations for pairs of leaves from different groups.
      num_groups <- length(node$partitions)
      for(group_a in 1:(num_groups-1)) {
        for(group_b in (group_a+1):num_groups) {
          updateDistancesWithCombinations(length_root_distances, topological_root_distances, node$partitions[[group_a]],
                                  node$partitions[[group_b]], index_offsets, node$root_distance, node$edges_to_root)
        }
      }

    }
  })

  if(!return.lambda.function)
    return(tip.weighting * (lambda * length_root_distances + (1-lambda) * topological_root_distances))
  else {
    return(function(l) {
      if(l<0 || l>1) stop("Pick lambda in [0,1]")
      return(tip.weighting * (l * length_root_distances + (1-l) * topological_root_distances)) })
  }
}




#' Metric function
#'
#' Comparison of two trees using the Kendall Colijn metric
#'
#' @author Jacob Almagro-Garcia \email{nativecoder@@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree.a an object of the class \code{phylo}
#' @param tree.b an object of the class \code{phylo} (with the same tip labels as tree.a)
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{return.lambda.function=TRUE}.
#' @param return.lambda.function If true, a function that can be invoked with different lambda values is returned.
#'  This function returns the vector of metric values for the given lambda.
#' @param emphasise.tips an optional list of tips whose entries in the tree vectors should be emphasised. Defaults to \code{NULL}.
#' @param emphasise.weight applicable only if a list is supplied to \code{emphasise.tips}, this value (default 2) is the number by which vector entries corresponding to those tips are emphasised.
#'
#' @return The distance between the two trees according to the metric for the given value of lambda, or a function that produces the distance given a value of lambda.
#'
#'
#' @import ape
#'
#'
#' @examples
#'
#' ## generate random trees
#' tree.a <- rtree(6)
#' tree.b <- rtree(6)
#' treeDist(tree.a,tree.b) # lambda=0
#' treeDist(tree.a,tree.b,1)  # lambda=1
#' dist.func <- treeDist(tree.a,tree.b,return.lambda.function=TRUE) # distance as a function of lambda
#' dist.func(0) # evaluate at lambda=0. Equivalent to treeDist(tree.a,tree.b).
#' ## We can see how the distance changes when moving from focusing on topology to length:
#' plot(sapply(seq(0,1,length.out=100), function(x) dist.func(x)), type="l",ylab="",xlab="")
#'
#' ## The distance may also change if we emphasise the position of certain tips:
#' plot(sapply(tree.a$tip.label, function(x) treeDist(tree.a,tree.b,emphasise.tips=x)),
#'      xlab="Tip number",ylab="Distance when  vector entries corresponding to tip are doubled")
#'
#'
#' @export
treeDist <- function(tree.a, tree.b, lambda=0, return.lambda.function=FALSE, emphasise.tips=NULL, emphasise.weight=2) {

    if(length(tree.a$tip.label) != length(tree.b$tip.label)) stop("Trees must have the same number of tips")

    if(setequal(tree.a$tip.label,tree.b$tip.label) == FALSE) stop("Trees must have the same tip label sets")

    metric_a <- treeVec(tree.a, lambda, return.lambda.function, emphasise.tips, emphasise.weight)
    metric_b <- treeVec(tree.b, lambda, return.lambda.function, emphasise.tips, emphasise.weight)
    if(!return.lambda.function) {
        return(sqrt(sum((metric_a - metric_b)^2)))
    }
    else {
        return(function(l) {
            return(sqrt(sum((metric_a(l) - metric_b(l))^2)))
        })
    }
}


#' Metric function for \code{multiPhylo} input
#'
#' Comparison of a list of trees using the Kendall Colijn metric. Output is given as a pairwise distance matrix. This is equivalent to the \code{$D} output from \code{treescape} but may be preferable for large datasets, and when principal co-ordinate analysis is not required. It includes an option to save memory at the expense of computation time.
#'
#' @author Jacob Almagro-Garcia \email{nativecoder@@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param trees an object of the class \code{multiPhylo} containing the trees to be compared
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{return.lambda.function=TRUE}.
#' @param return.lambda.function If true, a function that can be invoked with different lambda values is returned.
#'  This function returns the matrix of metric values for the given lambda.
#' @param save.memory A flag that saves a lot of memory but increases the execution time (not compatible with return.lambda.function=TRUE).
#' @param emphasise.tips an optional list of tips whose entries in the tree vectors should be emphasised. Defaults to \code{NULL}.
#' @param emphasise.weight applicable only if a list is supplied to \code{emphasise.tips}, this value (default 2) is the number by which vector entries corresponding to those tips are emphasised.
#'
#' @return The pairwise tree distance matrix or a function that produces the distance matrix given a value for lambda.
#'
#'
#' @import ape
#' @importFrom fields rdist
#' @importFrom stats as.dist
#'
#'
#' @examples
#'
#' ## generate 10 random trees, each with 6 tips
#' trees <- rmtree(10,6)
#'
#' ## pairwise distance matrix when lambda=0
#' multiDist(trees)
#'
#' ## pairwise distance matrix as a function of lambda:
#' m <- multiDist(trees, return.lambda.function=TRUE)
#'
#' ## evaluate at lambda=0. Equivalent to multiDist(trees).
#' m0 <- m(0)
#'
#' ## save memory by recomputing each tree vector for each pairwise tree comparison (for fixed lambda):
#' m0.5 <- multiDist(trees,0.5,save.memory=TRUE)
#'
#'
#' @export
multiDist <- function(trees, lambda=0,
                      return.lambda.function=FALSE, save.memory=FALSE,
                      emphasise.tips=NULL, emphasise.weight=2) {
  
  if(!inherits(trees, "multiPhylo")) stop("trees should be a multiphylo object")
  num_trees <- length(trees) 
  if(num_trees<2) {
    stop("multiDist expects at least two trees")
  }
  
  # make name labels well defined
  if(is.null(names(trees))) names(trees) <- 1:num_trees 
  else if(length(unique(names(trees)))!=num_trees){
    warning("duplicates detected in tree labels - using generic names")
    names(trees) <- 1:num_trees
    }
  lab <- names(trees)
  
  # check all trees have same tip labels
  for (i in 1:num_trees) {
    if (!setequal(trees[[i]]$tip.label,trees[[1]]$tip.label)) {
      stop(paste0("Tree ",lab[[i]]," has different tip labels from the first tree."))
    } 
  }
  

  # Working with numbers (no functions).
  if(!return.lambda.function) {

    # Here we speed up the computation by storing all vectors (a lot of memory for big trees).
    if(!save.memory) {
      # Compute the metric vector for all trees.
      tree_metrics <- t(sapply(trees, function(tree) {treeVec(tree, lambda, F, emphasise.tips, emphasise.weight)}))
      distances <- rdist(tree_metrics)
    }

    # To save memory we recompute the vectors for each tree comparison (way slower but we don't eat a ton of memory).
    else {
      distances <- matrix(0.0, num_trees, num_trees)
      
      sapply(1:(num_trees-1), function(i) {
        sapply((i+1):num_trees, function(j) {
          distances[i,j] <<- distances[j,i] <<- treeDist(trees[[i]], trees[[j]], lambda, F, emphasise.tips, emphasise.weight)
        })
      })
    }

    return(as.dist(distances))
  }

  # Working with functions.
  else {

    if(save.memory)
      warning("save.memory=TRUE is incompatible with return.lambda.function=TRUE, setting save.memory=FALSE")

    # Compute the list of metric functions for all trees.
    tree_metric_functions <- sapply(trees, function(tree) {treeVec(tree, lambda, T, emphasise.tips, emphasise.weight)})

    # Inner function that we'll return, computes the distance matrix given lambda.
    compute_distance_matrix_function <- function(l) {
      distances <- matrix(0.0, num_trees, num_trees)
      sapply(1:(num_trees-1), function(i) {
        sapply((i+1):num_trees, function(j) {
          distances[i,j] <<- distances[j,i] <<- sqrt(sum((tree_metric_functions[[i]](l) - tree_metric_functions[[j]](l))^2))
        })
      })
      return(as.dist(distances))
    }
    return(compute_distance_matrix_function)
  }
}

#' Metric function for comparing a reference \code{phylo} to \code{multiPhylo} input
#'
#' Comparison of a single reference tree to a list of trees using the Kendall Colijn metric. Output is given as a vector of distances from the reference tree.
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param refTree a tree of class \code{phylo}, the "reference tree".
#' @param trees an object of the class \code{multiPhylo} containing the trees to be compared to the reference tree
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{return.lambda.function=TRUE}.
#' @param return.lambda.function If true, a function that can be invoked with different lambda values is returned.
#'  This function returns the vector of metric values for the given lambda.
#' @param emphasise.tips an optional list of tips whose entries in the tree vectors should be emphasised. Defaults to \code{NULL}.
#' @param emphasise.weight applicable only if a list is supplied to \code{emphasise.tips}, this value (default 2) is the number by which vector entries corresponding to those tips are emphasised.
#'
#' @return The vector of distances, where entry i corresponds to the distance between the i'th tree and the reference tree, or a function that produces the vector of distances given a value for lambda.
#'
#'
#' @import ape
#' @importFrom stats as.dist
#'
#'
#' @examples
#'
#' ## generate a single reference tree with 6 tips
#' refTree <- rtree(6)
#'
#' ## generate 10 random trees, each with 6 tips
#' trees <- rmtree(10,6)
#'
#' ## find the distances from each of the 10 random trees to the single reference tree
#' refTreeDist(refTree,trees)
#'
#' @export
refTreeDist <- function(refTree, trees, lambda=0, return.lambda.function=FALSE, 
                        emphasise.tips=NULL, emphasise.weight=2) {
  
  if(!inherits(refTree, "phylo")) stop("refTree should be a phylo object")
  if(!inherits(trees, "multiPhylo")) stop("trees should be a multiphylo object")
  num_trees <- length(trees) 

  # make name labels well defined
  if(is.null(names(trees))) names(trees) <- 1:num_trees 
  else if(length(unique(names(trees)))!=num_trees){
    warning("duplicates detected in tree labels - using generic names")
    names(trees) <- 1:num_trees
  }
  lab <- names(trees)
  
  # check all trees have same tip labels
  for (i in 1:num_trees) {
    if (!setequal(trees[[i]]$tip.label,refTree$tip.label)) {
      stop(paste0("Tree ",lab[[i]]," has different tip labels from the reference tree."))
    } 
  }

  # Working with numbers (no functions).
  if(!return.lambda.function) {
    # compute reference tree vector, which will be used repeatedly
    refVec <- treeVec(refTree, lambda, F, emphasise.tips, emphasise.weight)
    
    # for each tree, compute its vector and store Euclidean distance from refVec
    distances <- sapply(trees, function(x) {
      tmpVec <- treeVec(x, lambda, F, emphasise.tips, emphasise.weight)
      sqrt(sum((refVec-tmpVec)^2))
      })
  return(as.vector(distances))
  }
  
  # Working with functions.
  else {
    # compute reference tree vector function 
    refVec <- treeVec(refTree, lambda, T, emphasise.tips, emphasise.weight)
    # Compute the list of metric functions for all trees.
    treeVecs <- sapply(trees, function(x) treeVec(x, lambda, T, emphasise.tips, emphasise.weight))
   
    # Inner function that we'll return, computes the distance matrix given lambda.
    compute_distance_vector_function <- function(l) {
      sapply(1:num_trees, function(x) sqrt(sum((refVec(l)-treeVecs[[x]](l))^2)))
    }
    return(compute_distance_vector_function)
  }
}