#' Geometric median tree function
#'
#' Finds the geometric median of a set of trees according to the Kendall Colijn metric.
#'
#' @author Jacob Almagro-Garcia \email{nativecoder@@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @param x A list of trees of the class multiPhylo, for which the median tree will be computed, \cr
#' OR a matrix of tree vectors as given by \code{treescape$vectors}.
#' @param groups an optional factor defining groups of trees; if provided, one median tree will be found for each group.
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised. This argument is ignored if \code{return.lambda.function=TRUE} or if the vectors are already supplied as the object \code{x}.
#' @param weights A vector of weights for the trees. Defaults to a vector of 1's so that all trees are equally weighted, but can be used to encode likelihood, posterior probabilities or other characteristics.
#' @param emphasise.tips an optional list of tips whose entries in the tree vectors should be emphasised. Defaults to \code{NULL}.
#' @param emphasise.weight applicable only if a list is supplied to \code{emphasise.tips}, this value (default 2) is the number by which vector entries corresponding to those tips are emphasised.
#' @param return.lambda.function If true, a function that can be invoked with different lambda values is returned.
#'  This function returns the vector of metric values for the given lambda. Ignored if the tree vectors are already supplied as the object \code{x}.
#' @param save.memory A flag that saves a lot of memory but increases the execution time (not compatible with return.lambda.function=TRUE). Ignored if the tree vectors are already supplied as the object \code{x}.
#'
#' @return A list of five objects: 
#' \itemize{
#' \item $centre is the "central vector", that is, the (weighted) mean of the tree vectors (which typically does not correspond to a tree itself); 
#' \item $distances gives the distance of each tree from the central vector; 
#' \item $mindist is the minimum of these distances; 
#' \item $treenumbers gives the numbers (and, if supplied, names) of the "median tree(s)", that is, the tree(s) which achieve this minimum distance to the centre; 
#' \item $trees if trees were supplied then this returns the median trees as a multiPhylo object. 
#' }
#' If groups are provided, then one list is returned for each group.
#' If \code{return.lambda.function=TRUE} then a function is returned that produces this list for a given value of lambda. 
#'
#'
#' @import ape
#'
#'
#' @examples
#'
#' ## EXAMPLE WITH WOODMICE DATA
#' data(woodmiceTrees)
#'
#' ## LOOKING FOR A SINGLE MEDIAN
#' ## get median tree(s)
#' res <- medTree(woodmiceTrees)
#' res
#'
#' ## plot first tree
#' med.tree <- res$trees[[1]]
#' plot(med.tree)
#'
#' ## LOOKING FOR MEDIANS IN SEVERAL CLUSTERS
#' ## identify 6 clusters
#' groves <- findGroves(woodmiceTrees, nf=3, nclust=6)
#'
#' ## find median trees
#' res.with.grp <- medTree(woodmiceTrees, groves$groups)
#'
#' ## there isone output per cluster
#' names(res.with.grp)
#'
#' ## get the first median of each
#' med.trees <- lapply(res.with.grp, function(e) ladderize(e$trees[[1]]))
#'
#' ## plot trees
#' par(mfrow=c(2,3))
#' for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i))
#'
#'
#' @export
medTree <- function(x, groups=NULL, lambda=0, weights=NULL, emphasise.tips=NULL, emphasise.weight=2,
                    return.lambda.function=FALSE, save.memory=FALSE) {

  ## CHECK input type ##
  if (inherits(x, "multiPhylo")) {
    type <- "multiPhylo_object"
    }
  else if (inherits(x, "matrix")) { 
    # x is a matrix of tree vectors
    type <- "tree_vectors"
    } 
  else stop("x should be a multiphylo object or a matrix of tree vectors")

if (type=="multiPhylo_object") {  
    ## DEFINE MAIN FUNCTION FINDING MEDIAN TREE ##
    findMedianPhylo <- function(trees, weights){
        ## checks, general variables
        num_trees <- length(trees)
        num_leaves <- length(trees[[1]]$tip.label)
        if(is.null(weights)) weights <- rep(1, num_trees)
        if(length(weights)!=num_trees) stop("Length of vector of weights must be the same as number of trees")

        ## Working with numbers (no functions).
        if(!return.lambda.function) {

            ## Here we speed up the computation by storing all vectors (a lot of memory for big trees).
            if(!save.memory) {

                ## Compute the metric vector for all trees.
                tree_metrics <- t(sapply(trees, function(tree) {treeVec(tree, lambda, emphasise.tips, emphasise.weight, return.lambda.function=F)}))

                ## Compute the centre metric vector by weighting the metric vector of each tree.
                centre <- (weights %*% tree_metrics)/num_trees

                ## Distances to the centre.
                distances <- apply(tree_metrics, 1, function(m){sqrt(sum((m-centre)^2))})

                ## Get the indices for the median tree(s).
                min_distance <- min(distances)
                median_trees <- which(min_distance == distances)

                return(list(centre=centre, distances=distances, mindist=min_distance, treenumbers=median_trees, trees=trees[median_trees]))
            }

            ## To save memory we recompute the vectors on the fly (way slower but we don't eat a ton of memory).
            ## We'll need a first pass to compute the centre and a second pass to compute distances.
            else {

                ## First pass: compute the centre.
                centre <- rep(0,(num_leaves*(num_leaves-1)/2) + num_leaves)
                for(i in 1:num_trees) {
                    centre <- centre + treeVec(trees[[i]], lambda, F) * weights[i]
                }
                centre <- centre/num_trees

                ## Second pass: compute the distances.
                distances <- rep(NA,num_trees)
                for(i in 1:num_trees) {
                    distances[i] <- sqrt(sum((treeVec(trees[[i]], lambda, F) - centre)^2))
                }

                ## Get the indices for the median tree(s).
                min_distance <- min(distances)
                median_trees <- which(min_distance == distances)

                return(list(centre=centre, distances=distances, mindist=min_distance, treenumbers=median_trees, trees=trees[median_trees]))
            }
        }

        ## Working with functions.
        else {

            if(save.memory)
                warning("save.memory=TRUE is incompatible with return.lambda.function=TRUE, setting save.memory=FALSE")

            ## Compute the list of metric functions for all trees.
            tree_metric_functions <- sapply(trees, function(tree) {treeVec(tree, lambda, emphasise.tips, emphasise.weight, return.lambda.function=T)})

            ## Inner function that we'll return, computes the distance matrix given lambda.
            compute_median_tree_function <- function(l) {

                ## Compute the tree metrics for the given lambda.
                tree_metrics <- t(sapply(tree_metric_functions, function(tmf){tmf(l)}))

                ## Compute the centre metric vector by weighting the metric vector of each tree.
                centre <- (weights %*% tree_metrics)/num_trees

                ## Distances to the centre.
                distances <- apply(tree_metrics, 1, function(m){sqrt(sum((m-centre)^2))})

                ## Get the indices for the median tree(s).
                min_distance <- min(distances)
                median_trees <- which(min_distance == distances)

                return(list(centre=centre, distances=distances, mindist=min_distance, treenumbers=median_trees, trees=trees[median_trees]))
            }

            return(compute_median_tree_function)
        }
    } # end findMedian


    ## APPLY FUNCTION TO TREES ##
    if(is.null(groups)){     ## no groups provided
        out <- findMedianPhylo(x, weights)
    } else { ## groups provided
        out <- tapply(x, groups, findMedianPhylo, weights)
    }
} # end if multiPhylo object

if (type=="tree_vectors"){
  ## Can define a much simpler version of the function to find a median tree ##
  findMedianVectors <- function(vectors, weights){
    ## checks, general variables
    num_trees <- length(vectors[,1])
    if(is.null(weights)) {weights <- rep(1,num_trees)}
    if(length(weights)!=num_trees) stop("Length of vector of weights must be the same as number of tree vectors")

    tree_metrics <- vectors
    
    ## Compute the centre metric vector by weighting the metric vector of each tree.
    centre <- (weights %*% tree_metrics)/num_trees
        
    ## Distances to the centre.
    distances <- apply(tree_metrics, 1, function(m){sqrt(sum((m-centre)^2))})
        
    ## Get the indices for the median tree(s).
    min_distance <- min(distances)
    median_trees <- which(min_distance == distances)
    
    ## Note we cannot return $trees because the trees were not supplied!    
    return(list(centre=centre, distances=distances, mindist=min_distance, treenumbers=median_trees))
    } # end findMedianVectors
  
  
  ## APPLY FUNCTION TO TREES ##
  if(is.null(groups)){     ## no groups provided
    out <- findMedianVectors(x, weights)
  } else { ## groups provided
    # need to first convert the vector matrix into a list:
    mylist <- lapply(1:length(x[,1]), function(a) x[a,])
    # and then coerce back into matrix within the function. 
    # Room for improvement here!
    out <- tapply(mylist, groups, function(a) {findMedianVectors(t(sapply(a, function(b) b)), weights)})
  }
}  
    ## RETURN ##
    return(out)
} ## end medTree