#' Stairs
#' 
#' Calculates the staircase-ness measure.
#' 
#' @author Michael Boyd \email{mboyd855@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @return Two numbers corresponding to the two staircase-ness measures for a tree. Defined in Norstrom et al., Evolutionary Bioinformatics online, 8:261 (2012) \doi{10.4137/EBO.S9738}, these are two related measures:
#' \itemize{ 
#' \item 1: the proportion of subtrees that are imbalanced (i.e. subtrees where the left child has more tip descendants than the right child, or vice versa)
#' \item 2: the average of all the min(l,r)/max(l,r) values of each subtree, where l and r are the number of tips in the left and right children of a subtree.
#' }
#' 
#' @import ape
#'   
#' @examples
#' ## Find the staircase-ness measures in a random tree with 20 tips:
#' stairs(rtree(20))
#'  
#' 
#' @export
stairs <- function(tree)  {
  N <- length(tree$tip.label)
  NDs <- treeImb(tree)[(N + 1):(2 * N - 1),]
  stair1 <- (1/(N - 1)) * sum(abs(NDs[2, ] - NDs[1, ]))
  stair2 <- (1/(N - 1)) * sum(pmin(NDs[2, ], NDs[1, ])/pmax(NDs[2, ], NDs[1, ]))
  return(c(stair1, stair2))
}