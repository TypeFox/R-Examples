#' Check if move is acceptable.
#' 
#' This function takes as input a new network proposal and checks that the
#' proposal does not exceed the maximum number of parents for a node, and that
#' there are no self loops (if self loops have been disallowed).
#' 
#' 
#' @param proposal The proposed network (K-by-q matrix with K segments and q
#' parent sets).
#' @param qmax Maximum number of parents allowed.
#' @param self.loops Flag indicating whether self loops are allowed.
#' @param target The current target node (only needed to find out which parent
#' would be the self loop).
#' @param fixed.edges Which edges in the network should be fixed for all 
#' segments (q-by-q matrix with entries 0 for fixed non-edge, 1 for fixed edge,
#' -1 for non-fixed edge).
#' @return Returns \code{TRUE} if the proposed move is acceptable, \code{FALSE}
#' otherwise.
#' @author Frank Dondelinger
#' @seealso \code{\link{make_structure_move}}
#' @export AcceptableMove
AcceptableMove <-
function(proposal, qmax, self.loops, target, fixed.edges) {
  # Checks if the proposed network is valid (does not exceed the maximum 
  # number of parents per node).
  #
  # Args:
  #  proposal: Proposed network (Kxq matrix with K segments and q 
  #            parent sets
  #  qmax:     Maximum number of parents allowed.
  #
  # Returns:
  #  True if no parent set exceeds the maximum number of parents, false
  #  otherwise
  # Check fan-in restriction
  acceptable = all(apply(proposal, 1, sum) <= qmax)
  
  # Check self-loops
  if(!self.loops) {
    acceptable = acceptable && all(proposal[,target] == 0) 
  }
  
  fixed.segs = t(matrix(c(fixed.edges, -1), dim(proposal)[2], dim(proposal)[1]))
  indices.fixed = fixed.segs >= 0
  
  acceptable = acceptable &&
    all(fixed.segs[indices.fixed] == proposal[indices.fixed])
  
  return(acceptable)
}

