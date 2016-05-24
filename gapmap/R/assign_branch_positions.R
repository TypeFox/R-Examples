#'Re-evaluate the position of branches
#'
#'This function reevaluate the position of branches based on the gaps calculated.
#'
#'
#' @param d dendrogram class object
#' @param verbose logical for whether in verbose mode or not
#' @param ... ignored
#' @export assign_branch_positions
#' @aliases assign_branch_positions
#' @return the reevaluated dendrogram object
#' @keywords internal
#' 

assign_branch_positions <- function(d, verbose=FALSE, ...){
  left = d[[1]]
  right = d[[2]]
  if(is.leaf(left) && is.leaf(right)){
    #set attribute for d
    xleft = attr(left, "xpos")
    xright = attr(right, "xpos")
    attr(d,"xleft") = xleft
    attr(d,"xright") = xright
    attr(d,"xmid") = (xleft +  xright)/2
    attr(d,"midpoint") = (xleft +  xright)/2
  }else if(!is.leaf(left) && is.leaf(right)){
    #set the lower branch posistion first
    left = assign_branch_positions(left, verbose=verbose)
    d[[1]] = left
    #get left middle point
    xleft = attr(left, "xmid")
    xright = attr(right, "xpos")
    attr(d,"xleft") = xleft
    attr(d,"xright") = xright
    attr(d,"xmid") = (xleft +  xright)/2
    attr(d,"midpoint") = attr(d,"xmid") - attr(left, "xleft")
  }else if(is.leaf(left) && !is.leaf(right)){
    #set the lower branch posistion first
    right = assign_branch_positions(right, verbose=verbose)
    d[[2]] = right
    xleft = attr(left, "xpos")
    #get right middle point
    xright = attr(right, "xmid")
    attr(d,"xleft") = xleft
    attr(d,"xright") = xright
    attr(d,"xmid") = (xleft +  xright)/2
    attr(d,"midpoint") = attr(d,"xmid") - attr(left, "xleft")
  }else{
    #both subtree
    left = assign_branch_positions(left, verbose=verbose)
    d[[1]] = left
    right = assign_branch_positions(right, verbose=verbose)
    d[[2]] = right
    
    xleft = attr(left, "xmid")
    xright = attr(right, "xmid")
    attr(d,"xleft") = xleft
    attr(d,"xright") = xright
    attr(d,"xmid") = (xleft +  xright)/2
    
    most_left = get_most_left_leaf(left)
    attr(d,"midpoint") = attr(d,"xmid") - attr(most_left, "xpos")#most left
  } 
  d
}
