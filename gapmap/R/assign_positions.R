#'Re-evaluate the position of leaves
#'
#'This function reevaluate the position of leaves based on the gaps calculated.
#'
#'
#' @param d dendrogram class object
#' @param runningX numerical position of leaf node on display
#' @param verbose logical for whether in verbose mode or not
#' @param ... ignored
#' @export assign_positions
#' @aliases assign_positions
#' @return the reevaluated dendrogram object
#' @keywords internal
#' 


#assign positions
assign_positions <- function(d, runningX=1, verbose=FALSE, ...){
  left = d[[1]]
  right = d[[2]]
  if(is.leaf(left) && is.leaf(right)){
    if(verbose)print("leaf_leaf")
    #both leaves
    attr(left, "xpos") = runningX
    if(verbose) print(paste0("    ","debug: set runningX for ",attr(left,"label") ," is ", runningX))
    runningX = runningX + attr(left, "right_gap") +1
    attr(right, "xpos") = runningX
    if(verbose) print(paste0("    ","debug: set runningX for ",attr(right,"label") ," is ", runningX))
  }else if(!is.leaf(left) && is.leaf(right)){
    if(verbose)print("sub_leaf")
    #set the positions on the left first
    left = assign_positions(left, runningX, verbose=verbose)
    #get the right most leaf of the subtree
    right_most = get_most_right_leaf(left)
    
    runningX = attr(right_most, "xpos") + attr(right_most,"right_gap") + 1
    if(verbose) print(paste0("sub_leaf: runningX for ", attr(right_most,"label"), " is ", runningX,", xpos is ", attr(right_most, "xpos"), ", right_gap is ", attr(right_most, "right_gap")))
    attr(right, "xpos") = runningX
    if(verbose) print(paste0("    ","debug: set runningX for ",attr(right,"label") ," is ", runningX))
  }else if(is.leaf(left) && !is.leaf(right)){
    if(verbose)print("leaf_sub")
    #set the position of the left leaf
    attr(left, "xpos") = runningX
    if(verbose) print(paste0("    ","debug: set runningX for ",attr(left,"label") ," is ", runningX))
    runningX = runningX + attr(left, "right_gap") +1
    #set the position on the right subtree
    right = assign_positions(right, runningX, verbose=verbose)
  }else{
    #both subtree
    if(verbose)print("sub_sub")
    #set the positions on the left first
    left = assign_positions(left, runningX, verbose=verbose)
    #get the right most leaf of the subtree
    right_most = get_most_right_leaf(left)
    runningX = attr(right_most, "xpos") + attr(right_most,"right_gap") + 1
    right = assign_positions(right, runningX, verbose=verbose)
  } 
  d[[1]] = left
  d[[2]] = right
  d
}