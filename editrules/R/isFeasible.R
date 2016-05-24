#' Check consistency of set of edits
#'
#' When variables are \code{\link[=eliminate]{eliminated}} one by one
#' from a set of edits, eventually either no edits are left or an 
#' \code{\link[=isObviouslyInfeasible]{obvious contradiction}} is encountered.
#' In the    case no records can obey all edits in the set which is therefore
#' \code{inFeasible}.
#'
#'
#' @note This function can potentially take a long time to complete, especially
#' when many connected (conditional) edits are present. Consider using \code{\link{blocks}}
#' to check feasibility of indendent blocks. 
#' 
#'
#' @param E an \code{\link{editmatrix}}, \code{\link{editarray}} or \code{\link{editset}}
#' @param warn logical: should a warning be emitted when system is infeasible?
#' @return TRUE or FALSE
#'
#' @seealso \code{\link{isObviouslyInfeasible}}, \code{\link{isObviouslyRedundant}}
#' @export
isFeasible <- function(E, warn=FALSE){
    ## TODO: make it return the subset of edits causing the contradiction.
    vars <- getVars(E)
    vars2 <- vars
    feasible <- any(!isObviouslyInfeasible(E))
    while( isTRUE(feasible) && length(vars) > 0 ){
        E <- eliminate(E,vars[1])
        vars <- vars[-1]
        ## TODO: cleanup editlists that have infeasible parts, currently they are included
        ## for all eliminations.
        feasible <- any(!isObviouslyInfeasible(E))
        if ( !feasible && warn )
            warning(
                paste("system becomes obviously infeasible after eliminating",
                paste(vars2[!(vars2 %in% vars)],collapse=", "))
            ) 
    }
    return(feasible)
}






