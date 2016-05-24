#' Field code forest algorithm
#'
#' Workhorse function for \code{\link{generateEdits}}
#'
#' @param E an editarray
#' @param vars variable names still to be eliminated from E
#' @param env an environment where all editmatrices will be stored
#'
#' @seealso \code{\link{generateEdits}}, \code{\link{editarray}}
#' 
#' @example ../examples/generateEdits.R
#' @keywords internal
fcf.env <- function(E,totreat,env){
    # add current edits to collection
    if (nrow(E)>0){
        env$i <- env$i  + 1
        U <- c(env$E,E)
        env$E <- U[!isSubset(U),]
    } else {
        # return if there are no more edits
        return()
    }
    
    # divide and conquer
    B <- blocks(E)
    for ( b in B){
        # variables to be treated in the block
        totreatb <- totreat[names(totreat) %in% getVars(b)]
        vrs <- names(totreatb)[totreatb]
        # variables which cannot be resolved need not be treated (prune the tree)
        totreatb[!resolves(b,vrs)] <- FALSE
        vrs <- names(totreatb)[totreatb]
        # order variables so the most connected variables are eliminated first,
        # eliminate variables and recurse.
        vrs <- vrs[order(colSums(contains(b,vrs)),decreasing=TRUE)]
        for ( k in vrs ){
            totreatb[k] <- FALSE
            fcf.env(reduce(eliminate(b,k)),totreatb[-which(names(totreatb)==k)],env)
        }
    }
}



#' Derive all essentially new implicit edits
#'
#' Implements the Field Code Forest (FCF) algorithm of Garfinkel et al (1986) to 
#' derive all essentially new implicit edits from an editarray. The FCF is really
#' a single, highly unbalanced tree. This algorithm traverses the tree, pruning many 
#' unnecessary branches, uses \code{\link{blocks}} to divide and conquer, and 
#' optimizes traversing order. See Van der Loo (2012) for a description
#' of the algorithms.
#'
#' @param E An \code{\link{editarray}}
#' @return A 3-element named \code{list}, where element \code{E} is an \code{\link{editarray}} containing all 
#' generated edits. \code{nodes} contains information on the number of nodes in the tree and vs the number of nodes
#' traversed and \code{duration} contains user, system and elapsed time inseconds.
#' The \code{\link[=editarray]{summary}} method for \code{\link{editarray}} prints this information.
#'
#'
#' @references
#' R.S. Garfinkel, A.S. Kunnathur and G.E. Liepins (1986). 
#'    Optimal imputation of erroneous data: categorical data, general edits.
#'    Operations Research 34, 744-751.
#'
#' M.P.J. Van der Loo (2012). Variable elimination and edit generation with a flavour of 
#' semigroup algebra (submitted)
#'
#' @export
generateEdits <- function(E){
    if ( !is.editarray(E) ) stop("Argument must be of class 'editarray' ")
    t0 <- proc.time()
    # initialize variables to treat
    vars <- getVars(E)
    totreat <- rep(TRUE,length(vars))
    names(totreat) <- vars
    # set up environment to collect generated edits
    e <- new.env()
    # node counter
    e$i <- 0
    e$E <- E[logical(0),]
    # call the workhorse
    fcf.env(E,totreat,e)
    duration <- getDuration(proc.time()-t0)
    # return edits
    return(
        list(
            edits=e$E, 
            nodes=c(nodesInTree = 2^length(vars), nodesTraversed = e$i), 
            duration=duration
        )
    )
}



# Check which variables of 'vars' can be resolved
resolves <- function(E,vars){
    if ( length(vars)==0) return(logical(0))
    ind <- getInd(E)[vars]
    Ic <- contains(E,vars)
    sapply(vars, function(v) all(colSums(E[Ic[,v],ind[[v]],drop=FALSE])>0))
}



