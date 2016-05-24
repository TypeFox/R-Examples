`paths` <- function(trans,start=1)
{
    ### Find all paths through a multi-state model from a starting state
    ### All direct transitions are specified through transition matrix
    ### trans
    ### Input:
    ###     start: numeric, specifying the starting state
    ###     trans: transition matrix
    ### Output:
    ###     matrix specifying all possible paths
    ### Details: recursive
    nstates <- trans[start,]
    nstates <- which(!is.na(nstates))
    if (length(nstates)==0) ## i.e. in absorbing state
        return(matrix(start,1,1))
    else {
        fpmat <- startmat <- matrix(start,1,1)
        for (nstate in nstates) {
            fpmat <- my.rbind(fpmat, my.cbind2(start,Recall(trans,start=nstate)))
        }
    }
    return(fpmat)
}
