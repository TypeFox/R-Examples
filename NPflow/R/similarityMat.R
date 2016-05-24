#'Computes the coclustering (or similarity) matrix
#'
#'
#'@param c a list of vector of length \code{n}. \code{c[[j]][i]} is 
#'the cluster allocation of observation \code{i=1...n} at iteration 
#'\code{j=1...N}.
#'
#'@param step provide coclustering every \code{step} iterations.
#'Default is 1.
#'
#'@return A matrix of size \code{n x n} whose term \code{[i,j]} 
#'is the proportion of MCMC iterations where observation \code{i} and 
#'observations \code{j} are allocated to the same cluster. 
#'
#'@author Boris Hejblum
#'
#'@export
#'

similarityMat <- function(c, step=1){
    n <- length(c[[1]])
    if(step>1){
        select <- c(TRUE, rep(FALSE, step-1))
        c <- c[select]
    }
    
    vclust2mcoclust <- function(v){
        m <- sapply(v, FUN=function(x){x==v})
        return(m)
    }
    list_mcoclust <- lapply(c, vclust2mcoclust)
    
    #coclustering matrix: out_coclust[i,j] is the proportion of 
    #iterations where c(i)=c(j)
    out_coclust <- Reduce('+', list_mcoclust)
    out_coclust <- out_coclust/length(c)
    
    return(out_coclust)
}