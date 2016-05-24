#'Point estimate of the partition using the F-measure as the cost function.
#'
#'Get a point estimate of the partition using the F-measure as the cost function.
#'
#'@param c a list of vector of length \code{n}. \code{c[[j]][i]} is
#'the cluster allocation of observation \code{i=1...n} at iteration
#'\code{j=1...N}.
#'
#'@param logposterior a vector of logposterior corresponding to each
#'partition from \code{c} used to break ties when minimizing the cost function
#'
#'@return a \code{list}:
#'  \itemize{
#'      \item{\code{c_est}:}{ a vector of length \code{n}. Point estimate of the partition}
#'      \item{\code{cost}:}{ a vector of length \code{N}. \code{cost[j]} is the cost
#'      associated to partition \code{c[[j]]}}
#'      \item{\code{similarity}:}{  matrix of size \code{n x n}. Similarity matrix
#'      (see \link{similarityMat})}
#'      \item{\code{opt_ind}:}{ the index of the optimal partition
#'      among the MCMC iterations.}
#'  }
#'
#'
#'@author Francois Caron, Boris Hejblum
#'
#'@export
#'
#'@seealso \code{\link{similarityMat}}
#'

cluster_est_Fmeasure <- function(c, logposterior){

    cat("Estimating posterior F-measures matrix...\n(this may take some time, complexity in O(n^2))\n")
    cmat <- sapply(c, "[")
    tempC <- Fmeasure_costC(cmat)
    cat("DONE!\n")

    Fmeas <- tempC$Fmeas
    cost <- tempC$cost


    opt_ind <- which(cost==min(cost)) #don't use which.min because we want the MCMC iteration maximizing logposterior in case of ties
    if(length(opt_ind)>1){
        opt_ind <- opt_ind[which.max(logposterior[opt_ind])]
    }
    c_est <- c[[opt_ind]]

    return(list("c_est"=c_est, "cost"=cost, "Fmeas"=Fmeas, "opt_ind"=opt_ind))
}