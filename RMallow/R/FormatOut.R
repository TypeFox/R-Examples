#' Formats the data in the "Solve" function for output.
#' 
#' Data formatting function.
#' 
#' 
#' @param R The modal sequences.
#' @param p Proportion of data in each cluster.
#' @param lambda Mallows' spread parameters for each cluster.
#' @param z Probability of cluster membership for each individual.
#' @param datas Matrix of partial sequences.
#' @param likelihood Vector of the log-likelihood of the model at each
#' iteration.
#' @return \item{R}{The modal sequences} \item{p}{Proportion in each cluster}
#' \item{lambda}{Spread parameters for each cluster} \item{datas}{Rankings
#' merged with their cluster membership, distance from each cluster center, and
#' probability of each cluster membership} \item{min.like}{Likelihood at each
#' iteration}

#' @author Erik Gregory

#' @keywords BubbleSort Kendall
FormatOut <-
function(R, p, lambda, z, datas, likelihood) {
  # Find which cluster the individual belongs to.
  clust <- apply(z, 1, which.max)
  # Change the forms of the sequences.
  seqs <- unlist(lapply(R, 
                        function(i) paste(i, collapse = " ")))
  datas <- data.frame(datas, clust)
  N <- ncol(datas) - 1
  dists <- AllKendall(datas[, 1:N], do.call("rbind", R))
  datas <- data.frame(datas, pvals = z, seq = seqs[clust], dists = dists)
  out <- list(R = R, p = p, lambda = lambda, 
              datas = datas, min.like = likelihood)
  return(out)
}
