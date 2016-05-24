#'Burning MCMC iterations from a Dirichlet Process Mixture Model.
#'
#'Utility function for burning MCMC iteration from a \code{DPMMclust} object.
#'
#'@name burn.DPMMclust
#'
#'@param x a \code{DPMMclust} object.
#'
#'@param burnin the number of MCMC iterations to burn (default is \code{0}).
#'
#'@param thin the spacing at which MCMC iterations are kept.
#'Default is \code{1}, i.e. no thining.
#'
#'@param ... further arguments passed to or from other methods.
#'
#'@return a \code{DPMMclust} object minus the burnt iterations
#'
#'@author Boris Hejblum
#'
#'@export burn.DPMMclust
#'
#'@seealso \code{\link{summary.DPMMclust}}
#'

burn.DPMMclust <- function(x, burnin=0, thin=1){

  if(burnin>=length(x[["mcmc_partitions"]])){
    stop("burnin argument is superior to the number of MCMC iterations sampled")
  }
  xburnt <- list()

  N <- x[["nb_mcmcit"]]

  if(thin>1){
    select <- c(TRUE, rep(FALSE, thin-1))

    xburnt[["mcmc_partitions"]] <- x[["mcmc_partitions"]][(burnin+1):N][select]
    xburnt[["alpha"]] <- x[["alpha"]][(burnin+1):N][select]
    xburnt[["U_SS_list"]] <- x[["U_SS_list"]][(burnin+1):N][select]
    xburnt[["weights_list"]] <- x[["weights_list"]][(burnin+1):N][select]
    xburnt[["logposterior_list"]] <- x[["logposterior_list"]][(burnin+1):N][select]

    if (x$clust_distrib=="gaussian"){
      xburnt[["listU_mu"]] <- x[["listU_mu"]][(burnin+1):N][select]
      xburnt[["listU_Sigma"]] <- x[["listU_Sigma"]][(burnin+1):N][select]
    }

  }else{
    xburnt[["mcmc_partitions"]] <- x[["mcmc_partitions"]][(burnin+1):N]
    xburnt[["alpha"]] <- x[["alpha"]][(burnin+1):N]
    xburnt[["U_SS_list"]] <- x[["U_SS_list"]][(burnin+1):N]
    xburnt[["weights_list"]] <- x[["weights_list"]][(burnin+1):N]
    xburnt[["logposterior_list"]] <- x[["logposterior_list"]][(burnin+1):N]

    if (x$clust_distrib=="gaussian"){
      xburnt[["listU_mu"]] <- x[["listU_mu"]]#[(burnin+1):N]
      xburnt[["listU_Sigma"]] <- x[["listU_Sigma"]][(burnin+1):N]
    }

  }

  xburnt[["hyperG0"]] <- x[["hyperG0"]]
  xburnt[["data"]] <- x[["data"]]
  xburnt[["nb_mcmcit"]] <- (N-burnin)/thin
  xburnt[["clust_distrib"]] <- x[["clust_distrib"]]

  class(xburnt) <- "DPMMclust"
  return(xburnt)
}
