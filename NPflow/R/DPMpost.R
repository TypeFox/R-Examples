#'Posterior estimation for Dirichlet process mixture of multivariate (potentially skew) distibutions models
#'
#'Partially collapse slice Gibbs sampling for Dirichlet process mixture of multivariate
#'normal, skew normal or skew t distributions.
#'
#'@details This function is a wrapper around \code{\link{DPMGibbsN}}, \code{\link{DPMGibbsN_parallel}},
#'\code{\link{DPMGibbsN_SeqPrior}}, \code{\link{DPMGibbsSkewN}}, \code{\link{DPMGibbsSkewN_parallel}},
#'\code{\link{DPMGibbsSkewT}}, \code{\link{DPMGibbsSkewT_parallel}},
#'\code{\link{DPMGibbsSkewT_SeqPrior}}, \code{\link{DPMGibbsSkewT_SeqPrior_parallel}}.
#'
#'@param data data matrix \code{d x n} with \code{d} dimensions in rows
#'and \code{n} observations in columns.
#'
#'@param hyperG0 prior mixing distribution.
#'
#'@param a shape hyperparameter of the Gamma prior
#'on the concentration parameter of the Dirichlet Process. Default is \code{0.0001}.
#'
#'@param b scale hyperparameter of the Gamma prior
#'on the concentration parameter of the Dirichlet Process. Default is \code{0.0001}. If \code{0}, 
#'then the concentration is fixed set to \code{a}.
#'
#'@param N number of MCMC iterations.
#'
#'@param doPlot logical flag indicating wether to plot MCMC iteration or not.
#'Default to \code{TRUE}.
#'
#'@param nbclust_init number of clusters at initialisation.
#'Default to 30 (or less if there are less than 30 observations).
#'
#'@param plotevery an integer indicating the interval between plotted iterations when \code{doPlot}
#'  is \code{TRUE}.
#'
#'@param diagVar logical flag indicating wether the variance of each cluster is
#'estimated as a diagonal matrix, or as a full matrix.
#'Default is \code{TRUE} (diagonal variance).
#'
#'@param verbose logical flag indicating wether partition info is
#'written in the console at each MCMC iteration.
#'
#'@param distrib the distribution used for the clustering. Current possibilities are
#'\code{"gaussian"}, \code{"skewnorm"} and \code{"skewt"}.
#'
#'@param ncores number of cores to use.
#'
#'@param type_connec The type of connection between the processors. Supported
#'cluster types are \code{"SOCK"}, \code{"FORK"}, \code{"MPI"}, and
#'\code{"NWS"}. See also \code{\link[parallel:makeCluster]{makeCluster}}.
#'
#'@param informPrior an optional informative prior such as the approximation computed
#'by \code{summary.DPMMclust}.
#'
#'@param ... additional arguments to be passed to \code{\link{plot_DPM}}.
#'Only used if \code{doPlot} is \code{TRUE}.
#'
#'@return a object of class \code{DPMclust} with the following attributes:
#'  \itemize{
#'      \item{\code{mcmc_partitions}:}{ a list of length \code{N}. Each
#'       element \code{mcmc_partitions[n]} is a vector of length
#'       \code{n} giving the partition of the \code{n} observations.}
#'      \item{\code{alpha}:}{a vector of length \code{N}. \code{cost[j]} is the cost
#' associated to partition \code{c[[j]]}}
#'       \item{\code{U_SS_list}:}{a list of length \code{N} containing the lists of
#'       sufficient statistics for all the mixture components at each MCMC iteration}
#'      \item{\code{weights_list}:}{a list of length \code{N} containing the weights of each
#'      mixture component for each MCMC iterations}
#'      \item{\code{logposterior_list}:}{a list of length \code{N} containing the logposterior values
#'       at each MCMC iterations}
#'      \item{\code{data}:}{the data matrix \code{d x n} with \code{d} dimensions in rows
#'and \code{n} observations in columns}
#'      \item{\code{nb_mcmcit}:}{ the number of MCMC itertations}
#'      \item{\code{clust_distrib}:}{the parametric distribution of the mixture component}
#'      \item{\code{hyperG0}:}{the prior on the cluster location}
#'  }
#'
#'@author Boris Hejblum
#'
#'@references Hejblum BP, Alkhassim C, Gottardo R, Caron F, Thiebaut R, Sequential Dirichlet
#'Process Mixtures of Multivariate Skew t-distributions for Model-based Clustering
#'of Flow Cytometry Data, in preparation.
#'
#'@export
#'
#'@examples
#' rm(list=ls())
#' library(ggplot2)
#' library(truncnorm)
#'
#' #Number of data
#' n <- 2000
#' set.seed(123)
#' set.seed(4321)
#'
#'
#' d <- 2
#' ncl <- 4
#'
#' # Sample data
#'
#' sdev <- array(dim=c(d,d,ncl))
#'
#' xi <- matrix(nrow=d, ncol=ncl, c(-1.5, 1.5, 1.5, 1.5, 2, -2.5, -2.5, -3))
#' psi <- matrix(nrow=d, ncol=4, c(0.3, -0.7, -0.8, 0, 0.3, -0.7, 0.2, 0.9))
#' nu <- c(100,25,8,5)
#' p <- c(0.15, 0.05, 0.5, 0.3) # frequence des clusters
#' sdev[, ,1] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.3))
#' sdev[, ,2] <- matrix(nrow=d, ncol=d, c(0.1, 0, 0, 0.3))
#' sdev[, ,3] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.2))
#' sdev[, ,4] <- .3*diag(2)
#'
#'
#' c <- rep(0,n)
#' w <- rep(1,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  w[k] <- rgamma(1, shape=nu[c[k]]/2, rate=nu[c[k]]/2)
#'  z[,k] <- xi[, c[k]] + psi[, c[k]]*rtruncnorm(n=1, a=0, b=Inf, mean=0, sd=1/sqrt(w[k])) +
#'                (sdev[, , c[k]]/sqrt(w[k]))%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
#'  #cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'
#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["b_xi"]] <- rowMeans(z)
#' hyperG0[["b_psi"]] <- rep(0,d)
#' hyperG0[["kappa"]] <- 0.001
#' hyperG0[["D_xi"]] <- 100
#' hyperG0[["D_psi"]] <- 100
#' hyperG0[["nu"]] <- d+1
#' hyperG0[["lambda"]] <- diag(apply(z,MARGIN=1, FUN=var))/3
#'
#'  # hyperprior on the Scale parameter of DPM
#'  a <- 0.0001
#'  b <- 0.0001
#'
#'
#'
#'  ## Data
#'  ########
#'  library(ggplot2)
#'  p <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,]), aes(x=X, y=Y))
#'        + geom_point()
#'        #+ ggtitle("Simple example in 2d data")
#'        +xlab("D1")
#'        +ylab("D2")
#'        +theme_bw())
#'  p #pdf(height=8.5, width=8.5)
#'
#'  c2plot <- factor(c)
#'  levels(c2plot) <- c("4", "1", "3", "2")
#'  pp <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,], "Cluster"=as.character(c2plot)))
#'        + geom_point(aes(x=X, y=Y, colour=Cluster, fill=Cluster))
#'        #+ ggtitle("Slightly overlapping skew-normal simulation\n")
#'        + xlab("D1")
#'        + ylab("D2")
#'        + theme_bw()
#'        + scale_colour_discrete(guide=guide_legend(override.aes = list(size = 6, shape=22))))
#'  pp #pdf(height=7, width=7.5)
#'
#'\dontrun{
#'  MCMCsample_st <- DPMpost(data=z, hyperG0=hyperG0, N=2000,
#'                           distrib="skewt",
#'                           gg.add=list(theme_bw(),
#'                           guides(shape=guide_legend(override.aes = list(fill="grey45"))))
#'  )
#'  s <- summary(MCMCsample_st, burnin = 1500, thin=5, lossFn = "Binder")
#'  s
#'  plot(s)
#'  #plot(s, hm=TRUE) #pdf(height=8.5, width=10.5) #png(height=700, width=720)
#'}
#'
#'
#'
#'
#'
#'
DPMpost <- function (data, hyperG0, a=0.0001, b=0.0001, N, doPlot=TRUE,
                     nbclust_init=30, plotevery=floor(N/10),
                     diagVar=TRUE, verbose=TRUE,
                     distrib=c("gaussian", "skewnorm", "skewt"),
                     ncores = 1,
                     type_connec = "SOCK",
                     informPrior=NULL,
                     ...
){

  if(ncores>1){
    if(ncores < parallel::detectCores()){
      stop("Number of requested cores is higher than what is available")
    }
  }


  if(ncores<2){
    if(is.null(informPrior)){
      res <- switch(distrib,
                    "gaussian"=DPMGibbsN(data, hyperG0, a, b, N, doPlot, nbclust_init, plotevery, diagVar,
                                         verbose, ...),
                    "skewnorm"=DPMGibbsSkewN(data, hyperG0, a, b, N, doPlot, nbclust_init, plotevery, diagVar,
                                             verbose, ...),
                    "skewt"=DPMGibbsSkewT(data, hyperG0, a, b, N, doPlot, nbclust_init, plotevery, diagVar,
                                          verbose, ...)
      )
    }else{
      res <- switch(distrib,
                    "gaussian"=DPMGibbsN_SeqPrior(data, informPrior, hyperG0, N,
                                                  nbclust_init, doPlot=doPlot, plotevery=plotevery,
                                                  diagVar=diagVar, verbose=verbose, ...),
                    "skewnorm"=stop("Skew normal ditributions with informative prior is not implemented yet.\n",
                                    "Contact the maintainer if you would like to see this feature implemented.\n",
                                    "In the meantime, try the skew t distribution with 'skewt' which is a generalization ",
                                    "of the skew normal distribution."),
                    "skewt"=DPMGibbsSkewT_SeqPrior(data, informPrior, hyperG0, N, nbclust_init,
                                                   doPlot=doPlot, plotevery=plotevery, diagVar=diagVar,
                                                   verbose=verbose, ...)
      )
    }
  }else{
    if(is.null(informPrior)){
      if(distrib=="skewnorm"){
        warning("Parallel implementation with skew normal ditributions is not available yet.\n",
                "Contact the maintainer if you would like to see this feature implemented.\n",
                "In the meantime, the non-parallel implementation is being run instead")
      }
      res <- switch(distrib,
                    "gaussian"=DPMGibbsN_parallel(ncores, type_connec, data, hyperG0, a, b, N,
                                                  doPlot, nbclust_init, plotevery, diagVar,
                                                  verbose, ...),
                    "skewnorm"=DPMGibbsSkewN(data, hyperG0, a, b, N, doPlot, nbclust_init, plotevery, diagVar,
                                             verbose, ...),
                    "skewt"=DPMGibbsSkewT_parallel(ncores, type_connec, data, hyperG0, a, b, N,
                                                   doPlot, nbclust_init, plotevery, diagVar,
                                                   verbose, ...)
      )
    }else{
      warning("Parallel implementation with an informative prior for gaussian ditributions is not available yet.\n",
              "Contact the maintainer if you would like to see this feature implemented.\n",
              "In the meantime, the non-parallel implementation is being run instead.")
      res <- switch(distrib,
                    "gaussian"=DPMGibbsN_SeqPrior(data, informPrior, hyperG0, N,
                                                  nbclust_init, doPlot=doPlot, plotevery=plotevery,
                                                  diagVar=diagVar, verbose=verbose, ...),
                    "skewnorm"=stop("Skew normal ditributions with informative prior is not implemented yet.\n",
                                    "Contact the maintainer if you would like to see this feature implemented.\n",
                                    "In the meantime, try the skew t distribution with 'skewt' which is a generalization ",
                                    "of the skew normal distribution."),
                    "skewt"=DPMGibbsSkewT_SeqPrior_parallel(ncores, type_connec, data, informPrior,
                                                            hyperG0, N, nbclust_init, doPlot=doPlot,
                                                            plotevery=plotevery, diagVar=diagVar,
                                                            verbose=verbose, ...)
      )
    }
  }

  return(res)
}
