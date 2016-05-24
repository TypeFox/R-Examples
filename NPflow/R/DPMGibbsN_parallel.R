#'Slice Sampling of the Dirichlet Process Mixture Model
#'with a prior on alpha
#'
#'@param Ncpus the number of processors available
#'
#'@param type_connec The type of connection between the processors. Supported
#'cluster types are \code{"SOCK"}, \code{"FORK"}, \code{"MPI"}, and
#'\code{"NWS"}. See also \code{\link[parallel:makeCluster]{makeCluster}}.
#'
#'@param z data matrix \code{d x n} with \code{d} dimensions in rows
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
#'@param plotevery an integer indicating the interval between plotted iterations when \code{doPlot}
#'is \code{TRUE}.
#'
#'@param nbclust_init number of clusters at initialisation.
#'Default to 30 (or less if there are less than 30 observations).
#'
#'@param diagVar logical flag indicating wether the variance of each cluster is
#'estimated as a diagonal matrix, or as a full matrix.
#'Default is \code{TRUE} (diagonal variance).
#'
#'@param use_variance_hyperprior logical flag indicating whether a hyperprior is added 
#'for the variance parameter. Default is \code{TRUE} which decrease the impact of the variance prior
#'on the posterior. \code{FALSE} is useful for using an informative prior.
#'
#'@param verbose logical flag indicating wether partition info is
#'written in the console at each MCMC iteration.
#'
#'@param monitorfile
#'a writable \link{connections} or a character string naming a file to write into,
#'to monitor the progress of the analysis.
#'Default is \code{""} which is no monitoring.  See Details.
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
#'       \item{\code{listU_mu}:}{a list of length \code{N} containing the matrices of
#'       mean vectors for all the mixture components at each MCMC iteration}
#'       \item{\code{listU_Sigma}:}{a list of length \code{N} containing the arrays of
#'       covariances matrices for all the mixture components at each MCMC iteration}
#'       \item{\code{U_SS_list}:}{a list of length \code{N} containing the lists of
#'       sufficient statistics for all the mixture components at each MCMC iteration}
#'      \item{\code{weights_list}:}{a list of length \code{N} containing the logposterior values
#'       at each MCMC iterations}
#'      \item{\code{logposterior_list}:}{a list of length \code{N} containing the logposterior values
#'       at each MCMC iterations}
#'      \item{\code{data}:}{the data matrix \code{d x n} with \code{d} dimensions in rows
#'and \code{n} observations in columns}
#'      \item{\code{nb_mcmcit}:}{ the number of MCMC itertations}
#'      \item{\code{clust_distrib}:}{the parametric distribution of the mixture component - \code{"gaussian"}}
#'      \item{\code{hyperG0}:}{the prior on the cluster location}
#'  }
#'
#'@author Boris Hejblum
#'
#'@seealso \code{\link{DPMGibbsN}}
#'
#'@export
#'
#'@examples
#'
#' # Scaling up: ----
#' rm(list=ls())
#' #Number of data
#' n <- 2000
#' set.seed(1234)
#'
#' # Sample data
#' d <- 3
#' nclust <- 5
#' m <- matrix(nrow=d, ncol=nclust, runif(d*nclust)*8)
#' # p: cluster probabilities
#' p <- runif(nclust)
#' p <- p/sum(p)
#'
#' # Covariance matrix of the clusters
#' sdev <- array(dim=c(d, d, nclust))
#' for (j in 1:nclust){
#'     sdev[, ,j] <- matrix(NA, nrow=d, ncol=d)
#'     diag(sdev[, ,j]) <- abs(rnorm(n=d, mean=0.3, sd=0.1))
#'     sdev[, ,j][lower.tri(sdev[, ,j], diag = FALSE)] <- rnorm(n=d*(d-1)/2,
#'     mean=0, sd=0.05)
#'     sdev[, ,j][upper.tri(sdev[, ,j], diag = FALSE)] <- (sdev[, ,j][
#'                                                         lower.tri(sdev[, ,j], diag = FALSE)])
#' }
#' c <- rep(0,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'     c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'     z[,k] <- m[, c[k]] + sdev[, , c[k]]%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
#'     #cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'
#' # hyperprior on the Scale parameter of DPM
#' a <- 0.001
#' b <- 0.001
#'
#' # Number of iterations
#' N <- 25
#'
#' # do some plots
#' doPlot <- TRUE
#'
#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["mu"]] <- rep(0, d)
#' hyperG0[["kappa"]] <- 0.01
#' hyperG0[["nu"]] <- d + 2
#' hyperG0[["lambda"]] <- diag(d)/10
#'
#'
#' nbclust_init <- 30
#' 
#'\dontrun{
#'  library(doParallel)
#'  MCMCsample <- DPMGibbsN_parallel(Ncpus=2, type_connec="FORK", z, hyperG0, a, b, N=1000, doPlot, 
#'                                   nbclust_init, plotevery=100, gg.add=list(theme_bw(),
#'                                   guides(shape=guide_legend(override.aes = list(fill="grey45")))),
#'                                   diagVar=FALSE)
#'}
#'
#'
DPMGibbsN_parallel <- function (Ncpus, type_connec,
                                z, hyperG0, a=0.0001, b=0.0001, N, doPlot=TRUE,
                                nbclust_init=30, plotevery=N/10,
                                diagVar=TRUE, use_variance_hyperprior=TRUE, verbose=TRUE, monitorfile="",
                                ...){
  
  dpmclus <- NULL
  
  if(!requireNamespace("doParallel", quietly=TRUE)){
    stop("Package 'doParallel' is not available.\n  -> Try running 'install.packages(\"doParallel\")'\n   or use non parallel version of the function: 'DPMGibbsN'")
  }else{
    requireNamespace("doParallel", quietly=TRUE)
    
    # declare the cores
    cl <- parallel::makeCluster(Ncpus, type = type_connec, outfile=monitorfile)
    doParallel::registerDoParallel(cl)
    
    if(doPlot){requireNamespace("ggplot2", quietly=TRUE)}
    
    p <- dim(z)[1]
    n <- dim(z)[2]
    U_xi <- matrix(0, nrow=p, ncol=n)
    U_psi <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    U_df = rep(10,n)
    U_B = array(0, dim=c(2, 2, n))
    U_nu <- rep(p,n)
    
    par_ind <- list()
    temp_ind <- 0
    if(Ncpus>1){
      nb_simult <- floor(n%/%(Ncpus))
      for(i in 1:(Ncpus-1)){
        par_ind[[i]] <- temp_ind + 1:nb_simult
        temp_ind <- temp_ind + nb_simult
      }
      par_ind[[Ncpus]] <- (temp_ind+1):n
    }
    else{
      cat("Only 1 core specified\n=> non-parallel version of the algorithm would be more efficient",
          file=monitorfile, append = TRUE)
      nb_simult <- n
      par_ind[[Ncpus]] <- (temp_ind+1):n
    }
    
    
    p <- nrow(z)
    n <- ncol(z)
    U_mu <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    listU_mu<-list()
    listU_Sigma<-list()
    
    par_ind <- list()
    temp_ind <- 0
    if(Ncpus>1){
      nb_simult <- floor(n%/%(Ncpus))
      for(i in 1:(Ncpus-1)){
        par_ind[[i]] <- temp_ind + 1:nb_simult
        temp_ind <- temp_ind + nb_simult
      }
      par_ind[[Ncpus]] <- (temp_ind+1):n
    }
    else{
      cat("Only 1 core specified\n=> non-parallel version of the algorithm would be more efficient",
          file=monitorfile, append = TRUE)
      nb_simult <- n
      par_ind[[Ncpus]] <- (temp_ind+1):n
    }
    
    # U_SS is a list where each U_SS[[k]] contains the sufficient
    # statistics associated to cluster k
    U_SS <- list()
    
    #store U_SS :
    U_SS_list <- list()
    #store clustering :
    c_list <- list()
    #store sliced weights
    weights_list <- list()
    
    #store log posterior probability
    logposterior_list <- list()
    
    
    m <- numeric(n) # number of obs in each clusters
    c <-numeric(n) # initial number of clusters
    
    
    # Initialisation----
    # each observation is assigned to a different cluster
    # or to 1 of the 50 initial clusters if there are more than
    # 50 observations
    
    i <- 1
    if(n<nbclust_init){
      for (k in 1:n){
        c[k] <- k
        #cat("cluster ", k, ":\n")
        U_SS[[k]] <- update_SS(z=z[, k, drop=FALSE], S=hyperG0, hyperprior = NULL)
        NiW <- rNiW(U_SS[[k]], diagVar=diagVar)
        
        U_mu[, k] <- NiW[["mu"]]
        U_SS[[k]][["mu"]] <- NiW[["mu"]]
        
        U_Sigma[, , k] <- NiW[["S"]]
        U_SS[[k]][["S"]] <- NiW[["S"]]
        
        m[k] <- m[k]+1
        U_SS[[k]][["weight"]] <- 1/n
      }
    } else{
      c <- sample(x=1:nbclust_init, size=n, replace=TRUE)
      for (k in unique(c)){
        obs_k <- which(c==k)
        #cat("cluster ", k, ":\n")
        U_SS[[k]] <- update_SS(z=z[, obs_k,drop=FALSE], S=hyperG0)
        NiW <- rNiW(U_SS[[k]], diagVar=diagVar)
        
        U_mu[, k] <- NiW[["mu"]]
        U_SS[[k]][["mu"]] <- NiW[["mu"]]
        
        U_Sigma[, , k] <- NiW[["S"]]
        U_SS[[k]][["S"]] <- NiW[["S"]]
        
        m[k] <- length(obs_k)
        U_SS[[k]][["weight"]] <- m[k]/n
      }
    }
    listU_mu[[i]]<-U_mu
    listU_Sigma[[i]]<-U_Sigma
    
    alpha <- c(log(n))
    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    weights_list[[i]] <- numeric(length(m))
    weights_list[[i]][unique(c)] <- table(c)/length(c)
    
    logposterior_list[[i]] <- logposterior_DPMG(z, mu=U_mu, Sigma=U_Sigma,
                                                hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b, diagVar)
    
    if(verbose){
      cat(i, "/", N, " samplings:\n", sep="")
      cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
      cl2print <- unique(c)
      cat("  ",length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
    }
    
    if(doPlot){
      plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma,
               m=m, c=c, i=i, alpha=alpha[length(alpha)], U_SS=U_SS, ...)
    }else if(verbose){
      cl2print <- unique(c)
      cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
    }
    
    
    
    # MCMC iterations
    
    for(i in 2:N){
      nbClust <- length(unique(c))
      alpha <- c(alpha,
                 sample_alpha(alpha_old=alpha[i-1], n=n,
                              K=nbClust, a=a, b=b)
      )
      slice <- sliceSampler_N_parallel(Ncpus=Ncpus, c=c, m=m, alpha=alpha[i],
                                       z=z, hyperG0=hyperG0,
                                       U_mu=U_mu, U_Sigma=U_Sigma, diagVar=diagVar,
                                       parallel_index=par_ind)
      
      m <- slice[["m"]]
      c <- slice[["c"]]
      weights_list[[i]] <- slice[["weights"]]
      U_mu<-slice[["U_mu"]]
      U_Sigma<-slice[["U_Sigma"]]
      
      # Update cluster locations
      fullCl <- which(m!=0)
      for(j in fullCl){
        obs_j <- which(c==j)
        #cat("cluster ", j, ":\n")
        if(use_variance_hyperprior){
          U_SS[[j]] <- update_SS(z=z[, obs_j,drop=FALSE], S=hyperG0, hyperprior = list("Sigma"=U_Sigma[,,j]))
        }else{
          U_SS[[j]] <- update_SS(z=z[, obs_j,drop=FALSE], S=hyperG0)
        }
        NiW <- rNiW(U_SS[[j]], diagVar=diagVar)
        
        U_mu[, j] <- NiW[["mu"]]
        U_SS[[j]][["mu"]] <- NiW[["mu"]]
        
        U_Sigma[, , j] <- NiW[["S"]]
        U_SS[[j]][["S"]] <- NiW[["S"]]
        
        U_SS[[j]][["weight"]] <- weights_list[[i]][j]
        #cat("sampled S =", NiW[["S"]], "\n\n\n")
      }
      
      listU_mu[[i]]<-U_mu
      listU_Sigma[[i]]<-U_Sigma
      U_SS_list[[i]] <- U_SS[which(m!=0)]
      c_list[[i]] <- c
      logposterior_list[[i]] <- logposterior_DPMG(z, mu=U_mu, Sigma=U_Sigma,
                                                  hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b, diagVar)
      
      if(verbose){
        cat(i, "/", N, " samplings:\n", sep="")
        cat("  logposterior = ", sum(logposterior_list[[i]]) , "\n", sep="")
        cl2print <- unique(c)
        cat("  ",length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
      }
      
      if(doPlot && i/plotevery==floor(i/plotevery)){
        plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma,
                 m=m, c=c, i=i, alpha=alpha[length(alpha)], U_SS=U_SS, ...)
      }
      
    }
    
    parallel::stopCluster(cl)
    
    dpmclus <- list("mcmc_partitions" = c_list,
                    "alpha"=alpha,
                    "listU_mu"=listU_mu,
                    "listU_Sigma"=listU_Sigma,
                    "U_SS_list"=U_SS_list,
                    "weights_list"=weights_list,
                    "logposterior_list"=logposterior_list,
                    "data"=z,
                    "nb_mcmcit"=N,
                    "clust_distrib"="gaussian",
                    "hyperG0"=hyperG0)
    class(dpmclus) <- "DPMMclust"
    
  }
  return(dpmclus)
}






