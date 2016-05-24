#'Parallel Implementation of Slice Sampling of Dirichlet Process Mixture of skew Normals
#'
#'If the \code{monitorfile} argument is a character string naming a file to
#'write into, in the case of a new file that does not exist yet, such a new
#'file will be created. A line is written at each MCMC iteration.
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
#'      \item{\code{mcmc_partitions}:}{a list of length \code{N}. Each
#'       element \code{mcmc_partitions[n]} is a vector of length
#'       \code{n} giving the partition of the \code{n} observations.}
#'      \item{\code{alpha}:}{a vector of length \code{N}. \code{cost[j]} is the cost
#' associated to partition \code{c[[j]]}}
#'       \item{\code{U_SS_list}:}{a list of length \code{N} containing the lists of
#'       sufficient statistics for all the mixture components at each MCMC iteration}
#'      \item{\code{weights_list}:}{}
#'      \item{\code{logposterior_list}:}{a list of length \code{N} containing the logposterior values
#'       at each MCMC iterations}
#'      \item{\code{data}:}{the data matrix \code{d x n} with \code{d} dimensions in rows
#'and \code{n} observations in columns}
#'      \item{\code{nb_mcmcit}:}{the number of MCMC itertations}
#'      \item{\code{clust_distrib}:}{the parametric distribution of the mixture component - \code{"skewnorm"}}
#'      \item{\code{hyperG0}:}{the prior on the cluster location}
#'  }
#'
#'@author Boris Hejblum
#'
#'@export
#'
#'@references Hejblum BP, Alkhassim C, Gottardo R, Caron F, Thiebaut R, Sequential Dirichlet
#'Process Mixtures of Multivariate Skew t-distributions for Model-based Clustering
#'of Flow Cytometry Data, in preparation.
#'
#'@examples
#' rm(list=ls())
#' #Number of data
#' n <- 2000
#' set.seed(1234)
#'
#' d <- 4
#' ncl <- 5
#'
#' # Sample data
#'
#' sdev <- array(dim=c(d,d,ncl))
#'
#' xi <- matrix(nrow=d, ncol=ncl, c(runif(n=d*ncl,min=0,max=3)))
#' psi <- matrix(nrow=d, ncol=ncl, c(runif(n=d*ncl,min=-1,max=1)))
#' p <- runif(n=ncl)
#' p <- p/sum(p)
#' sdev0 <- diag(runif(n=d, min=0.05, max=0.6))
#' for (j in 1:ncl){
#'      sdev[, ,j] <- invwishrnd(n = d+2, lambda = sdev0)
#' }
#'
#'
#' c <- rep(0,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- xi[, c[k]] + psi[, c[k]]*abs(rnorm(1)) + sdev[, , c[k]]%*%matrix(rnorm(d, mean = 0,
#'                                                                         sd = 1), nrow=d, ncol=1)
#'  #cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'
#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["b_xi"]] <- rep(0,d)
#' hyperG0[["b_psi"]] <- rep(0,d)
#' hyperG0[["kappa"]] <- 0.001
#' hyperG0[["D_xi"]] <- 100
#' hyperG0[["D_psi"]] <- 100
#' hyperG0[["nu"]] <- d + 1
#' hyperG0[["lambda"]] <- diag(d)/10
#'
#'  # hyperprior on the Scale parameter of DPM
#'  a <- 0.0001
#'  b <- 0.0001
#'
#'  # do some plots
#'  doPlot <- TRUE
#'  nbclust_init <- 30
#'
#'
#'  z <- z*200
#'  ## Data
#'  ########
#' library(ggplot2)
#'  p <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,]), aes(x=X, y=Y))
#'        + geom_point()
#'        + ggtitle("Simple example in 2d data")
#'        +xlab("D1")
#'        +ylab("D2")
#'        +theme_bw())
#'  p
#'
#'
#'  ## alpha priors plots
#'  #####################
#'  prioralpha <- data.frame("alpha"=rgamma(n=5000, shape=a, scale=1/b),
#'                          "distribution" =factor(rep("prior",5000),
#'                          levels=c("prior", "posterior")))
#'  p <- (ggplot(prioralpha, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..),
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="red")
#'        + ggtitle(paste("Prior distribution on alpha: Gamma(", a,
#'                  ",", b, ")\n", sep=""))
#'       )
#'  p
#'
#'
#'
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  ##############################################
#'  \dontrun{
#'  MCMCsample_sn_par <- DPMGibbsSkewN_parallel(Ncpus=parallel::detectCores()-1,
#'                                              type_connec="SOCK", z, hyperG0,
#'                                              a, b, N=5000, doPlot, nbclust_init,
#'                                              plotevery=25, gg.add=list(theme_bw(),
#'                                 guides(shape=guide_legend(override.aes = list(fill="grey45")))))
#'  plot_ConvDPM(MCMCsample_sn_par, from=2)
#'
#'  postalpha <- data.frame("alpha"=MCMCsample_sn_par$alpha[50:500],
#'                          "distribution" = factor(rep("posterior",500-49),
#'                          levels=c("prior", "posterior")))
#'  p <- (ggplot(postalpha, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..), binwidth=.1,
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="blue")
#'        + ggtitle("Posterior distribution of alpha\n")
#'        # Ignore NA values for mean
#'        # Overlay with transparent density plot
#'        + geom_vline(aes(xintercept=mean(alpha, na.rm=T)),
#'                     color="red", linetype="dashed", size=1)
#'      )
#'  p
#'
#'  p <- (ggplot(drop=FALSE, alpha=.6)
#'        + geom_density(aes(x=alpha, fill=distribution),
#'                       color=NA, alpha=.6,
#'                       data=prioralpha)
#'        + geom_density(aes(x=alpha, fill=distribution),
#'                       color=NA, alpha=.6,
#'                       data=postalpha)
#'        + ggtitle("Prior and posterior distributions of alpha\n")
#'        + scale_fill_discrete(drop=FALSE)
#'      )
#'  p
#'}
#'
#'# k-means
#'#########
#'
#'  plot(x=z[1,], y=z[2,], col=kmeans(t(z), centers=4)$cluster,
#'       xlab = "d = 1", ylab= "d = 2", main="k-means with K=4 clusters")
#'
#'  KM <- kmeans(t(z), centers=4)
#'  dataKM <- data.frame("X"=z[1,], "Y"=z[2,],
#'                     "Cluster"=as.character(KM$cluster))
#'  dataCenters <- data.frame("X"=KM$centers[,1],
#'                            "Y"=KM$centers[,2],
#'                            "Cluster"=rownames(KM$centers))
#'
#'  p <- (ggplot(dataKM)
#'        + geom_point(aes(x=X, y=Y, col=Cluster))
#'        + geom_point(aes(x=X, y=Y, fill=Cluster, order=Cluster),
#'                     data=dataCenters, shape=22, size=5)
#'        + scale_colour_discrete(name="Cluster")
#'        + ggtitle("K-means with K=4 clusters\n"))
#'  p
#'
#'
#'
#'
DPMGibbsSkewN_parallel <- function (Ncpus, type_connec,
                                    z, hyperG0, a=0.0001, b=0.0001, N, doPlot=FALSE,
                                    nbclust_init=30, plotevery=N/10,
                                    diagVar=TRUE, use_variance_hyperprior=TRUE, verbose=FALSE,
                                    monitorfile="",
                                    ...){

  dpmclus <- NULL

  if(!requireNamespace("doParallel", quietly=TRUE)){
    stop("Package 'doParallel' is not available.\n  -> Try running 'install.packages(\"doParallel\")'\n   or use non parallel version of the function: 'DPMGibbsN'")
  }else{
    requireNamespace("doParallel", quietly=TRUE)

    # declare the cores
    cl <- parallel::makeCluster(Ncpus, type = type_connec, outfile=monitorfile)
    doParallel::registerDoParallel(cl)

    if(doPlot){requireNamespace("ggplot2", quietly = TRUE)}

    p <- dim(z)[1]
    n <- dim(z)[2]
    U_xi <- matrix(0, nrow=p, ncol=n)
    U_psi <- matrix(0, nrow=p, ncol=n)
    U_Sigma = array(0, dim=c(p, p, n))
    U_B = array(0, dim=c(2, 2, n))

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
    c <- numeric(n) # cluster label of ech observation
    ltn <- rtruncnorm(n, a=0, b=Inf, mean=0, sd=1) # latent truncated normal

    # Initialisation----
    # each observation is assigned to a different cluster
    # or to 1 of the 50 initial clusters if there are more than
    # 50 observations

    i <- 1

    if(ncol(z)<nbclust_init){
      for (k in 1:n){
        c[k] <- k
        U_SS[[k]] <- update_SSsn(z=z[, k], S=hyperG0, ltn=ltn[k],
                                 hyperprior = NULL)
        NNiW <- rNNiW(U_SS[[k]], diagVar)
        U_xi[, k] <- NNiW[["xi"]]
        U_SS[[k]][["xi"]] <- NNiW[["xi"]]
        U_psi[, k] <- NNiW[["psi"]]
        U_SS[[k]][["psi"]] <- NNiW[["psi"]]
        U_Sigma[, , k] <- NNiW[["S"]]
        U_SS[[k]][["S"]] <- NNiW[["S"]]
        U_B[, ,k] <- U_SS[[k]][["B"]]
        m[k] <- m[k]+1
      }
    } else{
      c <- sample(x=1:nbclust_init, size=n, replace=TRUE)
      for (k in unique(c)){
        obs_k <- which(c==k)
        U_SS[[k]] <- update_SSsn(z=z[, obs_k], S=hyperG0, ltn=ltn[obs_k],
                                 hyperprior = NULL)
        NNiW <- rNNiW(U_SS[[k]], diagVar)
        U_xi[, k] <- NNiW[["xi"]]
        U_SS[[k]][["xi"]] <- NNiW[["xi"]]
        U_psi[, k] <- NNiW[["psi"]]
        U_SS[[k]][["psi"]] <- NNiW[["psi"]]
        U_Sigma[, , k] <- NNiW[["S"]]
        U_SS[[k]][["S"]] <- NNiW[["S"]]
        U_B[, ,k] <- U_SS[[k]][["B"]]
        m[k] <- length(obs_k)
      }
    }



    alpha <- c(log(n))


    U_SS_list[[i]] <- U_SS
    c_list[[i]] <- c
    weights_list[[1]] <- numeric(length(m))
    weights_list[[1]][unique(c)] <- table(c)/length(c)

    logposterior_list[[i]] <- logposterior_DPMSN(z, xi=U_xi, psi=U_psi, Sigma=U_Sigma, B=U_B,
                                                 hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b, diagVar)

    if(doPlot){
      plot_DPMsn(z=z, c=c, i=i, alpha=alpha[i], U_SS=U_SS_list[[i]], ellipses=TRUE, ...)
    }
    if(verbose){
      cat(i, "/", N, " samplings:\n", sep="",
          file=monitorfile, append = TRUE)
      cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="",
          file=monitorfile, append = TRUE)
      cl2print <- unique(c)
      cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n",
          file=monitorfile, append = TRUE)
    }


    if(N>1){
      for(i in 2:N){
        nbClust <- length(unique(c))

        alpha <- c(alpha,
                   sample_alpha(alpha_old=alpha[i-1], n=n,
                                K=nbClust, a=a, b=b)
        )

        slice <- sliceSampler_SkewN_parallel(Ncpus=Ncpus,
                                             c=c, m=m,
                                             alpha=alpha[i],
                                             z=z, hyperG0=hyperG0,
                                             U_xi=U_xi,
                                             U_psi=U_psi,
                                             U_Sigma=U_Sigma,
                                             diagVar=diagVar,
                                             parallel_index=par_ind)
        m <- slice[["m"]]
        c <- slice[["c"]]
        weights_list[[i]] <- slice[["weights"]]
        ltn <- slice[["latentTrunc"]]



        # Update cluster locations
        fullCl <- which(m!=0)
        for(j in fullCl){
          obs_j <- which(c==j)
          if(use_variance_hyperprior){
            U_SS[[j]] <- update_SSsn(z=z[, obs_j, drop=FALSE], S=hyperG0,
                                     ltn=ltn[obs_j],
                                     hyperprior = list("Sigma"=U_Sigma[,,j])
            )
          }else{
            U_SS[[j]] <- update_SSsn(z=z[, obs_j, drop=FALSE], S=hyperG0,
                                     ltn=ltn[obs_j]
            )
          }
          
          NNiW <- rNNiW(U_SS[[j]], diagVar)
          U_xi[, j] <- NNiW[["xi"]]
          U_SS[[j]][["xi"]] <- NNiW[["xi"]]
          U_psi[, j] <- NNiW[["psi"]]
          U_SS[[j]][["psi"]] <- NNiW[["psi"]]
          U_Sigma[, , j] <- NNiW[["S"]]
          U_SS[[j]][["S"]] <- NNiW[["S"]]
          U_B[, ,j] <- U_SS[[j]][["B"]]
        }


        U_SS_list[[i]] <- U_SS[which(m!=0)]
        c_list[[i]] <- c

        logposterior_list[[i]] <- logposterior_DPMSN(z, xi=U_xi, psi=U_psi, Sigma=U_Sigma, B=U_B,
                                                     hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b, diagVar)

        if(doPlot && i/plotevery==floor(i/plotevery)){
          plot_DPMsn(z=z, c=c, i=i, alpha=alpha[i], U_SS=U_SS_list[[i]], ellipses=TRUE, ...)
        }
        if(verbose){
          cat(i, "/", N, " samplings:\n", sep="",
              file=monitorfile, append = TRUE)
          cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="",
              file=monitorfile, append = TRUE)
          cl2print <- unique(c)
          cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n",
              file=monitorfile, append = TRUE)
        }
      }
    }

    parallel::stopCluster(cl)

    dpmclus <- list("mcmc_partitions" = c_list,
                    "alpha"=alpha,
                    "U_SS_list"=U_SS_list,
                    "weights_list"=weights_list,
                    "logposterior_list"=logposterior_list,
                    "data"=z,
                    "nb_mcmcit"=N,
                    "clust_distrib"="skewnorm",
                    "hyperG0"=hyperG0)
    class(dpmclus) <- "DPMMclust"
  }
  return(dpmclus)
}






