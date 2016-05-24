#'Slice Sampling of Dirichlet Process Mixture of skew Students's t-distibutions
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
#'@param verbose logical flag indicating whether partition info is
#'written in the console at each MCMC iteration.
#'
#'@param ... additional arguments to be passed to \code{\link{plot_DPMst}}.
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
#'      \item{\code{nb_mcmcit}:}{the number of MCMC itertations}
#'      \item{\code{clust_distrib}:}{the parametric distribution of the mixture component - \code{"skewt"}}
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
#'
#' #Number of data
#' n <- 2000
#' set.seed(4321)
#'
#'
#' d <- 2
#' ncl <- 4
#'
#' # Sample data
#' library(truncnorm)
#'
#' sdev <- array(dim=c(d,d,ncl))
#'
#' #xi <- matrix(nrow=d, ncol=ncl, c(-1.5, 1.5, 1.5, 1.5, 2, -2.5, -2.5, -3))
#' #xi <- matrix(nrow=d, ncol=ncl, c(-0.5, 0, 0.5, 0, 0.5, -1, -1, 1))
#' xi <- matrix(nrow=d, ncol=ncl, c(-0.2, 0.5, 2.4, 0.4, 0.6, -1.3, -0.9, -2.7))
#' psi <- matrix(nrow=d, ncol=4, c(0.3, -0.7, -0.8, 0, 0.3, -0.7, 0.2, 0.9))
#' nu <- c(100,25,8,5)
#' p <- c(0.15, 0.05, 0.5, 0.3) # frequence des clusters
#' sdev[, ,1] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.3))
#' sdev[, ,2] <- matrix(nrow=d, ncol=d, c(0.1, 0, 0, 0.3))
#' sdev[, ,3] <- matrix(nrow=d, ncol=d, c(0.3, 0.15, 0.15, 0.3))
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
#'                 (sdev[, , c[k]]/sqrt(w[k]))%*%matrix(rnorm(d, mean = 0, sd = 1), nrow=d, ncol=1)
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
#'  # do some plots
#'  doPlot <- TRUE
#'  nbclust_init <- 30
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
#'\dontrun{
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  ##############################################
#'  MCMCsample_st <- DPMGibbsSkewT(z, hyperG0, a, b, N=2000,
#'                                 doPlot, nbclust_init, plotevery=100,
#'                                 gg.add=list(theme_bw(),
#'                                  guides(shape=guide_legend(override.aes = list(fill="grey45")))),
#'                                diagVar=FALSE)
#'  s <- summary(MCMCsample_st, burnin = 900, thin=2, lossFn = "Binder")
#'  print(s)
#'  plot(s, hm=TRUE) #pdf(height=8.5, width=10.5) #png(height=700, width=720)
#'  plot_ConvDPM(MCMCsample_st, from=2)
#'  #cluster_est_binder(MCMCsample_st$mcmc_partitions[900:1000])
#'
#'  postalpha <- data.frame("alpha"=s$alpha,
#'                          "distribution" = factor(rep("posterior",length(s$alpha)),
#'                          levels=c("prior", "posterior")))
#'
#'  postK <- data.frame("K"=sapply(lapply(postalpha$alpha, "["),
#'                                 function(x){sum(x/(x+0:(n-1)))}))
#'
#'  p <- (ggplot(postalpha, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..), binwidth=.1,
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="blue")
#'        + ggtitle("Posterior distribution of alpha\n")
#'        # Ignore NA values for mean
#'        # Overlay with transparent density plot
#'        + geom_vline(aes(xintercept=mean(alpha, na.rm=TRUE)),
#'                     color="red", linetype="dashed", size=1)
#'      )
#'  p
#'
#'  p <- (ggplot(postK, aes(x=K))
#'        + geom_histogram(aes(y=..density..),
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="blue")
#'        + ggtitle("Posterior distribution of predicted K\n")
#'        # Ignore NA values for mean
#'        # Overlay with transparent density plot
#'        + geom_vline(aes(xintercept=mean(K, na.rm=TRUE)),
#'                     color="red", linetype="dashed", size=1)
#'        #+ scale_x_continuous(breaks=c(0:6)*2, minor_breaks=c(0:6)*2+1)
#'        + scale_x_continuous(breaks=c(1:12))
#'      )
#'  p
#'
#'  p <- (ggplot(drop=FALSE, alpha=.6)
#'        + geom_density(aes(x=alpha, fill=distribution),
#'                       color=NA, alpha=.6,
#'                       data=postalpha)
#'        + geom_density(aes(x=alpha, fill=distribution),
#'                       color=NA, alpha=.6,
#'                       data=prioralpha)
#'        + ggtitle("Prior and posterior distributions of alpha\n")
#'        + scale_fill_discrete(drop=FALSE)
#'        + theme_bw()
#'      )
#'  p
#'
#'  postK <- data.frame("K"=sapply(lapply(s$mcmc_partitions, unique), length))
#'  (ggplot(postK, aes(x=K))
#'        + geom_histogram(aes(y=..density..),
#'                         colour="black", fill="grey45", binwidth=1)
#'        + ggtitle("Posterior distribution of K\n")
#'        # Ignore NA values for mean
#'        # Overlay with transparent density plot
#'        + geom_vline(aes(xintercept=mean(K, na.rm=TRUE)),
#'                     color="red", linetype="dashed", size=1)
#'        + theme_bw()
#'        + scale_x_continuous(breaks=c(3:11))
#'        + xlim(3,11)
#'      )
#'}
#'
#'
#'
#'
#'  # k-means
#'
#'  plot(x=z[1,], y=z[2,], col=kmeans(t(z), centers=4)$cluster,
#'       xlab = "d = 1", ylab= "d = 2", main="k-means with K=4 clusters")
#'
#'  KM <- kmeans(t(z), centers=4)
#'  KMclust <- factor(KM$cluster)
#'  levels(KMclust) <- c("2", "4", "1", "3")
#'  dataKM <- data.frame("X"=z[1,], "Y"=z[2,],
#'                     "Cluster"=as.character(KMclust))
#'  dataCenters <- data.frame("X"=KM$centers[,1],
#'                            "Y"=KM$centers[,2],
#'                            "Cluster"=c("2", "4", "1", "3"))
#'
#'  p <- (ggplot(dataKM)
#'        + geom_point(aes(x=X, y=Y, col=Cluster))
#'        + geom_point(aes(x=X, y=Y, fill=Cluster, order=Cluster),
#'                     data=dataCenters, shape=22, size=5)
#'        + scale_colour_discrete(name="Cluster", guide=guide_legend(override.aes=list(size=6,
#'                                                                                     shape=22)))
#'        + ggtitle("K-means with K=4 clusters\n")
#'        + theme_bw()
#'  )
#'  p
#'
#'
#'
#'
DPMGibbsSkewT <- function (z, hyperG0, a=0.0001, b=0.0001, N, doPlot=TRUE,
                           nbclust_init=30, plotevery=N/10,
                           diagVar=TRUE, use_variance_hyperprior=TRUE, verbose=TRUE,
                           ...){
  
  if(nbclust_init > ncol(z)){
    stop("'nbclust_init' is larger than the number of observations")
  }
  
  if(doPlot){requireNamespace("ggplot2", quietly = TRUE)}
  
  p <- dim(z)[1]
  n <- dim(z)[2]
  U_xi <- matrix(0, nrow=p, ncol=n)
  U_psi <- matrix(0, nrow=p, ncol=n)
  U_Sigma <- array(0, dim=c(p, p, n))
  U_df <- rep(10,n)
  U_B <- array(0, dim=c(2, 2, n))
  U_nu <- rep(p,n)
  
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
  sc <- rep(1,n)
  
  # Initialisation----
  # each observation is assigned to a different cluster
  # or to 1 of the 50 initial clusters if there are more than
  # 50 observations
  
  i <- 1
  
  if(ncol(z)<nbclust_init){
    for (k in 1:n){
      c[k] <- k
      #cat("cluster ", k, ":\n")
      U_SS[[k]] <- update_SSst(z=z[, k, drop=FALSE], S=hyperG0, ltn=ltn[k], scale=sc[k], df=U_df[k])
      NNiW <- rNNiW(U_SS[[k]], diagVar)
      U_xi[, k] <- NNiW[["xi"]]
      U_SS[[k]][["xi"]] <- NNiW[["xi"]]
      U_psi[, k] <- NNiW[["psi"]]
      U_SS[[k]][["psi"]] <- NNiW[["psi"]]
      U_Sigma[, , k] <- NNiW[["S"]]
      U_SS[[k]][["S"]] <- NNiW[["S"]]
      U_B[, ,k] <- U_SS[[k]][["B"]]
      m[k] <- m[k]+1
      U_SS[[k]][["weight"]] <- 1/n
    }
  } else{
    c <- sample(x=1:nbclust_init, size=n, replace=TRUE)
    for (k in unique(c)){
      obs_k <- which(c==k)
      U_SS[[k]] <- update_SSst(z=z[, obs_k, drop=FALSE], S=hyperG0, ltn=ltn[obs_k], scale=sc[obs_k], df=U_df[k])
      NNiW <- rNNiW(U_SS[[k]], diagVar)
      U_xi[, k] <- NNiW[["xi"]]
      U_SS[[k]][["xi"]] <- NNiW[["xi"]]
      U_psi[, k] <- NNiW[["psi"]]
      U_SS[[k]][["psi"]] <- NNiW[["psi"]]
      U_Sigma[, , k] <- NNiW[["S"]]
      U_SS[[k]][["S"]] <- NNiW[["S"]]
      U_B[, ,k] <- U_SS[[k]][["B"]]
      m[k] <- length(obs_k)
      U_SS[[k]][["weight"]] <- m[k]/n
    }
  }
  
  
  alpha <- c(log(n))
  
  
  U_SS_list[[i]] <- U_SS
  c_list[[i]] <- c
  weights_list[[1]] <- numeric(length(m))
  weights_list[[1]][unique(c)] <- table(c)/length(c)
  
  logposterior_list[[i]] <- logposterior_DPMST(z, xi=U_xi, psi=U_psi, Sigma=U_Sigma, df=U_df, B=U_B,
                                               hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b, diagVar)
  
  if(doPlot){
    plot_DPMst(z=z, c=c, i=i, alpha=alpha[i], U_SS=U_SS_list[[i]], ellipses=TRUE, ...)
  }
  if(verbose){
    cat(i, "/", N, " samplings:\n", sep="")
    cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
    cl2print <- unique(c)
    cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
  }
  
  acc_rate <- 0
  
  
  if(N>1){
    for(i in 2:N){
      nbClust <- length(unique(c))
      
      alpha <- c(alpha,
                 sample_alpha(alpha_old=alpha[i-1], n=n,
                              K=nbClust, a=a, b=b)
      )
      
      slice <- sliceSampler_skewT(c=c, m=m, alpha=alpha[i],
                                  z=z, hyperG0=hyperG0,
                                  U_xi=U_xi, U_psi=U_psi,
                                  U_Sigma=U_Sigma, U_df=U_df,
                                  scale=sc, diagVar)
      m <- slice[["m"]]
      c <- slice[["c"]]
      weights_list[[i]] <- slice[["weights"]]
      ltn <- slice[["latentTrunc"]]
      U_xi <- slice[["xi"]]
      U_psi <- slice[["psi"]]
      U_Sigma <- slice[["Sigma"]]
      U_df <- slice[["df"]]
      
      
      # Update cluster locations
      fullCl <- which(m!=0)
      fullCl_nb <- length(fullCl)
      for(k in 1:fullCl_nb){
        j <- fullCl[k]
        obs_j <- which(c==j)
        #cat("cluster ", j, ":\n")
        if(use_variance_hyperprior){
          U_SS[[j]] <- update_SSst(z=z[, obs_j, drop=FALSE], S=hyperG0,
                                   ltn=ltn[obs_j], scale=sc[obs_j],
                                   df=U_df[j],
                                   hyperprior = list("Sigma"=U_Sigma[,,j])
          )
        }else{
          U_SS[[j]] <- update_SSst(z=z[, obs_j, drop=FALSE], S=hyperG0,
                                   ltn=ltn[obs_j], scale=sc[obs_j],
                                   df=U_df[j]
          )
        }
        U_nu[j] <- U_SS[[j]][["nu"]]
        NNiW <- rNNiW(U_SS[[j]], diagVar)
        U_xi[, j] <- NNiW[["xi"]]
        U_SS[[j]][["xi"]] <- NNiW[["xi"]]
        U_psi[, j] <- NNiW[["psi"]]
        U_SS[[j]][["psi"]] <- NNiW[["psi"]]
        U_Sigma[, , j] <- NNiW[["S"]]
        U_SS[[j]][["S"]] <- NNiW[["S"]]
        U_B[, ,j] <- U_SS[[j]][["B"]]
        U_SS[[j]][["weight"]] <- weights_list[[i]][j]
      }
      
      update_scale <- sample_scale(c=c, m=m, z=z, U_xi=U_xi,
                                   U_psi=U_psi, U_Sigma=U_Sigma,
                                   U_df=U_df, ltn=ltn,
                                   weights=weights_list[[i]],
                                   scale=sc)
      U_df_list <- update_scale[["df"]]
      sc <- update_scale[["scale"]]
      acc_rate <- acc_rate + update_scale[["acc_rate"]]
      
      for(k in 1:fullCl_nb){
        j <- fullCl[k]
        U_df[j] <- U_df_list[[k]]
        U_SS[[j]][["df"]] <- U_df[j]
      }
      
      U_SS_list[[i]] <- c(U_SS[which(m!=0)])
      c_list[[i]] <- c
      
      logposterior_list[[i]] <- logposterior_DPMST(z, xi=U_xi, psi=U_psi, Sigma=U_Sigma, df=U_df, B=U_B,
                                                   hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b, diagVar)
      
      if(doPlot && i/plotevery==floor(i/plotevery)){
        plot_DPMst(z=z, c=c, i=i, alpha=alpha[i], U_SS=U_SS_list[[i]], ellipses=TRUE, ...)
      }
      if(verbose){
        cat(i, "/", N, " samplings:\n", sep="")
        cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
        cl2print <- unique(c)
        cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
      }
    }
  }
  acc_rate <- acc_rate/N
  
  dpmclus <- list("mcmc_partitions" = c_list,
                  "alpha"=alpha,
                  "U_SS_list"=U_SS_list,
                  "weights_list"=weights_list,
                  "logposterior_list"=logposterior_list,
                  "data"=z,
                  "nb_mcmcit"=N,
                  "clust_distrib"="skewt",
                  "acc_rate"=acc_rate,
                  "hyperG0"=hyperG0)
  class(dpmclus) <- "DPMMclust"
  return(dpmclus)
}






