#'Slice Sampling of Dirichlet Process Mixture of skew normal ditributions
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
#'@param diagVar logical flag indicating wether the variance of a cluster is a diagonal matrice.
#'Default is \code{FALSE} (full matrix).
#'
#'@param use_variance_hyperprior logical flag indicating whether a hyperprior is added 
#'for the variance parameter. Default is \code{TRUE} which decrease the impact of the variance prior
#'on the posterior. \code{FALSE} is useful for using an informative prior.
#'
#'@param verbose logical flag indicating whether partition info is
#'written in the console at each MCMC iteration.
#'
#'@param ... additional arguments to be passed to \code{\link{plot_DPMsn}}.
#'Only used if \code{doPlot} is \code{TRUE}.
#'
#'@return a object of class \code{DPMclust} with the following attributes:
#'  \itemize{
#'      \item{\code{mcmc_partitions}:}{a list of length \code{N}. Each
#'       element \code{mcmc_partitions[n]} is a vector of length
#'       \code{n} giving the partition of the \code{n} observations.}
#'      \item{\code{alpha}:}{ a vector of length \code{N}. \code{cost[j]} is the cost
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
#'
#' #Number of data
#' n <- 1000
#' set.seed(123)
#'
#' d <- 2
#' ncl <- 4
#'
#' # Sample data
#'
#' sdev <- array(dim=c(d,d,ncl))
#'
#' #xi <- matrix(nrow=d, ncol=ncl, c(-1.5, 1, 1.5, 1, 1.5, -2, -2, -2))
#' xi <- matrix(nrow=d, ncol=ncl, c(-0.5, 0, 0.5, 0, 0.5, -1, -1, 1))
#' ##xi <- matrix(nrow=d, ncol=ncl, c(-0.3, 0, 0.5, 0.5, 0.5, -1.2, -1, 1))
#' psi <- matrix(nrow=d, ncol=4, c(0.4, -0.6, 0.8, 0, 0.3, -0.7, -0.3, -1.2))
#' p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters
#' sdev[, ,1] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.3))
#' sdev[, ,2] <- matrix(nrow=d, ncol=d, c(0.1, 0, 0, 0.3))
#' sdev[, ,3] <- matrix(nrow=d, ncol=d, c(0.3, 0.15, 0.15, 0.3))
#' sdev[, ,4] <- .3*diag(2)
#'
#'
#' c <- rep(0,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- xi[, c[k]] + psi[, c[k]]*abs(rnorm(1)) + sdev[, , c[k]]%*%matrix(rnorm(d, mean = 0,
#'                                                                        sd = 1), nrow=d, ncol=1)
#'  #cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'
#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["b_xi"]] <- rep(0,d)
#' hyperG0[["b_psi"]] <- rep(0,d)
#' hyperG0[["kappa"]] <- 0.0001
#' hyperG0[["D_xi"]] <- 100
#' hyperG0[["D_psi"]] <- 100
#' hyperG0[["nu"]] <- d + 1
#' hyperG0[["lambda"]] <- diag(d)
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
#'        + ggtitle("Simple example in 2d data")
#'        +xlab("D1")
#'        +ylab("D2")
#'        +theme_bw())
#'  p
#'
#'  c2plot <- factor(c)
#'  levels(c2plot) <- c("3", "2", "4", "1")
#'  pp <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,], "Cluster"=as.character(c2plot)))
#'        + geom_point(aes(x=X, y=Y, colour=Cluster, fill=Cluster))
#'        + ggtitle("Slightly overlapping skew-normal simulation\n")
#'        + xlab("D1")
#'        + ylab("D2")
#'        + theme_bw()
#'        + scale_colour_discrete(guide=guide_legend(override.aes = list(size = 6, shape=22))))
#'  pp
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
#'\dontrun{
#'  # Gibbs sampler for Dirichlet Process Mixtures
#'  ##############################################
#'
#'  MCMCsample_sn <- DPMGibbsSkewN(z, hyperG0, a, b, N=2500,
#'                                 doPlot, nbclust_init, plotevery=100,
#'                                 gg.add=list(theme_bw(),
#'                                  guides(shape=guide_legend(override.aes = list(fill="grey45")))),
#'                                diagVar=FALSE)
#'
#'  s <- summary(MCMCsample_sn, burnin = 2000, thin=10)
#'  #cluster_est_binder(MCMCsample_sn$mcmc_partitions[1000:1500])
#'  print(s)
#'  plot(s)
#'  #plot_ConvDPM(MCMCsample_sn, from=2)
#'
#'
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
#'
#'  p <- (ggplot(dataKM)
#'        + geom_point(aes(x=X, y=Y, col=Cluster))
#'        + geom_point(aes(x=X, y=Y, fill=Cluster, order=Cluster),
#'                     data=dataCenters, shape=22, size=5)
#'        + scale_colour_discrete(name="Cluster",
#'                                guide=guide_legend(override.aes=list(size=6, shape=22)))
#'        + ggtitle("K-means with K=4 clusters\n")
#'        + theme_bw()
#'  )
#'  p
#'
#'  postalpha <- data.frame("alpha"=MCMCsample_sn$alpha[501:1000],
#'                          "distribution" = factor(rep("posterior",1000-500),
#'                          levels=c("prior", "posterior")))
#'
#'  postK <- data.frame("K"=sapply(lapply(postalpha$alpha, "["),
#'                                 function(x){sum(x/(x+0:(1000-1)))}))
#'
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
#'  p <- (ggplot(postK, aes(x=K))
#'        + geom_histogram(aes(y=..density..),
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="blue")
#'        + ggtitle("Posterior distribution of predicted K\n")
#'        # Ignore NA values for mean
#'        # Overlay with transparent density plot
#'        + geom_vline(aes(xintercept=mean(K, na.rm=T)),
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
#'        + xlim(0,100)
#'      )
#'  p
#'
#'#Skew Normal
#'n=100000
#' xi <- 0
#' d <- 0.995
#' alpha <- d/sqrt(1-d^2)
#' z <- rtruncnorm(n,a=0, b=Inf)
#' e <- rnorm(n, mean = 0, sd = 1)
#' x <- d*z + sqrt(1-d^2)*e
#'o <- 1
#' y <- xi+o*x
#'nu=1.3
#'w <- rgamma(n,scale=nu/2, shape=nu/2)
#'yy <- xi+o*x/w
#'snd <- data.frame("Y"=y,"YY"=yy)
#'p <- (ggplot(snd)+geom_density(aes(x=Y), fill="blue", alpha=.2)
#'      + theme_bw()
#'      + ylab("Density")
#'      + ggtitle("Y~SN(0,1,10)\n")
#'      + xlim(-1,6)
#'      + ylim(0,0.8)
#'      )
#'p
#'
#'p <- (ggplot(snd)+geom_density(aes(x=YY), fill="blue", alpha=.2)
#'      + theme_bw()
#'      + ylab("Density")
#'      + ggtitle("Y~ST(0,1,10,1.3)\n")
#'      + xlim(-2,40)
#'      + ylim(0,0.8)
#'      )
#'p
#'
#'p <- (ggplot(snd)
#'      + geom_density(aes(x=Y, fill="blue"), alpha=.2)
#'      + geom_density(aes(x=YY, fill="red"), alpha=.2)
#'      + theme_bw()
#'      +theme(legend.text = element_text(size = 13), legend.position="bottom")
#'      + ylab("Density")
#'      + xlim(-2,40)
#'      + ylim(0,0.8)
#'      + scale_fill_manual(name="", labels=c("Y~SN(0,1,10)       ", "Y~ST(0,1,10,1.3)"),
#'      guide="legend", values=c("blue", "red"))
#'      )
#'p
#'
#'#Variations
#'n=100000
#' xi <- -1
#' d <- 0.995
#' alpha <- d/sqrt(1-d^2)
#' z <- rtruncnorm(n,a=0, b=Inf)
#' e <- rnorm(n, mean = 0, sd = 1)
#' x <- d*z + sqrt(1-d^2)*e
#'snd <- data.frame("X"=x)
#'p <- (ggplot(snd)+geom_density(aes(x=X), fill="blue", alpha=.2)
#'      + theme_bw()
#'      + ylab("Density")
#'      + ggtitle("X~SN(10)\n")
#'      + xlim(-1.5,4)
#'      + ylim(0,1.6)
#'      )
#'p
#'
#'o <- 0.5
#' y <- xi+o*x
#'snd <- data.frame("Y"=y)
#'p <- (ggplot(snd)+geom_density(aes(x=Y), fill="blue", alpha=.2)
#'      + theme_bw()
#'      + ylab("Density")
#'      + ggtitle("Y~SN(-1,1,10)\n")
#'      + xlim(-1.5,4)
#'      + ylim(0,1.6)
#'      )
#'p
#'
#'
#'
#'
#'#Simple toy example
#'###################
#'
#' n <- 500
#' set.seed(12345)
#'
#'
#' d <- 2
#' ncl <- 4
#'
#' # Sample data
#'
#' sdev <- array(dim=c(d,d,ncl))
#'
#' xi <- matrix(nrow=d, ncol=ncl, c(-1.5, 1, 1.5, 1, 1.5, -2, -2, -2))
#' psi <- matrix(nrow=d, ncol=4, c(0.4, -0.6, 0.8, 0, 0.3, -0.7, -0.3, -1.2))
#' p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters
#' sdev[, ,1] <- matrix(nrow=d, ncol=d, c(0.3, 0, 0, 0.3))
#' sdev[, ,2] <- matrix(nrow=d, ncol=d, c(0.1, 0, 0, 0.3))
#' sdev[, ,3] <- matrix(nrow=d, ncol=d, c(0.3, 0.15, 0.15, 0.3))
#' sdev[, ,4] <- .3*diag(2)
#'
#'#' # Set parameters of G0
#' hyperG0 <- list()
#' hyperG0[["b_xi"]] <- rep(0,d)
#' hyperG0[["b_psi"]] <- rep(0,d)
#' hyperG0[["kappa"]] <- 0.0001
#' hyperG0[["D_xi"]] <- 100
#' hyperG0[["D_psi"]] <- 100
#' hyperG0[["nu"]] <- d + 1
#' hyperG0[["lambda"]] <- diag(d)
#'
#'
#' c <- rep(0,n)
#' z <- matrix(0, nrow=d, ncol=n)
#' for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- xi[, c[k]] + psi[, c[k]]*abs(rnorm(1)) + sdev[, , c[k]]%*%matrix(rnorm(d, mean = 0,
#'                                                                         sd = 1), nrow=d, ncol=1)
#'  cat(k, "/", n, " observations simulated\n", sep="")
#' }
#'
#'  MCMCsample_sn_sep <- DPMGibbsSkewN(z, hyperG0, a, b, N=600,
#'                                     doPlot, nbclust_init, plotevery=100,
#'                                     gg.add=list(theme_bw(),
#'                                guides(shape=guide_legend(override.aes = list(fill="grey45")))),
#'                                     diagVar=TRUE)
#'
#'  s <- summary(MCMCsample_sn, burnin = 400)
#'
#'}
#'
DPMGibbsSkewN <- function (z, hyperG0, a=0.0001, b=0.0001, N, doPlot=TRUE,
                           nbclust_init=30, plotevery=N/10,
                           diagVar=TRUE, use_variance_hyperprior=TRUE, verbose=TRUE,
                           ...){
  
  if(doPlot){requireNamespace("ggplot2", quietly=TRUE)}
  
  p <- dim(z)[1]
  n <- dim(z)[2]
  U_xi <- matrix(0, nrow=p, ncol=n)
  U_psi <- matrix(0, nrow=p, ncol=n)
  U_Sigma = array(0, dim=c(p, p, n))
  U_B = array(0, dim=c(2, 2, n))
  
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
      #cat("cluster ", k, ":\n")
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
      #cat("cluster ", k, ":\n")
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
    cat(i, "/", N, " samplings:\n", sep="")
    cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
    cl2print <- unique(c)
    cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
  }
  
  
  if(N>1){
    for(i in 2:N){
      
      nbClust <- length(unique(c))
      
      alpha <- c(alpha,
                 sample_alpha(alpha_old=alpha[i-1], n=n,
                              K=nbClust, a=a, b=b)
      )
      
      slice <- sliceSampler_SkewN(c=c, m=m, alpha=alpha[i],
                                  z=z, hyperG0=hyperG0,
                                  U_xi=U_xi, U_psi=U_psi,
                                  U_Sigma=U_Sigma, diagVar)
      m <- slice[["m"]]
      c <- slice[["c"]]
      weights_list[[i]] <- slice[["weights"]]
      ltn <- slice[["latentTrunc"]]
      
      
      
      # Update cluster locations
      fullCl <- which(m!=0)
      for(j in fullCl){
        obs_j <- which(c==j)
        #cat("cluster ", j, ":\n")
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
        cat(i, "/", N, " samplings:\n", sep="")
        cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
        cl2print <- unique(c)
        cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
      }
    }
  }
  
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
  return(dpmclus)
}






