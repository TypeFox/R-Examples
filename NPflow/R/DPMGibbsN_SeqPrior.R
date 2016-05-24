
#'Slice Sampling of Dirichlet Process Mixture of Gaussian distibutions
#'
#'@param z data matrix \code{d x n} with \code{d} dimensions in rows
#'and \code{n} observations in columns.
#'
#'@param prior_inform an informative prior such as the approximation computed by \code{summary.DPMMclust}.
#'
#'@param hyperG0 a non informative prior component for the mixing distribution.
#'Only used if \code{add.vagueprior} is \code{TRUE}.
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
#'@param add.vagueprior logical flag indicating wether a non informative component should
#' be added to the informative prior. Default is \code{TRUE}.
#'
#'@param weightnoninfo a real between 0 and 1 giving the weights of the non informative component
#'in the prior.
#'
#'@param diagVar logical flag indicating wether the variance of each cluster is
#'estimated as a diagonal matrix, or as a full matrix.
#'Default is \code{TRUE} (diagonal variance).
#'
#'@param verbose logical flag indicating wether partition info is
#'written in the console at each MCMC iteration.
#'
#'@param ... additional arguments to be passed to \code{\link{plot_DPM}}.
#'Only used if \code{doPlot} is \code{TRUE}.
#'
#'@return a object of class \code{DPMclust} with the following attributes:
#'  \itemize{
#'      \item{\code{mcmc_partitions}:}{ a list of length \code{N}. Each
#'       element \code{mcmc_partitions[n]} is a vector of length
#'       \code{n} giving the partition of the \code{n} observations.}
#'      \item{\code{alpha}:}{ a vector of length \code{N}. \code{cost[j]} is the cost
#' associated to partition \code{c[[j]]}}
#'       \item{\code{listU_mu}:}{ a list of length \code{N} containing the matrices of
#'       mean vectors for all the mixture components at each MCMC iteration}
#'       \item{\code{listU_Sigma}:}{ a list of length \code{N} containing the arrays of
#'       covariances matrices for all the mixture components at each MCMC iteration}
#'       \item{\code{U_SS_list}:}{ a list of length \code{N} containing the lists of
#'       sufficient statistics for all the mixture components at each MCMC iteration}
#'      \item{\code{weights_list}:}{}
#'      \item{\code{logposterior_list}:}{ a list of length \code{N} containing the logposterior values
#'       at each MCMC iterations}
#'      \item{\code{data}:}{ the data matrix \code{d x n} with \code{d} dimensions in rows
#'and \code{n} observations in columns.}
#'      \item{\code{nb_mcmcit}:}{ the number of MCMC itertations}
#'      \item{\code{clust_distrib}:}{ the parametric distribution of the mixture component - \code{"gaussian"}}
#'      \item{\code{hyperG0}:}{ the prior on the cluster location}
#'  }
#'
#'@author Boris Hejblum, Chariff Alkhassim
#'
#'@seealso \code{\link{postProcess.DPMMclust}} \code{\link{DPMGibbsN}}
#'
#'@references Hejblum BP, Alkhassim C, Gottardo R, Caron F, Thiebaut R, Sequential Dirichlet
#'Process Mixtures of Multivariate Skew t-distributions for Model-based Clustering
#'of Flow Cytometry Data, in preparation.
#'
#'@export
#'
#'@examples
#'
#'rm(list=ls())
#'#Number of data
#'n <- 1500
#'# Sample data
#'#m <- matrix(nrow=2, ncol=4, c(-1, 1, 1.5, 2, 2, -2, 0.5, -2))
#'m <- matrix(nrow=2, ncol=4, c(-.8, .7, .5, .7, .5, -.7, -.5, -.7))
#'p <- c(0.2, 0.1, 0.4, 0.3) # frequence des clusters
#'
#'sdev <- array(dim=c(2,2,4))
#'sdev[, ,1] <- matrix(nrow=2, ncol=2, c(0.3, 0, 0, 0.3))
#'sdev[, ,2] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.3))
#'sdev[, ,3] <- matrix(nrow=2, ncol=2, c(0.3, 0.15, 0.15, 0.3))
#'sdev[, ,4] <- .3*diag(2)
#'c <- rep(0,n)
#'z <- matrix(0, nrow=2, ncol=n)
#'for(k in 1:n){
#'  c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
#'  z[,k] <- m[, c[k]] + sdev[, , c[k]]%*%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
#'  #cat(k, "/", n, " observations simulated\n", sep="")
#'}
#'
#'d<-2
#'# Set parameters of G0
#'hyperG0 <- list()
#'hyperG0[["mu"]] <- rep(0,d)
#'hyperG0[["kappa"]] <- 0.001
#'hyperG0[["nu"]] <- d+2
#'hyperG0[["lambda"]] <- diag(d)/10
#'
#'# hyperprior on the Scale parameter of DPM
#'a <- 0.0001
#'b <- 0.0001
#'
#'# Number of iterations
#'N <- 30
#'
#'# do some plots
#'doPlot <- TRUE
#'nbclust_init <- 20
#'
#'
#'
#'## Data
#'########
#'library(ggplot2)
#'p <- (ggplot(data.frame("X"=z[1,], "Y"=z[2,]), aes(x=X, y=Y))
#'       + geom_point()
#'       + ggtitle("Toy example Data"))
#'p
#'
#'
#'
#' # Gibbs sampler for Dirichlet Process Mixtures
#' ##############################################
#' \dontrun{
#' MCMCsample <- DPMGibbsN(z, hyperG0, a, b, N=1500, doPlot, nbclust_init, plotevery=200,
#'                         gg.add=list(theme_bw(),
#'                                  guides(shape=guide_legend(override.aes = list(fill="grey45")))),
#'                         diagVar=FALSE)
#'
#' s <- summary(MCMCsample, posterior_approx=TRUE, burnin = 1000, thin=5)
#' F1 <- FmeasureC(pred=s$point_estim$c_est, ref=c)
#' F1
#'
#'
#' MCMCsample2 <- DPMGibbsN_SeqPrior(z, prior_inform=s$param_posterior,
#'                                   hyperG0, N=1500,
#'                                   add.vagueprior = TRUE,
#'                                   doPlot=TRUE, plotevery=100,
#'                                   nbclust_init=nbclust_init,
#'                                   gg.add=list(theme_bw(),
#'                                  guides(shape=guide_legend(override.aes = list(fill="grey45")))),
#'                                   diagVar=FALSE)
#'
#'
#' s2 <- summary(MCMCsample2, burnin = 500, thin=5)
#' F2 <- FmeasureC(pred=s2$point_estim$c_est, ref=c)
#' F2
#' }


DPMGibbsN_SeqPrior <- function (z, prior_inform, hyperG0, N, nbclust_init,
                                add.vagueprior = TRUE, weightnoninfo=NULL,
                                doPlot=TRUE, plotevery=N/10,
                                diagVar=TRUE, verbose=TRUE,
                                ...){

  if(nbclust_init > ncol(z)){
    stop("'nbclust_init' is larger than the number of observations")
  }

  if(doPlot){requireNamespace("ggplot2", quietly=TRUE)}

  p <- nrow(z)
  n <- ncol(z)
  U_mu <- matrix(0, nrow=p, ncol=n)
  U_Sigma = array(0, dim=c(p, p, n))
  listU_mu<-list()
  listU_Sigma<-list()

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


  priorG1 <- prior_inform
  nonnullpriors_ind <- which(priorG1$weights!=0)
  priorG1$weights <- priorG1$weights[nonnullpriors_ind]
  priorG1$parameters <- priorG1$parameters[nonnullpriors_ind]
  nbmix_prior <- length(priorG1[["weights"]])
  if(add.vagueprior){
    nbmix_prior <- nbmix_prior + 1
    priorG1[["parameters"]][[nbmix_prior]] <- hyperG0

    if(is.null(weightnoninfo)){
      priorG1$weights <- c(priorG1$weights, 1/length(priorG1$weights))
      priorG1$weights <- priorG1$weights/sum(priorG1$weights)
    }else{
      #TODO
      priorG1$weights <- c(rep((1-weightnoninfo)/(nbmix_prior-1), (nbmix_prior-1)), weightnoninfo)
    }
  }


  a <- prior_inform$alpha_param$shape
  b <- prior_inform$alpha_param$rate

  #         a <- 0.00001
  #         b <- 0.00001


  # Initialisation----
  # each observation is assigned to cluster
  i <- 1
  c <- sample(1:nbclust_init, size=n, replace=TRUE)

  for (k in unique(c)){
    obs_k <- which(c==k)
    hyper_num <- sample(x=1:nbmix_prior, size=1)#, prob=priorG1$weights)
    priormix <- priorG1[["parameters"]][[hyper_num]]

    U_SS[[k]] <- update_SS(z=z[, obs_k], S=priormix)
    NiW <- rNiW(U_SS[[k]], diagVar)
    U_mu[, k] <- NiW[["mu"]]
    U_SS[[k]][["mu"]] <- NiW[["mu"]]
    U_Sigma[, , k] <- NiW[["S"]]
    U_SS[[k]][["S"]] <- NiW[["S"]]
    m[k] <- length(obs_k)
    U_SS[[k]][["weight"]] <-m[k]/n
  }

  listU_mu[[i]]<-U_mu
  listU_Sigma[[i]]<-U_Sigma

  if(is.null(hyperG0[["alpha"]])){
    alpha <- nbmix_prior/log(n)
  }else{
    alpha <- hyperG0[["alpha"]]
  }

  alpha <- nbmix_prior/log(n)
  U_SS_list[[i]] <- U_SS[which(m!=0)]
  c_list[[i]] <- c
  weights_list[[1]] <- numeric(length(m))
  weights_list[[1]][sort(unique(c))] <- table(c)/length(c)

  logposterior_list[[i]] <- logposterior_DPMG(z, mu=U_mu, Sigma=U_Sigma,
                                             hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b)


  if(doPlot){
    plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma,
             m=m, c=c, i=i, alpha=alpha[i], U_SS=U_SS_list[[i]], ...)
  }
  if(verbose){
    cat(i, "/", N, " samplings:\n", sep="")
    cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
    cl2print <- unique(c)
    cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
  }

  # MCMC
  for(i in 2:N){
    nbClust <- length(unique(c))
    alpha <- c(alpha,sample_alpha(alpha_old=alpha[i-1], n=n,
                                  K=nbClust, a=a, b=b))


    slice <- sliceSampler_N_SeqPrior(c=c, m=m, alpha=alpha[i],
                                     z=z, priorG1=priorG1,
                                     U_mu=U_mu, U_Sigma=U_Sigma, diagVar=diagVar)

    m <- slice[["m"]]
    c <- slice[["c"]]
    U_mu<-slice[["U_mu"]]
    U_Sigma<-slice[["U_Sigma"]]
    weights_list[[i]] <- slice[["weights"]]

    # Update cluster locations
    fullCl <- which(m!=0)
    fullCl_nb <- length(fullCl)

    U_SS_prior <- list()
    p <- matrix(nrow=nbmix_prior, ncol=fullCl_nb)
    vrais <- rep(NA,fullCl_nb)



    for(k in 1:fullCl_nb){
      j <- fullCl[k]
      obs_j <- which(c==j)
      U_SS_prior[[k]] <- list()
      for(l in 1:nbmix_prior){
        U_SS_prior[[k]][[l]] <- update_SS(z=z[, obs_j, drop=FALSE],
                                          S=priorG1[["parameters"]][[l]])

      }


      p[,k] <- mmNiWpdfC(Mu=U_mu[,j, drop=FALSE], Sigma=list(U_Sigma[,,j]),
                         U_Mu0=sapply(U_SS_prior[[k]], "[[", "mu"),
                         U_Kappa0=sapply(U_SS_prior[[k]], "[[", "kappa"),
                         U_Nu0=sapply(U_SS_prior[[k]], "[[", "nu"),
                         U_Sigma0=lapply(U_SS_prior[[k]], "[[", "lambda"),
                         Log=TRUE)

      vrais[k] <- sum(mmvnpdfC(x=z[,obs_j, drop=FALSE], mean=U_mu[,j, drop=FALSE],
                               varcovM=list(U_Sigma[,,j]), Log=TRUE))

    }

    p0 <- mmNiWpdfC(Mu=U_mu[, fullCl, drop=FALSE],
                    Sigma=lapply(fullCl, function(m) U_Sigma[, ,m]),
                    U_Mu0=sapply(priorG1[["parameters"]], "[[", "mu"),
                    U_Kappa0=sapply(priorG1[["parameters"]], "[[", "kappa"),
                    U_Nu0=sapply(priorG1[["parameters"]], "[[", "nu"),
                    U_Sigma0=lapply(priorG1[["parameters"]], "[[", "lambda"),
                    Log=TRUE)

    pfin_log <- apply(X=(p0 - p), MARGIN=1, FUN=function(r){vrais + r})


    if(is.null(dim(pfin_log))){

      #only one sampled cluster non empty
      logexptrick_const <- max(pfin_log)
      wfin_log_const <- pfin_log - logexptrick_const
      w2fin <- exp(wfin_log_const)*priorG1[["weights"]]
      w2fin_sums <- sum(w2fin)
      wfin <- matrix(w2fin/w2fin_sums, nrow=1)
    }else if(ncol(pfin_log)==1){
      #only one prior component
      logexptrick_const <- apply(X=pfin_log, MARGIN=1, FUN=max)
      wfin_log_const <- apply(X=pfin_log, MARGIN=2, FUN=function(cv){cv - logexptrick_const})
      w2fin <- apply(X=exp(wfin_log_const), MARGIN=1, FUN=function(r){r*priorG1[["weights"]]})
      wfin <- matrix(w2fin, ncol=1)
    }else{
      logexptrick_const <- apply(X=pfin_log, MARGIN=1, FUN=max)
      wfin_log_const <- apply(X=pfin_log, MARGIN=2, FUN=function(cv){cv - logexptrick_const})
      w2fin <- apply(X=exp(wfin_log_const), MARGIN=1, FUN=function(r){r*priorG1[["weights"]]})
      w2fin_sums <- colSums(w2fin)
      #w2fin_sums0_ind <- which(w2fin_sums==0)
      #             if(length(w2fin_sums0_ind)>0){
      #                 w2fin[nbmix_prior,w2fin_sums0_ind] <- 1
      #                 w2fin_sums[w2fin_sums0_ind] <- 1
      #             }
      wfin <- apply(X=w2fin, MARGIN=1, FUN=function(r){r/w2fin_sums})
      #browser()
      #any(rowSums(wfin)!=1) #should all be 1
    }

    if(any(is.nan(wfin))){
      na_rws <- unique(which(is.nan(wfin), arr.ind=TRUE)[,"row"])
      for(ro in na_rws){
        wfin[ro,] <- priorG1[["weights"]]
      }
    }


    for(k in 1:fullCl_nb){
      j <- fullCl[k]
      obs_j <- which(c==j)
      #cat("cluster ", j, ":\n")

      #sample prior mixture element to update
      hyper_num <- sample(x=1:nbmix_prior, size=1, prob=wfin[k,])
      priormix <- priorG1[["parameters"]][[hyper_num]]
      #cat(hyper_num, " ")

      U_SS[[j]] <- update_SS(z=z[, obs_j, drop=FALSE], S=priormix)


      NiW <- rNiW(U_SS[[j]], diagVar)
      U_mu[, j] <- NiW[["mu"]]
      U_SS[[j]][["mu"]] <- NiW[["mu"]]
      U_Sigma[, , j] <- NiW[["S"]]
      U_SS[[j]][["S"]] <- NiW[["S"]]
      U_SS[[j]][["weight"]] <-weights_list[[i]][j]

    }

    listU_mu[[i]]<-U_mu
    listU_Sigma[[i]]<-U_Sigma
    U_SS_list[[i]] <- U_SS[which(m!=0)]
    c_list[[i]] <- c

    logposterior_list[[i]] <- logposterior_DPMG(z, mu=U_mu, Sigma=U_Sigma,
                                               hyper=hyperG0, c=c, m=m, alpha=alpha[i], n=n, a=a, b=b)


    if(doPlot && i/plotevery==floor(i/plotevery)){
      plot_DPM(z=z, U_mu=U_mu, U_Sigma=U_Sigma, m=m, c=c, i=i,
               alpha=alpha[i], U_SS=U_SS, ...)
    }
    if(verbose){
      cat(i, "/", N, " samplings:\n", sep="")
      cat("  logposterior = ", sum(logposterior_list[[i]]), "\n", sep="")
      cl2print <- unique(c)
      cat(length(cl2print), "clusters:", cl2print[order(cl2print)], "\n\n")
    }
  }
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
                  #"acc_rate"=acc_rate,
                  "hyperG0"=hyperG0)
  class(dpmclus) <- "DPMMclust"
  return(dpmclus)
}






