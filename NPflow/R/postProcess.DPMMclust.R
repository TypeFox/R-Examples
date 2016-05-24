#'Post-processing Dirichlet Process Mixture Models results to get
#'a mixture distribution of the posterior locations
#'
#'@param x a \code{DPMMclust} object.
#'
#'@param burnin integer giving the number of MCMC iterations to burn (defaults is half)
#'
#'@param thin integer giving the spacing at which MCMC iterations are kept.
#'Default is \code{1}, i.e. no thining.
#'
#'@param lossFn character string specifying the loss function to be used.
#'Either "F-measure" or "Binder" (see Details). Default is "F-measure".
#'
#'@param gs optionnal vector of length \code{n} containing the gold standard
#'partition of the \code{n} observations to compare to the point estimate.
#'
#'@param K integer giving the number of mixture components. Default is \code{10}.
#'
#'@param ... further arguments passed to or from other methods
#'
#'@return a \code{list}:
#'  \itemize{
#'      \item{\code{burnin}:}{an integer passing along the \code{burnin} argument}
#'      \item{\code{thin}:}{an integer passing along the \code{thin} argument}
#'      \item{\code{lossFn}:}{a character string passing along the \code{lossFn} argument}
#'      \item{\code{point_estim}:}{}
#'      \item{\code{loss}:}{}
#'      \item{\code{index_estim}:}{}
#'  }
#'
#'@details The cost of a point estimate partition is calculated using either a pairwise
#' coincidence loss function (Binder), or 1-Fmeasure (F-measure).
#'
#'@author Boris Hejblum
#'
#'@importFrom stats uniroot
#'
#'@export
#'
#'@importFrom gplots heatmap.2
#'
#'@seealso \code{\link{similarityMat}} \code{\link{summary.DPMMclust}}
#'
postProcess.DPMMclust <- function(x, burnin=0, thin=1, gs=NULL, lossFn="F-measure", K=10, ...){

  x_invar <- burn.DPMMclust(x, burnin = burnin, thin=thin)

  EM_init_nb_max <- 10
  elem <- which(lapply(x_invar$U_SS_list,FUN=length)==K)
  len <- length(elem)
  EM_init <- list()
  if(len>=EM_init_nb_max){
    EM_init_nb <- EM_init_nb_max
    randind <- elem[sample(1:len,EM_init_nb,replace=FALSE)]
    for (el in 1:EM_init_nb){
      EM_init[[el]] <- x_invar$U_SS_list[[randind[el]]]
    }
  }
  else{
    EM_init_nb <- len
    cpt <- 1
    for (el in elem){
      EM_init[[cpt]] <- x_invar$U_SS_list[[el]]
      cpt <- cpt+1
    }
  }


  if(x$clust_distrib=="skewt"){

    xi_list <- list()
    psi_list <- list()
    S_list <- list()
    w_list <- list()

    #m_final <- list()
    #S_final <- list()

    for(i in 1:length(x_invar$U_SS_list)){
      xi_list <- c(xi_list, sapply(x_invar$U_SS_list[[i]], "[", "xi"))
      psi_list <- c(psi_list, sapply(x_invar$U_SS_list[[i]], "[", "psi"))

      S_list <- c(S_list, sapply(x_invar$U_SS_list[[i]], "[", "S"))

      if(is.null(x_invar$U_SS_list[[1]][["weights"]])){
        #for compatibility with older DPMclust objects
        w_list <- c(w_list, x_invar$weights_list[[i]][unique(x_invar$mcmc_partitions[[i]])])
      }else{
        w_list <- c(w_list,sapply(x_invar$U_SS_list[[i]], "[", "weights"))
      }
    }

    mle_g <- MLE_gamma(x_invar$alpha)

    if(K>1){

      MAPprior <- x_invar$hyperG0
      #MAPprior$lambda <-10*MAPprior$lambda
      param_post_list <- list()

      chr_str <- paste(paste("/", as.character(EM_init_nb),sep=""), "computed",sep=" ")
      for (j in 1:EM_init_nb){
        param_post_list[[j]] <- MAP_sNiW_mmEM(xi_list, psi_list, S_list,
                                               hyperG0 = MAPprior, K=K,
                                               init=EM_init[[j]],verbose=FALSE,...)

        cat("EM ", j,chr_str, "\n", sep="")
      }
      param_post <- param_post_list[[which.max(sapply(lapply(param_post_list, "[[", "loglik"), FUN=max))]]


      parameters <- list()
      for (i in 1:length(param_post$U_xi)){
        parameters[[i]] <- list("b_xi" = param_post[["U_xi"]][[i]],
                                "b_psi" = param_post[["U_psi"]][[i]],
                                "B" = param_post[["U_B"]][[i]],
                                "lambda" = param_post[["U_Sigma"]][[i]],
                                "nu" = param_post[["U_df"]][[i]]
        )
      }
    }
    else{
      param_post <- MLE_sNiW(xi_list, psi_list, S_list, ...)
      parameters <- list()
      parameters[[1]] <- list("b_xi" = param_post[["U_xi"]],
                              "b_psi" = param_post[["U_psi"]],
                              "B" = param_post[["U_B"]],
                              "lambda" = param_post[["U_Sigma"]],
                              "nu" = param_post[["U_df"]]
      )
      param_post$weights <- 1
    }


  }else if (x$clust_distrib=="gaussian"){

    mle_g <- MLE_gamma(x_invar$alpha)

    mu_list <- list()
    S_list <- list()
    w_list <- list()


    for(i in 1:length(x_invar$U_SS_list)){
      mu_list <- c(mu_list, sapply(x_invar$U_SS_list[[i]], "[", "mu"))
      S_list <- c(S_list, sapply(x_invar$U_SS_list[[i]], "[", "S"))

      if(is.null(x_invar$U_SS_list[[1]][["weights"]])){
        #for compatibility with older DPMclust objects
        w_list <- c(w_list, x_invar$weights_list[[i]][unique(x_invar$mcmc_partitions[[i]])])
      }else{
        w_list <- c(w_list,sapply(x_invar$U_SS_list[[i]], "[", "weights"))
      }
    }

    param_post_list <- list()
    for (j in 1:EM_init_nb){

      param_post_list[[j]] <- MLE_NiW_mmEM(mu_list, S_list, x_invar$hyperG0, K, maxit=100, tol=1E-1, doPlot=TRUE)
      cat("EM ", j, "/10 computed", "\n", sep="")
    }
    param_post <- param_post_list[[which.max(sapply(lapply(param_post_list, "[[", "loglik"), FUN=max))]]


    parameters <- list()
    for (i in 1:length(param_post$U_mu)){
      parameters[[i]] <- list("mu" = as.vector(param_post[["U_mu"]][[i]]),
                              "kappa" = param_post[["U_kappa"]][[i]],
                              "lambda" = param_post[["U_lambda"]][[i]],
                              "nu" = param_post[["U_nu"]][[i]]
      )
    }

  }
  else {stop("clust_distrib is neither 'skewt' nor 'gaussian'\n other distributions are not implemented yet")}

  return(list("parameters"=parameters, "weights"=param_post$weights,
              "alpha_param"=mle_g))

}


