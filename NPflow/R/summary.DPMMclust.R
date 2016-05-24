#'Summarizing Dirichlet Process Mixture Models
#'
#'Summary methods for \code{DPMMclust} objects.
#'
#'@param object a \code{DPMMclust} object.
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
#'partition of the \code{n} observations to compare to the point estimate
#'
#'@param posterior_approx logical flag wether a parametric approximation of the posterior should be
#'computed. Default is \code{FALSE}
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
#'@export
#'
#'@importFrom pheatmap pheatmap
#'
#'@importFrom gplots heatmap.2
#'
#'@importFrom fastcluster hclust
#'
#'@seealso \code{\link{similarityMat}} \code{\link{similarityMatC}}
#'

summary.DPMMclust <- function(object, burnin=0, thin=1, gs=NULL, lossFn="F-measure",
                              posterior_approx=FALSE, ...){

  x_invar <- burn.DPMMclust(object, burnin = burnin, thin=thin)



  #if(!posterior_approx){
  if(lossFn == "F-measure"){
    point_estim <- cluster_est_Fmeasure(x_invar$mcmc_partitions,
                                        logposterior = sapply(x_invar$logposterior_list, sum))
  }else if(lossFn == "Binder"){
    point_estim <- cluster_est_binder(x_invar$mcmc_partitions,
                                      logposterior = sapply(x_invar$logposterior_list, sum))
  }else if(lossFn == "MBinderN"){
    if (x_invar$clust_distrib!="gaussian"){stop("Clusters distribution must be Gaussian")}
    point_estim <- cluster_est_Mbinder_norm(x_invar$mcmc_partitions,x_invar$listU_mu,
                                            x_invar$listU_Sigma,
                                            logposterior = sapply(x_invar$logposterior_list, sum), ...)
  }else{
    stop("Specified loss function not available.\n
         Specify either 'F-measure' or 'Binder' for the lossFn argument.")
  }

  index_estim <- point_estim$opt_ind
  loss <- NA
  if(!is.null(gs)){
    loss <- evalClustLoss(c=point_estim$c_est, gs=gs, lossFn=lossFn)
  }
  #}
  #point_estim=NULL
  #index_estim <- NA
  #loss <- NA
  #Posterior approximation
  K<-length(unique(point_estim$c_est))
  if(posterior_approx){
    if(K>1){
      param_post <- postProcess.DPMMclust(x_invar, plot=FALSE, K=K, ...)
    }else{
      param_post <- postProcess.DPMMclust(x_invar, K=1, ...)
    }
  }
  else{
    param_post <- NULL
  }

  s <- c(x_invar, list("burnin"=burnin,
                       "thin"=thin,
                       "point_estim"=point_estim,
                       "loss"=loss,
                       "index_estim"=index_estim,
                       "param_posterior"=param_post))
  class(s) <- "summaryDPMMclust"

  invisible(s)
}


#' Methods for a summary of a 'DPMMclust' object
#'@rdname methods.summaryDPMMclust
#'@aliases summaryDPMMclust print.summaryDPMMclust plot.summaryDPMMclust
#'@param x a \code{summaryDPMMclust} object.
#'@param hm logical flag to plot the heatmap of the similarity matrix.
#'Default is \code{FALSE}.
#'@param nbsim_densities the number of simulated observations to be used
#'to plot the density lines of the clusters in the point estimate partition plot
#'@param ... further arguments passed to or from other methods
#'@author Boris Hejblum
#'@export
print.summaryDPMMclust <- function(x,...){

  cat(class(x), "object with", x$nb_mcmcit, "MCMC iterations:\n")
  cat(rep("-",40),"\n", sep="")
  cat(rep("-",40),"\n", sep="")
  cat("Burnin =", x$burnin, "iterations\n\n")
  cat("Point estimate of the partition with a ", x$clust_distrib," mixture:\n", sep="")
  t<-table(x$point_estim$c_est)
  t2print <- paste(formatC(as.vector(t/length(x$point_estim$c_est)),
                           digits=2),"%", sep="")
  names(t2print) <- names(t)
  print(t2print, quote=F)
  cat("\nLoss of the point estimate partition:", x$loss, "\n")
  cat(rep("-",40),"\n", sep="")

}


#'@rdname methods.summaryDPMMclust
#'@param gg.add a list of instructions to add to the ggplot2 instruction.
#'See \code{\link[ggplot2]{+.gg}}. Default is \code{list(theme())}, which adds nothing to the plot.
#'@importFrom stats dist
#'@importFrom grDevices colorRampPalette
#'@export
plot.summaryDPMMclust <- function(x, hm=FALSE, nbsim_densities=5000, gg.add=list(theme_bw()),...){

  if(length(x$logposterior_list[[1]])>1){
    plot_ConvDPM(x, shift=x$burnin, thin=x$thin)
  }

  ind <- x$index_estim

  cat("Plotting point estimate (may take a few sec)... ")
  if(x$clust_distrib=="gaussian"){
    plot_DPM(z=x$data,
             c=x$point_estim$c_est,
             i=ind+x$burnin,
             alpha=x$alpha[ind],
             U_SS=x$U_SS_list[[ind]],
             ellipses=TRUE,
             gg.add=gg.add,
             ...
    )
  }else if(x$clust_distrib=="skewnorm"){
    plot_DPMsn(z=x$data,
               c=x$point_estim$c_est,
               i=ind+x$burnin,
               alpha=x$alpha[ind],
               U_SS=x$U_SS_list[[ind]],
               ellipses=TRUE,
               gg.add= gg.add,
               nbsim_dens=nbsim_densities,
               ...
    )
  }else if(x$clust_distrib=="skewt"){
    plot_DPMst(z=x$data,
               c=x$point_estim$c_est,
               i=ind+x$burnin,
               alpha=x$alpha[ind],
               U_SS=x$U_SS_list[[ind]],
               ellipses=TRUE,
               gg.add=gg.add,
               nbsim_dens=nbsim_densities,
               nice=TRUE,
               ...
    )
  }
  cat("DONE!\n")

  if(hm){
    cat("Plotting heatmap of similarity (may take a few min)...\n")
    if(is.null(x$point_estim$similarity)){
      stop("In order to plot the similarity matrix, the 'Binder' loss function should be used")
    }
    tree <- fastcluster::hclust(stats::dist(x$point_estim$similarity, method = "euclidean"),
                                method = "complete")
    ord_index <- tree$order
    pheatmap::pheatmap(x$point_estim$similarity[ord_index, ord_index], scale="none", border_color=NA,
                       color=grDevices::colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF",
                                                #"#9ECAE1", "#FEB24C",
                                                "#FD8D3C", "#BD0026", "#800026"))(200),
                       show_rownames=FALSE, show_colnames=FALSE,
                       cluster_rows=FALSE, cluster_cols =TRUE,
                       main="Posterior similarity matrix")
    cat("DONE! Wait for plot rendering...\n")
  }
}
