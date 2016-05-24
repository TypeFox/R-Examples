

#' Simulate the behaviour of phylogenetic signal statistics with a given phylogeny
#'
#' This function simulates different phylogenetic signal statistics
#' for a given phylogenetic tree along a gradient of Brownian Motion influence.
#'
#' @param tree a \code{phylo}, \code{phylo4} or \code{phylo4d} object.
#' @param methods a character vector giving the methods to compute phylogenetic signal (see Details).
#' @param nsim a numeric value. Number of simulated traits at each step in the gradient.
#' @param reps a numeric value. Number of repetitions for the estimation of p.values with randomization.
#' @param W an optional matrix of phylogenetic weights to compute Moran's I. By default the matrix
#' is computed with the function \code{\link[adephylo]{proxTips}} with patristic distances.
#' @param model the model to use for traits simulation (only "\code{BM}", default, is available).
#' @param pb a logical. Should a progress bar be printed? (default \code{TRUE}).
#'
#' @details By default, the \code{methods} argument is set to "\code{all}" and all the available methods are used.
#' The user can specify which method(s) to use. Possible values are
#' "\code{I}", "\code{Cmean}", "\code{Lambda}", "\code{K}" and "\code{K.star}",
#' see \code{\link{phyloSignal}} for further details.
#'
#' @return An object of class \code{phylosim}.
#'
#' @seealso \code{\link{phyloSimSignal}}.
#'
#' @examples
#' \dontrun{
#' data(navic)
#' psim <- phyloSim(navic)
#' plot(psim)
#' plot.phylosim(psim, what = "pval", stacked.methods = TRUE)
#'}
#'
#' @export
phyloSim <- function(tree, methods = c("all", "I", "Cmean", "Lambda", "K", "K.star"),
                     nsim = 99, reps = 999, W = NULL, model = "BM", pb = TRUE){
  if (inherits(tree, "phylo4")){
    tree <- as(tree, "phylo")
  }
  if (!inherits(tree, "phylo")){
    stop("x has to be a phylo, phylo4 or phylo4d object")
  }
  methods <- match.arg(methods, several.ok = TRUE)
  pb <- txtProgressBar(0, nsim, style = 3)
  sim.list <- vector("list", nsim)
  for(i in 1:nsim){
    sim.list[[i]] <- phyloSignal(rTraitContWeight(tree, model = model,
                                                  weight = seq(from = 0, to = 1, by = 0.01),
                                                  as.p4d = TRUE),
                                 methods = methods,
                                 reps=reps,
                                 W = W)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  sim.stat <- lapply(sim.list, function(x) x[["stat"]])
  stat.names <- colnames(sim.stat[[1]])
  sim.stat <- array(unlist(sim.stat),
                    dim = c(nrow(sim.stat[[1]]), ncol(sim.stat[[1]]), length(sim.stat)))
  sim.pvalue <- lapply(sim.list, function(x) x[["pvalue"]])
  sim.pvalue <- array(unlist(sim.pvalue),
                      dim = c(nrow(sim.pvalue[[1]]), ncol(sim.pvalue[[1]]), length(sim.pvalue)))
  res <- list(tree = tree,
              sim.stat = sim.stat,
              sim.pvalue = sim.pvalue,
              stat.names = stat.names,
              model = model)
  class(res) <- "phylosim"
  return(res)
}



#' Plot \code{phylosim} object
#'
#' This function plots a \code{phylosim} object to visualize the behaviour of
#' phylogenetic signal statistics for a given phylogenetic tree
#'
#' @param x a \code{phylosim} object.
#' @param what what to represent on the plot.
#' Can be the statistics used to measure the signal ("\code{stat}") or the p-values ("\code{pval}").
#' @param stacked.methods a logical. If different methods have been used, should they
#' be plotted on the same graphic (\code{TRUE}) or not (\code{FALSE}, default).
#' @param quantiles a vector of two numeric values between 0 and 1
#' giving the minimum and the maximum quantiles to plot.
#' Set to \code{NULL} to not plot quantiles.
#' @param col a vector of colors for the different methods.
#' @param legend a logical. If \code{stacked.methods} is set to \code{TRUE},
#' should a legend be printed to differentiate the different methods?
#' @param ... further arguments to be passed to or from other methods.
#'
#'@seealso \code{\link{phyloSim}}.
#'
#' @examples
#' \dontrun{
#' data(navic)
#' psim <- phyloSim(navic)
#' plot(psim)
#' plot.phylosim(psim, what = "pval", stacked.methods = TRUE)
#'}
#'@method plot phylosim
#'@export
plot.phylosim <- function(x, what = c("stat", "pval"), stacked.methods = FALSE,
                          quantiles = c(0.05, 0.95), col = 1:5, legend = TRUE, ...){
  what <- match.arg(what)
  if(what == "stat"){
    sim.stat <- x$sim.stat
    nstat <- dim(sim.stat)[2]
    
    sim.stat.mean <- apply(sim.stat, c(1, 2), mean, na.rm=T)
    ylim <- c(min(sim.stat.mean, na.rm = TRUE), max(sim.stat.mean, na.rm = TRUE))
    xlim <- c(0, 100)
    
    if(!is.null(quantiles)){
      sim.stat.qmin <- apply(sim.stat, c(1, 2), quantile, probs=min(quantiles), na.rm = TRUE)
      sim.stat.qmax <- apply(sim.stat, c(1, 2), quantile, probs=max(quantiles), na.rm = TRUE)
      ylim <- c(min(sim.stat.qmin, na.rm = TRUE), max(sim.stat.qmax, na.rm = TRUE))
    }
    
    gcol <- rgb(t(col2rgb(col)), alpha=40, maxColorValue=255)
    
    if(stacked.methods){
      plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab="Brownian Motion (%)", ylab="Statistics values")
      if(legend){
        legend(x="topleft", legend=x$stat.names, lty=c(1,1), col=col)
      }
      if(!is.null(quantiles)){
        for(i in 1:nstat){
          polygon(c(0:100, 100:0), c(sim.stat.qmin[,i], rev(sim.stat.qmax[,i])), border=NA, col=gcol[i])
        }
      }
      for(i in 1:nstat){
        lines(0:100, sim.stat.mean[,i], col=col[i], lwd=2)
      }
    } else {
      par(mfrow=c(1, nstat))
      for(i in 1:nstat){
        plot(0, 0, type="n", xlim=xlim, ylim=ylim,
             xlab="Brownian Motion (%)", ylab="Statistic values", main=x$stat.names[i])
        if(!is.null(quantiles)){
          polygon(c(0:100, 100:0), c(sim.stat.qmin[,i], rev(sim.stat.qmax[,i])), border=NA, col=gcol[i])
        }
        lines(0:100, sim.stat.mean[,i], col=col[i], lwd=2)
      }
      par(mfrow=c(1, 1))
    }
    
  }
  
  if(what=="pval"){    # PLOT PVALUES

    sim.pvalue <- x$sim.pvalue
    nstat <- dim(sim.pvalue)[2]
    
    sim.pvalue.frq <- apply(sim.pvalue, c(1, 2), function(x) sum(na.omit(x)<0.05)/length(na.omit(x)))
    ylim <- c(0, 1)
    xlim <- c(0, 100)
      
    gcol <- rgb(t(col2rgb(col)), alpha=40, maxColorValue=255)
    
    if(stacked.methods){
      plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab="Brownian Motion (%)", ylab="Freq. pvalues < 0.05")
      if(legend){
        legend(x="topleft", legend=x$stat.names, lty=c(1,1), col=col)
      }
      for(i in 1:nstat){
        lines(0:100, sim.pvalue.frq[,i], col=col[i], lwd=2)
      }
    } else {
      par(mfrow=c(1, nstat))
      for(i in 1:nstat){
        plot(0, 0, type="n", xlim=xlim, ylim=ylim,
             xlab="Brownian Motion (%)", ylab="Freq. pvalues < 0.05", main=x$stat.names[i])
        lines(0:100, sim.pvalue.frq[,i], col=col[i], lwd=2)
      }
      par(mfrow=c(1, 1))
    }
  }
}