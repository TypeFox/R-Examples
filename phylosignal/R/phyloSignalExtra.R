

#' Computes phylogenetic signal for bootstrapped replicates of a phylogeny.
#'
#' This function computes phylogenetic signal statistics and p-values
#' for bootsrapped replicates a phylogenetic tree and produce boxplots to represent the results.
#' This can be useful to check the impact of phylogenetic reconstruction uncertainty on phylogenetic signal.
#'
#' @param p4d a \code{phylo4d} object.
#' @param multiphylo a \code{multiphylo} object containing bootstrapped trees of \code{p4d}.
#' @param methods a character vector giving the methods to compute phylogenetic signal (see \code{\link{phyloSignal}}).
#' @param reps an integer. The number of repetitions for the estimation of p.values with randomization.
#' @param W an optional matrix of phylogenetic weights to compute Moran's I. By default the matrix
#' is computed with the function \code{\link[adephylo]{proxTips}} with patristic distances.
#' @param pb a logical. Should a progress bar be printed? (default \code{TRUE}).
#' 
#' @details
#' Time consumption can be important if there are many bootraped trees and tested traits.
#' 
#' @return
#' The data generated are returned invisibly as a list.
#' 
#' @export
phyloSignalBS <- function(p4d, multiphylo, methods = c("all", "I", "Cmean", "Lambda", "K", "K.star"), reps = 999, W = NULL, pb = TRUE){
  
  methods <- match.arg(methods, several.ok = TRUE)
  if("all" %in% methods){
    methods <- c("I", "Cmean", "Lambda", "K", "K.star")
  }
  
  n.tree <- length(multiphylo)
  X <- tdata(p4d, type = "tip")
  X <- as.matrix(X)
  
  pb <- txtProgressBar(0, n.tree, style = 3)
  
  signal.bs <- vector("list", n.tree)
  for(i in 1:n.tree){
    p4d.i <- phylo4d(multiphylo[[i]], tip.data = X)
    signal.bs[[i]] <- phyloSignal(p4d.i, methods = methods, reps = reps, W = W)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  par(mfrow=c(1, length(methods)+1))
  for(i in methods){
    signal.bs.s <- t(sapply(signal.bs, function(x) x$stat[, as.character(i)]))
    colnames(signal.bs.s) <- colnames(X)
    signal.bs.p <- t(sapply(signal.bs, function(x) x$pvalue[, as.character(i)]))
    colnames(signal.bs.p) <- colnames(X)
    signal.bs.p <- apply(signal.bs.p, 2, function(x) sum(x<0.05)/length(x)*100)
    
    gcol <- grey.colors(101, 1 , 0)
    boxplot(signal.bs.s, col = gcol[signal.bs.p + 1], horizontal = T, las = 1, xlab = i)
  }
  
  plt.init <- par("plt")
  par(plt = c(par("plt")[1] + 0.05, par("plt")[2] - (length(methods) * -0.075 + 0.575), 0.2, 0.8))
  plot.new()
  breaks <- 0:100
  scale.xlim <- c(0, 1)
  scale.ylim <- c(0, 100)
  plot.window(xlim = scale.xlim, ylim = scale.ylim)
  for(i in 1:100){
    polygon(c(0, 0, 1, 1), c(breaks[i], breaks[i + 1], breaks[i + 1], breaks[i]),
            col = gcol[i], border = NA)
  }
  axis(2, las = 2)
  mtext("% Significant Test", side = 4)
  par(plt = plt.init)
  par(mfrow=c(1, 1))
  
  invisible(signal.bs)
}


#' Computes phylogenetic signal at each internal node of a phylogeny
#'
#' This function computes phylogenetic signal statistics and p-values for a given trait
#' and a given method at each internal node of a phylogenetic tree. 
#'
#' @param p4d a \code{phylo4d} object.
#' @param trait a character string giving the trait to use to compute the signal. By default the first trait is taken from \code{p4d}.
#' @param method a character vector giving the method to use to compute phylogenetic signal
#' (default is "\code{Cmean}"; see \code{\link{phyloSignal}}).
#' @param reps an integer. The number of repetitions for the estimation of p.values with randomization.
#' @param W an optional matrix of phylogenetic weights to compute Moran's I. By default the matrix
#' is computed with the function \code{\link[adephylo]{proxTips}} with patristic distances.
#' 
#' @return
#' A \code{phylo4d} object with phylogenetic signal statistics and p-values as nodes associated data.
#' 
#' @export
phyloSignalINT <- function(p4d, trait = names(tipData(p4d))[1], method = "Cmean", reps = 999, W = NULL){
  
  trait <- match.arg(trait, names(tipData(p4d)), several.ok = FALSE)
  method <- match.arg(method, c("I", "Cmean", "Lambda", "K", "K.star"), several.ok = FALSE)
  
  int.nodes <- (nTips(p4d) + 1):(nTips(p4d) + nNodes(p4d))
  new.data <- matrix(NA, nrow = nNodes(p4d), ncol = 2)
  colnames(new.data) <- c(paste("stat", method, trait, sep = "."), paste("pvalue", method, trait, sep = "."))
  rownames(new.data) <- int.nodes
  
  for(i in int.nodes){
    p4d.i <- subset(p4d, node.subtree = i)
    signal.i <- phyloSignal(p4d.i, methods = method, reps = reps, W = W)
    new.data[as.character(i), 1] <- signal.i$stat[trait, method]
    new.data[as.character(i), 2] <- signal.i$pvalue[trait, method]
  }  
  nodeData(p4d) <- data.frame(nodeData(p4d), new.data)
  return(p4d)
}
