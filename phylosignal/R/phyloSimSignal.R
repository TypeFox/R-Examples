

#' Phylogenetic signal estimation as a fraction of a Brownian Motion process.
#'
#' This function crosses traits data with simulations to estimate the phylogenetic signal as a
#' fraction of a Brownian Motion process. This is experimental.
#'
#' @param p4d a \code{phylo4d} object.
#' @param phylosim a \code{phylosim} object.
#' @param quantiles a vector of two numeric values between 0 and 1 giving the minimum
#' and the maximum quantiles to estimate statistics fluctuations.
#'
#'@return An object of class \code{phylosimsignal} with a \code{print} and a \code{plot} method.
#'
#'@seealso \code{\link{phyloSignal}}, \code{\link{phyloSim}}.
#'
#'@export
phyloSimSignal <- function(p4d, phylosim, quantiles = c(0.05, 0.95)){
  
  stat.signal <- phyloSignal(p4d, methods = "all", reps = 999)
  stat.signal$stat <- stat.signal$stat[,phylosim$stat.names]
  stat.signal$pvalue <- stat.signal$pvalue[,phylosim$stat.names]
  
  sim.stat.mean <- apply(phylosim$sim.stat, c(1, 2), mean)
  sim.stat.qmin <- apply(phylosim$sim.stat, c(1, 2), quantile, probs = min(quantiles))
  sim.stat.qmax <- apply(phylosim$sim.stat, c(1, 2), quantile, probs = max(quantiles))
  
  nstat <- length(phylosim$stat.names)
  ntrait <- nrow(stat.signal$stat)
  
  res.mean <- res.qmin <- res.qmax <- data.frame(matrix(ncol=nstat, nrow=ntrait))
  
  names(res.mean) <- names(res.qmin) <-names(res.qmax) <- phylosim$stat.names
  rownames(res.mean) <- rownames(res.qmin) <- rownames(res.qmax) <- rownames(stat.signal$stat)
  
  for(i in 1:ntrait){
    for(j in 1:nstat){
      res.mean[i,j] <- which.min(abs(sim.stat.mean[,j] - stat.signal$stat[i,j]))-1
      res.qmin[i,j] <- which.min(abs(sim.stat.qmin[,j] - stat.signal$stat[i,j]))-1
      res.qmax[i,j] <- which.min(abs(sim.stat.qmax[,j] - stat.signal$stat[i,j]))-1
    }
  }
  res <- list(signal.mean=res.mean, signal.qmin=res.qmin, signal.qmax=res.qmax,
              stat.signal=stat.signal, quantiles=quantiles, tree=p4d, phylosim=phylosim)
  class(res) <- "phylosimsignal"
  return(res)
}


#' Print signal estimation as a fraction of a Brownian Motion process.
#' 
#' @param x an object of class \code{phylosimsignal}.
#' @param ... further arguments to be passed to or from other methods.
#' 
#' @method print phylosimsignal
#' @export
print.phylosimsignal <- function(x, ...){
  
  res <- paste(as.matrix(x$signal.mean), " (",
               as.matrix(x$signal.qmax), " - ",
               as.matrix(x$signal.qmin),")", sep="")
  res <- matrix(res, ncol=ncol(x$signal.mean), nrow=nrow(x$signal.mean))
  res <- data.frame(res)
  names(res) <- names(x$signal.mean)
  rownames(res) <- rownames(x$signal.mean)
  print(res)
}


#' Plot signal estimation as a fraction of a Brownian Motion process.
#' 
#' @param x an object of class \code{phylosimsignal}.
#' @param methods a character vector giving the methods
#' (included in the \code{phylosimsignal} object) to plot.
#' @param traits a character vector giving the traits
#' (included in the \code{phylosimsignal} object) to plot.
#' @param stacked.methods If different methods have been used, should they
#' be plotted on the same graphic (\code{TRUE}) or not (\code{FALSE}, default).
#' @param stacked.traits If different traits have been used, should they
#' be plotted on the same graphic (\code{TRUE}) or not (\code{FALSE}, default).
#' @param print.quantiles logical stating whether quantiles should be plotted.
#' @param col a vector of colors for the different methods.
#' @param legend a logical. If \code{stacked.methods} is set to \code{TRUE},
#' should a legend be printed to differentiate the different methods?
#' @param ... further arguments to be passed to or from other methods.
#' 
#' @method plot phylosimsignal
#' @export
plot.phylosimsignal <- function(x, methods=NULL, traits=NULL,
                                stacked.methods=FALSE, stacked.traits=FALSE,
                                print.quantiles=TRUE, col=1:5, legend=TRUE, ...){
  
  x <- subsetPhyloSimSignal(x, methods = methods, traits = traits)
  sim.stat <- x$phylosim$sim.stat
  nstat <- dim(sim.stat)[2]
  ntrait <- nrow(x$signal.mean)
  
  sim.stat.mean <- apply(sim.stat, c(1, 2), mean, na.rm=T)
  xlim <- c(0, 100)
  ylim <- c(min(sim.stat.mean), max(sim.stat.mean))
  
  if(print.quantiles){
    sim.stat.qmin <- apply(sim.stat, c(1, 2), quantile, probs=min(x$quantiles), na.rm=T)
    sim.stat.qmax <- apply(sim.stat, c(1, 2), quantile, probs=max(x$quantiles), na.rm=T)
    ylim <- c(min(sim.stat.qmin), max(sim.stat.qmax))
  }
  
  gcol <- rgb(t(col2rgb(col)), alpha=40, maxColorValue=255)
  save.par <- par()$mfrow
  par(mfrow=c(1, 1))  
  
  if(stacked.traits){
    if(stacked.methods){
      plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab="Brownian Motion (%)", ylab="Statistics values")
      if(legend){
        legend(x="topleft", legend=x$phylosim$stat.names, lty=c(1,1), col=col)
      }
      if(print.quantiles){
        for(i in 1:nstat){
          polygon(c(0:100, 100:0), c(sim.stat.qmin[,i], rev(sim.stat.qmax[,i])), border=NA, col=gcol[i])
        }
      }
      for(i in 1:nstat){
        lines(0:100, sim.stat.mean[,i], col=seq(1,5)[i], lwd=2)
        for(j in 1:ntrait){
          lines(c(-1, x$signal.mean[j,i]), c(x$stat.signal$stat[j,i], x$stat.signal$stat[j,i]),
                lty="dashed", col=col[i])
          lines(c(x$signal.mean[j,i], x$signal.mean[j,i]), c(-100, x$stat.signal$stat[j,i]),
                lty="dashed", col=col[i])
        }
      }
    } else {
      par(mfrow=c(1, nstat))
      for(i in 1:nstat){
        plot(0, 0, type="n", xlim=xlim, ylim=ylim,
             xlab="Brownian Motion (%)", ylab="Statistic values", main=x$phylosim$stat.names[i])
        if(print.quantiles){
          polygon(c(0:100, 100:0), c(sim.stat.qmin[,i], rev(sim.stat.qmax[,i])), border=NA, col=gcol[i])
        }
        lines(0:100, sim.stat.mean[,i], col=col[i], lwd=2)
        for(j in 1:ntrait){
          lines(c(-1, x$signal.mean[j,i]), c(x$stat.signal$stat[j,i], x$stat.signal$stat[j,i]),
                lty="dashed", col=col[i])
          lines(c(x$signal.mean[j,i], x$signal.mean[j,i]), c(-100, x$stat.signal$stat[j,i]),
                lty="dashed", col=col[i])
        }
      }
      par(mfrow=c(1, 1))
    }
    
  } else { # IF TRAITS ARE NOT STACKED
    
    if(stacked.methods){
      for(j in 1:ntrait){
        plot(0, 0, type="n", xlim=xlim, ylim=ylim,
             xlab="Brownian Motion (%)", ylab="Statistics values", main=rownames(x$signal.mean)[j])
        if(legend){
          legend(x="topleft", legend=x$phylosim$stat.names, lty=c(1,1), col=col)
        }
        if(print.quantiles){
          for(i in 1:nstat){
            polygon(c(0:100, 100:0), c(sim.stat.qmin[,i], rev(sim.stat.qmax[,i])), border=NA, col=gcol[i])
          }
        }
        for(i in 1:nstat){
          lines(0:100, sim.stat.mean[,i], col=seq(1,5)[i], lwd=2)
          lines(c(-1, x$signal.mean[j,i]), c(x$stat.signal$stat[j,i], x$stat.signal$stat[j,i]),
                lty="dashed", col=col[i])
          lines(c(x$signal.mean[j,i], x$signal.mean[j,i]), c(-100, x$stat.signal$stat[j,i]),
                lty="dashed", col=col[i])
        }
      }
    } else {
      par(mfrow=c(1, nstat))
      
      for(j in 1:ntrait){
        for(i in 1:nstat){
          main.t <- paste(rownames(x$signal.mean)[j], " (", x$phylosim$stat.names[i],")", sep="")
          plot(0, 0, type="n", xlim=xlim, ylim=ylim,
               xlab="Brownian Motion (%)", ylab="Statistic values", main=main.t)
          if(print.quantiles){
            polygon(c(0:100, 100:0), c(sim.stat.qmin[,i], rev(sim.stat.qmax[,i])), border=NA, col=gcol[i])
          }
          lines(0:100, sim.stat.mean[,i], col=col[i], lwd=2)
          lines(c(-1, x$signal.mean[j,i]), c(x$stat.signal$stat[j,i], x$stat.signal$stat[j,i]),
                lty="dashed", col=col[i])
          lines(c(x$signal.mean[j,i], x$signal.mean[j,i]), c(-100, x$stat.signal$stat[j,i]),
                lty="dashed", col=col[i])
        }
      }
    }
  }
  par(mfrow=save.par)
}