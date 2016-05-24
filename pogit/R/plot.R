#' Plot an object of class \code{pogit}
#' 
#' This function provides traceplots, autocorrelation plots and density plots 
#' of the MCMC samples for an object of class "\code{pogit}" to graphically assess
#' convergence of the MCMC simulations. It also displays the (model averaged) 
#' posterior means and 95\%-HPD intervals for the regression effects. 
#'
#' @param x an object of class \code{pogit}
#' @param type type of plot: "\code{traceplot}" (default) for traceplots of the 
#'  MCMC draws, "\code{acf}" for autocorrelation plots of the MCMC draws, 
#'  "\code{density}" for density plots and "\code{hpd}" to display 
#'  (model averaged) posterior means with 95\%-HPD intervals for the regression 
#'  effects. 
#' @param burnin logical. If \code{TRUE} (default), burn-in draws (as specified 
#'  in \code{x}) are discarded. 
#' @param thin logical. If \code{TRUE} (default), thinning (as specified in 
#'  \code{x}) is considered for diagnostic MCMC plots. 
#' @param lag.max maximum lag for autocorrelation plot; if \code{NULL} (default), 
#'  the default of \code{\link[stats]{acf}} is used. 
#' @param ci logical. If \code{TRUE} (default), the confidence interval in the
#'  autocorrelation plot is shown (see \code{\link[stats]{acf}} for details).
#' @param maxPlots maximum number of plots on a single page; if \code{NULL} (default),
#'  the number of plots dispayed on a single page is specified according to the
#'  used model. 
#' @param ... further arguments (not used)
#' 
#' @author Michaela Dvorzak <m.dvorzak@@gmx.at>
#' @import ggplot2
#' @import stats
#' @import utils
#' @import grDevices
#' @export
#' 
#' @examples
#' ## see examples for pogitBvs, logitBvs, poissonBvs and negbinBvs

plot.pogit <- function(x, type = "traceplot", burnin = TRUE, thin = TRUE,  
                       lag.max = NULL, ci = TRUE, maxPlots = NULL, ...){
  
  stopifnot(class(x) == "pogit")
  
  #### logit / pogit
  if (x$family %in% c("logit", "pogit")){
    if (burnin && thin){
      chainL <- lapply(x$samplesL, thinMCMC, start = x$mcmc$burnin + 1, 
                       x$mcmc)[c("alpha", "thetaAlpha")]
      
    } else if (burnin && !thin){
      chainL <- lapply(x$samplesL, thinMCMC, start = x$mcmc$burnin + 1, 
                       mcmc = list(nmc = x$mcmc$nmc, thin = 1))[c("alpha", "thetaAlpha")]
      
    } else if (!burnin && thin){
      chainL <- lapply(x$samplesL, thinMCMC, start = NULL, 
                       x$mcmc)[c("alpha", "thetaAlpha")]
      
    } else { # if (!burnin && !thin)
      chainL <- lapply(x$samplesL, thinMCMC, start = NULL, 
                       mcmc = list(nmc = x$mcmc$nmc, thin = 1))[c("alpha", "thetaAlpha")]
      
    } 
    if (x$model.logit$ri == 1) chainL$thetaAlpha <- abs(chainL$thetaAlpha)
    chainL <- do.call(cbind, chainL)
  } else chainL <- NULL
  
  
  ##### poisson / pogit
  if (x$family %in% c("poisson", "pogit")){
    if (burnin && thin){
      chainP <- lapply(x$samplesP, thinMCMC, start = x$mcmc$burnin + 1, 
                       x$mcmc)[c("beta", "thetaBeta")]
      
    } else if (burnin && !thin){
      chainP <- lapply(x$samplesP, thinMCMC, start = x$mcmc$burnin + 1, 
                       mcmc = list(nmc=x$mcmc$nmc, thin=1))[c("beta", "thetaBeta")]
      
    } else if (!burnin && thin){
      chainP <- lapply(x$samplesP, thinMCMC, start = NULL, 
                       x$mcmc)[c("beta","thetaBeta")]
      
    } else { # if (!burnin && !thin)
      chainP <- lapply(x$samplesP, thinMCMC, start = NULL, 
                       mcmc = list(nmc = x$mcmc$nmc, thin = 1))[c("beta", "thetaBeta")]
    } 
    if (x$model.pois$ri == 1) chainP$thetaBeta <- abs(chainP$thetaBeta)
    chainP <- do.call(cbind, chainP)
  } else chainP <- NULL
  
  
  #### negbin
  if (x$family == "negbin"){
    if (burnin && thin){
      chainNB <- lapply(x$samplesNB, thinMCMC, start = x$mcmc$burnin + 1, 
                        x$mcmc)[c("beta", "rho")]
      
    } else if (burnin && !thin){
      chainNB <- lapply(x$samplesNB, thinMCMC, start = x$mcmc$burnin + 1, 
                        mcmc = list(nmc=x$mcmc$nmc, thin=1))[c("beta", "rho")]
      
    } else if (!burnin && thin){
      chainNB <- lapply(x$samplesNB, thinMCMC, start = NULL, 
                        x$mcmc)[c("beta","rho")]
      
    } else { # if (!burnin && !thin)
      chainNB <- lapply(x$samplesNB, thinMCMC, start = NULL, 
                        mcmc = list(nmc = x$mcmc$nmc, thin = 1))[c("beta", "rho")]
    } 
    chainNB <- do.call(cbind, chainNB)
  } else chainNB <- NULL
  
  chains <- cbind(chainL, chainP, chainNB)
  bi <- x$mcmc$burnin
  nPlots <- ncol(chains)
  
  
  if (!(type %in% c("traceplot", "acf", "density", "hpd"))){
    stop(paste(strwrap(paste("invalid 'type' argument: this plot type is not 
                             implemented"), exdent = 1), collapse = "\n"))
  }
  
  if (is.null(maxPlots)){
    maxPlots <- min(min(ncol(chainL), ncol(chainP), ncol(chainNB)), 9)
  }
  if (!is.numeric(maxPlots) || maxPlots < 1){
    stop("maximum number of plots per page must be at least 1")
  }
  
  if (type=="acf" && !is.null(lag.max) && lag.max < 0){
    stop("'lag.max' must be at least 0")
  }
  
  st <- seq(1, nPlots, by = maxPlots)
  sp <- c(st[-1] - 1, nPlots)
  nPages <- length(st) 
  
  # plotting layout
  prepLayout <- function(np){
    d2 <- floor(sqrt(np))
    d1 <- ceiling(np/d2)      
    pos <- expand.grid(1:d1, 1:d2)
    layout <- pos[1:np,]
    return(layout)
  } 
  
  # ggplot basic style elements
  ggElems <- theme_bw() + theme(axis.title.x = element_blank(), 
                                axis.title.y = element_blank())
  
  
  if (type %in% c("traceplot", "acf", "density")){
    doPlot <- sapply(1:nPages, function(page){
      layout <- prepLayout(sp[page] - st[page] + 1)
      p.rows <- layout[, 1]
      p.cols <- layout[, 2]
      
      if (type=="traceplot"){
        it <- values <- NULL # Setting the variables to NULL first
        mat <- stack(as.data.frame(chains[, st[page]:sp[page], drop = FALSE]))
        mat$it <- rep(seq_len(nrow(chains)), sp[page] - st[page] + 1)
        mat$ind <- factor(mat$ind, levels = colnames(chains[, st[page]:sp[page], drop=FALSE]))
        
        ggBase <- ggplot(mat, aes(x = it, y = values))
        plotElems <- ggBase + ggElems + geom_path(colour = "lightseagreen") + 
          facet_wrap(~ind, nrow = max(p.rows), ncol = max(p.cols), scales = "free", drop = TRUE)
        
        if (!burnin && !thin){
          plotElems <- plotElems + geom_vline(xintercept = bi, colour = "lightcoral", 
                                              linetype = 2)
        }
        
        format_xaxis_tp <- function(y){
          if (burnin && thin) {
            return(x$mcmc$burnin + (y - 1)*x$mcmc$thin + x$mcmc$thin)
          } else if (burnin && !thin){
            return(x$mcmc$burnin + y)
          } else if (!burnin && thin){
            return((y - 1)*x$mcmc$thin + x$mcmc$thin)
          } else { # if (!burnin && !thin)
            return(y)
          }
        }
        plotElems <- plotElems + scale_x_continuous(labels = format_xaxis_tp)
        
      } else if (type == "acf"){
        getACF <-  apply(chains, MARGIN = 2, function(x){
          return(acf(x, lag.max = lag.max, plot = FALSE)$acf)
        })
        
        mat <- stack(as.data.frame(getACF[, st[page]:sp[page], drop = FALSE]))
        mat$lag <- rep(seq_len(nrow(getACF)) - 1, sp[page] - st[page] + 1)
        mat$ind <- factor(mat$ind, levels = colnames(getACF[, st[page]:sp[page], drop = FALSE]))
        
        pci <- ci > 0
        clim <- if (pci) qnorm((1 + 0.95)/2)/sqrt(nrow(chains)) else c(0, 0)
        ylim <- range(c(-clim, clim, mat$values))
        
        ggBase <- ggplot(mat, aes(x = lag, y = values)) + geom_hline(aes(yintercept = 0))
        if (pci) ggBase <- ggBase + geom_hline(yintercept = c(clim, -clim), linetype = 2, 
                                               size = 0.4, colour = "blue")
        
        plotElems <- ggBase + ggElems + 
          geom_segment(mapping = aes(xend = lag, yend = 0), colour="gray30", na.rm = TRUE) + 
          facet_wrap(~ind, nrow = max(p.rows), ncol = max(p.cols), scales = "free", drop = TRUE)
        
        plotElems <- plotElems + scale_x_continuous(breaks = pretty(unique(mat$lag), 7)) + 
          scale_y_continuous(limits = ylim)
        
      } else if (type == "density"){
        mat <- stack(as.data.frame(chains[, st[page]:sp[page], drop = FALSE]))
        mat$it <- rep(seq_len(nrow(chains)), sp[page] - st[page] + 1)
        mat$ind <- factor(mat$ind, levels = colnames(chains[, st[page]:sp[page], drop = FALSE]))
        
        ggBase <- ggplot(mat, aes(x = values))
        plotElems <- ggBase + ggElems + 
          geom_density(alpha = 0.3, colour = "lightseagreen", fill = "lightseagreen") + 
          facet_wrap(~ind, nrow = max(p.rows), ncol = max(p.cols), scales = "free", drop = TRUE)
      }
      
      if (dev.interactive() && nPages > 1 && page!=1){ 
        oask <- devAskNewPage(ask = TRUE)
        on.exit(devAskNewPage(oask))
      }
      
      print(plotElems)
      #invisible(plotElems)
    })
  } else if (type == "hpd"){
    res   <- as.data.frame(summary(x)$modTable)
    if (x$family == "negbin"){
      ciHpd <- res[-nrow(res), c(1, (ncol(res) - 1):ncol(res))]
    } else ciHpd <- res[, c(1, (ncol(res) - 1):ncol(res))]
    colnames(ciHpd) <- c("est", "hpdl", "hpdu")
    vars <- factor(colnames(cbind(chainP, chainL, chainNB)))
    if (x$family == "negbin"){
      vars <- vars[-nrow(res)]
      vars <- droplevels(vars)
    }
    ciHpd$var <- factor(vars, levels = colnames(cbind(chainP, chainL, chainNB)))
    
    #hh <- gsub("\\.", "[", ciHpd$var)
    #lab <- paste(hh, "]", sep="")
    #ciHpd$lab <- lab
    
    var <- est <- hpdl <- hpdu <- NULL # Setting the variables to NULL first
    ylim <- range(ciHpd[, c(1:3)])
    ggBase <- ggplot(ciHpd, aes(x = var, y = est, ymin = hpdl, ymax = hpdu)) + 
      geom_hline(aes(yintercept = 0))
    plotElems <- ggBase + ggElems + geom_point() + geom_linerange(colour = "lightseagreen")
    plotElems <- plotElems + scale_y_continuous(limits = ylim) + 
      scale_x_discrete(labels = abbreviate)
    if (x$family == "pogit"){
      plotElems <- plotElems + 
        geom_vline(xintercept = x$model.pois$d + x$model.pois$ri + 1.5, linetype = 2)
    }
    print(plotElems)
  }
}



