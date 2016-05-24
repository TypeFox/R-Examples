#' @title Plot method for an object of class \code{mbl}
#' @description  Plots the content of an object of class \code{mbl}
#' @aliases plot.mbl
#' @usage \method{plot}{mbl}(x, g = c("validation", "pca"), param = "rmse", pcs = c(1,2), ...)
#' @param x an object of class \code{mbl} (as returned by \code{mbl}). 
#' @param g a character \code{vector} indicating what results shall be plotted. Options are: "validation" (for plotting the validation results) and/or "pca" (for plotting the principal components).
#' @param param one of the following options "rmse", "st.rmse" or "r2". The respective validation statistic is then plotted. It is only available if \code{"validation"} is specified in the \code{g} argument. 
#' @param pcs a vector of length one or two indicating the principal components to be plotted. Default is \code{c(1, 2)}. It is only available if \code{"pca"} is specified in the \code{g} argument. 
#' @param ... some arguments to be passed to methods as graphical parameters.
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @seealso \code{\link{mbl}}
#' @examples
#' \dontrun{
#' require(prospectr)
#' 
#' data(NIRsoil)
#' 
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' Xu <- Xu[!is.na(Yu),]
#' Yu <- Yu[!is.na(Yu)]
#' 
#' Xr <- Xr[!is.na(Yr),]
#' Yr <- Yr[!is.na(Yr)] 
#'
#' ctrl <- mblControl(sm = "cor", ws = 51, 
#'                    pcSelection = list("cumvar", 0.999), 
#'                    valMethod = c("NNv"), 
#'                    scaled = TRUE, center = TRUE)
#'
#' ex1 <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
#'            mblCtrl = ctrl,
#'            dissUsage = "none", 
#'            k = seq(30, 250, 30), 
#'            method = "wapls1",
#'            plsC = c(7, 20))
#'
#' plot(ex1)
#' }
#' @export
###########################################################################
## History:
## 23.04.2014 Leo     Plot function when the data is not centred now 
##                    draws the circles around the actual centre
## 12.04.2015         When the circle was plotted there was an small
##                    gap in it. This was fixed by modifiying the pntCirc 
##                    function




plot.mbl <- function(x, 
                     g = c("validation", "pca"), 
                     param = "rmse", 
                     pcs = c(1, 2), ...){
  
  op.o <- par()$ask
  if(length(g) != 1){
    one.fig <- prod(par("mfcol")) == 1
    op <- par(ask = T)
    on.exit(par(op))
  }
  
  in.call <- match.call()
  if(!is.null(in.call$main))
    main <- in.call$main
  else
    main <- "Memory-based learning results"
  if(!is.null(in.call$col.axis))
    col.axis <- in.call$col.axis
  else
    col.axis <- grey(0.3)
  if(!is.null(in.call$pch))
    pch <- in.call$pch
  else
    pch <- 16
  
  object <- x
  #pm <- par()$mfrow
  
  if("validation" %in% g)
  {
    col <- NULL
    if(!is.null(object$nnValStats)){
      nnValStats <- cbind(object$nnValStats, val = "NNv")
      col <- c(col, "dodgerblue")
    } else {nnValStats <- NULL}
    if(!is.null(object$localCrossValStats)){
      localCrossValStats <- cbind(object$localCrossValStats, r2 = NA, val = "loc_crossval")
      col <- c(col, "green4")
    } else {localCrossValStats <- NULL}
    if(!is.null(object$YuPredictionStats)){
      YuPredictionStats <- cbind(object$YuPredictionStats, val = "Yu prediction")
      col <- c(col, "red")
    } else {YuPredictionStats <- NULL}
    
    tpl <- rbind(nnValStats, localCrossValStats, YuPredictionStats)
    
    if(is.null(tpl)){
      message("No validation results to plot")
    }else{
      #par(mfrow = c(1, length(g)))
      dtn <- colnames(tpl) 
      opt <- c("rmse", "st.rmse", "r2")
      dt <- !is.element(dtn, opt[!is.element(opt, param)])
      toPlot <- reshape(tpl[,dt], timevar = "val", idvar = "k", direction = "wide")
      
      if(param == "r2"){
        toPlot <- toPlot[,!colnames(toPlot) == "r2.loc_crossval"]
        col <- col[!col == "green4"]
      }
      if(is.null(in.call$main)){
        matplot(toPlot[,1], toPlot[,-1], 
                type="b", xlab = dtn[1], 
                ylab = param, pch = 1, ylim = c(min(toPlot[,-1]), 1.1 * max(toPlot[,-1])),
                main = main,
                col = col,
                col.axis = col.axis,...)
      }else{
        matplot(toPlot[,1], toPlot[,-1], 
                type="b", xlab = dtn[1], 
                ylab = param, pch = 1, ylim = c(min(toPlot[,-1]), 1.1 * max(toPlot[,-1])),
                col = col,
                col.axis = col.axis,...) 
      }
      grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
           lwd = par("lwd"), equilogs = TRUE)
      mtext("Validation results", col = grey(0.3))
      # Adding a legend
      legend("topright", legend = colnames(toPlot[,-1,drop = FALSE]), pch = 1,
             col = col, border = "red", box.lty = 3, box.col = "grey")
      #is.element("ggplot2", installed.packages()[,"Package"])
    }
  }
  
  if("pca" %in% g)
  {
    if(object$pcAnalysis$n.componentsUsed == 1)
    {
      rng <- range(object$pcAnalysis$scores_Xr[,pcs[1]], object$pcAnalysis$scores_Xu[,pcs[1]])
      rng <- 1.2 * c(-max(abs(rng)), max(abs(rng)))
      tp <- c(object$pcAnalysis$scores_Xr[,1], object$pcAnalysis$scores_Xu[,1])
      tp <- data.frame(index = 1:length(tp), tp = tp, set = c(rep("Xr", length(object$pcAnalysis$scores_Xr[,1])), rep("Xu", length(object$pcAnalysis$scores_Xu[,1]))))
      tp <- tp[order(tp$tp),]
      tp$index <- 1: length(tp$index)
      plot(tp[tp$set == "Xr",1:2], ylim = rng, 
           col = rainbow(1, s = 1, v = 0, alpha = 0.3), 
           pch = pch, ylab = "pc1 (standardized)", xlab = "index (ordered values)",
           col.axis = col.axis, ...)   
      
      mtext("Prinipal component analyisis", col = grey(0.3))
      points(tp[tp$set == "Xu",1:2], xlim = rng, ylim = rng, col = heat.colors(1, alpha = 0.6), pch = 16)
      grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
           lwd = par("lwd"), equilogs = TRUE)
      legend("topright", legend = c("Xr", "Xu"),
             col = c(rainbow(1, s = 1, v = 0, alpha = 0.3), heat.colors(1, alpha = 0.4)), pch = 16, cex = 0.8, box.lty = 3, box.col = "grey")
    }else{
      if(object$cntrlParam$center){
        rng <- 1.2 * range(object$pcAnalysis$scores_Xr[,pcs], object$pcAnalysis$scores_Xu[,pcs])
        rng1 <- rng2 <-  c(-max(abs(rng)), max(abs(rng)))
      }else{
        rng1 <- range(object$pcAnalysis$scores_Xr[,pcs[1]], object$pcAnalysis$scores_Xu[,pcs[1]])
        rng2 <- range(object$pcAnalysis$scores_Xr[,pcs[2]], object$pcAnalysis$scores_Xu[,pcs[2]])
      }
      
      xl <- paste(colnames(object$pcAnalysis$scores_Xr[,pcs[1],drop=F]), " (standardized)", sep = "") 
      yl <- paste(colnames(object$pcAnalysis$scores_Xr[,pcs[2],drop=F]), " (standardized)", sep = "") 
      
      plot(object$pcAnalysis$scores_Xr[,pcs], xlab = xl, ylab = yl, xlim = rng1, ylim = rng2, 
           col = rainbow(1, s = 1, v = 0, alpha = 0.3), pch = pch,
           col.axis = col.axis, ...)
      
      mtext("Prinipal component analyisis", col = grey(0.3))
      
      points(object$pcAnalysis$scores_Xu[,pcs], xlim = rng1, ylim = rng2, col = heat.colors(1, alpha = 0.4), pch = pch)
      grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
           lwd = par("lwd"), equilogs = TRUE)
      legend("topright", legend = c("Xr", "Xu"),
             col = c(rainbow(1, s = 1, v = 0, alpha = 0.3), heat.colors(1, alpha = 0.4)), pch = pch, cex = 0.8, box.lty = 3, box.col = "grey")
      
      pntCirc <- function(r) {
        n <- 100
        a <- matrix(0, n + 1, 2)
        for (i in 1:n) {
          pnts <- (c(cos(2 * pi/n * i) * r, sin(2 * pi/n * i) * r))
          a[i, ] <- pnts
        }
        a[i+1, ] <- a[1,]
        return(a)
      }
        
      for(i in 1:floor(max(object$pcAnalysis$scores_Xr[,pcs]))){
        crc <- pntCirc(i)
        crc <- rbind(crc, crc[1,])
        if(!object$cntrlParam$center)
          crc <- sweep(x= crc, FUN = "+", MARGIN = 2, STATS = colMeans(object$pcAnalysis$scores_Xr[,pcs]))
        lines(crc, col = "#9ED400", lty = 5, lwd = 0.5)
      }
    }
  }
  
  par(ask = op.o)
  op <- par(ask = op.o)
  on.exit(par(op.o))
  #par(mfrow = pm)
  #title(main = "Memory-based learning results")
  #dev.flush()
}
