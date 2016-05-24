#  File R/degreedist.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
degreedist.egodata <- function(egodata, freq = FALSE, prob = !freq, 
                               by = NULL, brgmod = FALSE){
  if (class(egodata) != "egodata"){
    stop("The egodata object passed to degreedist.egodata must be of class egodata.")
  }
  color <- "#83B6E1"
  beside <- TRUE
  ylabel <- "Frequency"
  egoIDcol <- egodata$egoIDcol
  degtable <- rep(0, nrow(egodata$egos))
  degtable[as.numeric(names(table(egodata$alters[egoIDcol])))] <- table(egodata$alters[egoIDcol])
  degtable.wt <- degtable * egodata$egoWt
  maxdeg <- max(degtable.wt)
  deg.ego <- summary(egodata ~ degree(0:maxdeg, by = by))
  names(deg.ego) <- 0:maxdeg
  maxfreq <- max(deg.ego)
  if(!is.null(by)){
    levs <- sort(unique(c(egodata$egos[[by]], egodata$alters[[by]])))
    deg.ego <- matrix(0, nrow = length(levs), ncol = maxdeg + 1)
    rownames(deg.ego) <- levs
    colnames(deg.ego) <- 0:maxdeg
    for(i in 1:length(levs)){
      vals <- table(degtable.wt[egodata$egos[by] == levs[i]])
      toreplace <- as.numeric(names(vals)) + 1
      deg.ego[i, toreplace] <- vals
    }
    ncolors <- dim(deg.ego)[1]
    if(ncolors == 2){
      color <- c("#eff3ff", "#377FBC")
    } else if(ncolors < 10){
      color <- RColorBrewer::brewer.pal(ncolors,"Blues")
    } else if(ncolors >= 10){
      color <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(ncolors)
    }
    
    ltext <- levs
    lfill <- c(color, 0)
    lborder <- c(rep("black", times = ncolors), 0)
    ltitle <- by
    maxfreq <- max(colSums(deg.ego))
  }
  if(prob){
    if(is.null(by)){
      scaledeg <- sum(deg.ego)
      deg.ego <- deg.ego/scaledeg
      maxfreq <- max(deg.ego, na.rm = TRUE)
    } else {
      scaledeg <- rowSums(deg.ego)
      deg.ego <- deg.ego/scaledeg
      deg.ego <- deg.ego
      maxfreq <- max(max(deg.ego, na.rm = TRUE))
      beside <- TRUE
    }
  }


  if(brgmod) {
    brgdraws <- simulate.ergm.ego(ergm.ego(egodata ~ edges), nsim = 50)
    deg.brg <- summary(brgdraws ~ degree(0:maxdeg))
    brgmeans <- apply(deg.brg, MARGIN = 2, FUN = mean)
    brgsd <- apply(deg.brg, MARGIN = 2, FUN = sd)
    upper <- brgmeans + 2 * brgsd
    lower <- brgmeans - 2 * brgsd
    
    if(prob){
      if(is.null(by)){
        brgmeans <- brgmeans/scaledeg
        upper <- upper/scaledeg
        lower <- lower/scaledeg
        ylabel <- "Probability"
      } else {
        upper <- upper/sum(brgmeans)
        lower <- lower/sum(brgmeans)
        brgmeans <- brgmeans/sum(brgmeans)
        ylabel <- "Probability (within attr level)"
      }
      
    }
    maxfreq <- max(maxfreq, upper, na.rm = TRUE)
  }
  
  baraxis <- barplot(deg.ego, xlab = "Degree", ylab = ylabel,
                     col = color, beside = beside, plot = TRUE,
                     ylim = c(0, maxfreq))
  
  if(brgmod){
    baraxis <- if(is.null(by)){
      baraxis - 0.15
    } else {
      colMeans(baraxis)
    }
    points(x = baraxis, y = brgmeans, col = "firebrick",
           lwd = 1, pch = 18, cex = 1.25)
    suppressWarnings(arrows(x0 = baraxis, y0 = upper,
                            x1 = baraxis, y1 = lower,
                            code = 3, length = 0.1, 
                            angle = 90, col = "firebrick"))
  } 
  if(!is.null(by)){
    legend(x="topright", legend = ltext, title = ltitle, fill = lfill, 
           border = lborder, bty = "n")
  }
}


mixingmatrix.egodata <- function(egodata, attrname, rowprob = FALSE){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  levs <- sort(unique(c(egos[[attrname]], alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  
  ties <- merge(egos[c(egoIDcol,attrname)], alters[c(egoIDcol,attrname)], 
                by = egoIDcol, suffixes = c(".ego",".alter"))
  ties$wt <- egodata$egoWt[ties$vertex.names]
  mxmat <- matrix(0, nrow = length(levs), ncol = length(levs))
  
  for(i in 1:length(levs)){
    for(j in 1:length(levs)){
      mxmat[i,j] <- sum(ties$wt[ties[,2]==i & ties[,3]==j])
    }
  }
  dimnames(mxmat) <- list(ego = levs,  
                          alter = levs)
  if(rowprob){
    mxmat <- mxmat/rowSums(mxmat)
    mxmat <- round(mxmat, digits = 3)
  }
  mxmat
}
  
  
