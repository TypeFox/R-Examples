## COMMUNAL package
## part 2: plot the cluster metrics

################   plotting Functions   ####################################
#Here each measure is scaled across ALL data_range, ks, and algs, meaning the final plot
#is an overall measure

getScaledAverages <- function(data_range, geneRange, algs, metrics, ks){
  #Number of cluster runs does not match number of variables
  stopifnot(length(data_range)==length(geneRange))
  lowerIsBetter <- c("widestgap", "wb.ratio", "max.diameter", "within.cluster.ss",
                     "average.within", "Connectivity", "entropy", "g3")
  
  #4D array to store all measures; will allow scaling
  mult_array <- array(data=NA, 
                      dim=c(length(metrics), length(data_range), 
                                     length(ks), length(algs)),
                      dimnames=list(metrics, geneRange, ks, algs) )
  
  for(j in 1:length(data_range)){
    #get measures for one data_range instance
    measures <- data_range[[j]]$measures[metrics,as.character(ks),algs,drop=F]
    #make sure all measures in group are plotted  so that higher is better
    for (i in 1:length(metrics)) { 
      flip <- ifelse(metrics[i] %in% lowerIsBetter, -1, 1)
      mult_array[i,j, ,] <- flip*as.numeric(measures[i,,]) 
    }
  }; rm(i,j)
  
  mult_array_scaled <- array(NA, dim=dim(mult_array))
  #scale ALL values for single measure (across data range) together
  for(i in 1:length(metrics)){
    mult_array_scaled[i,,,] <- array(data=(mult_array[i,,,]-mean(mult_array[i,,,], na.rm=T))/sd(mult_array[i,,,], na.rm=T), 
                                     dim=dim(mult_array[i,,,])) 
  }
  dimnames(mult_array_scaled) <- dimnames(mult_array)

  #take average over all algs and measures
  average <- data.frame(matrix(NA, nrow=length(ks), ncol=length(data_range)))  
  rownames(average) <- ks; colnames(average) <- geneRange  
  for(j in 1:length(data_range)){
    ## if only one metric passed, then the first dim of array collapes, so need to decrement index
    ## if only one alg passed, 
    ## if only one of each passed, skip apply() and just take mean
    if(length(metrics) > 1) {
      average[,j] <- apply(mult_array_scaled[,j,,], 2, mean, na.rm=T)
    } else if (length(metrics)==1 & length(algs) > 1) {
      average[,j] <- apply(mult_array_scaled[,j,,], 1, mean, na.rm=T)
    } else if (length(metrics)==1 & length(algs) == 1) {
      average[,j] <- mult_array_scaled[,j,,]
    }
  }
  if(any(is.na(average))) {
    print(average); 
    stop( "NAs in averages matrix (see above). Decrease range of K.")
  }
  
  return(average)
}


# Find the highest non-edge 'peak' in a vector (most concave point)
#  or, if no local 'peak', return highest value
findSteepestPeak <- function(x){
    ## get all peaks
    neg_if_peak <- diff(diff(x)>=0)
    if(any(neg_if_peak<0, na.rm=T)){
        ## the further negative, the steeper
        ## for all peaks, pick the steepest    
        diffs <- diff(diff(x))
        x.peak <- which(diffs == min(diffs[neg_if_peak<0]))+1
        ## returns c(position_in_vector, value)
        return(c(x.peak, x[x.peak]))
    } else {
        return(c(which.max(x), x[which.max(x)]))
    }
}


# v 1.1: added max peaks point (light blue) and colorbar
plotRange <- function(averageMtx, geneRange, ks, filename, colorbar=T){
  geneRangeLabels <- factor(unlist(geneRange))
  names(geneRangeLabels) <- unlist(geneRange)
  rgl::open3d()
  for(i in 1:length(geneRangeLabels)){
    rgl::lines3d(x=ks, y=rep(geneRangeLabels[i], length(ks)), z=averageMtx[,i], lwd=3)
    ## connect grid to lines
    rgl::lines3d(x=rep(ks[1], 2), y=rep(geneRangeLabels[i], 2), 
                 z=c(max(averageMtx), averageMtx[1,i]), 
                 col="grey30", lwd=0.3)
  }
  for(i in seq(min(ks), max(ks), 2)){
    rgl::lines3d(x=c(i,i), y=c(1, length(geneRangeLabels)), z=rep(min(averageMtx),2), lwd=0.5)
  }
  steepPeaks <- apply(averageMtx, 2, findSteepestPeak)
  rgl::points3d(x=ks[steepPeaks[1,]], y=geneRangeLabels, z=steepPeaks[2,], col=2, size=9)
  
  highPeaks <- apply(averageMtx, 2, function(x) c(which.max(x), x[which.max(x)]))
  rgl::points3d(x=ks[highPeaks[1,]], y=geneRangeLabels, z=highPeaks[2,], col=5, size=11)
  
  
  z <- data.matrix(averageMtx)
  z <- data.matrix(averageMtx)
  # height color lookup table
  colorlut <- rainbow(length(z),alpha=0, start=length(z)/(4*length(z)) )
  # stretch Z along the length of z
  z.stretch <- (z-min(z))*(length(z)-1)/(max(z)-min(z))+1
  colors <- colorlut[z.stretch]
  
  rgl::surface3d(x=ks, y=1:length(geneRangeLabels), z=z, color=colors, alpha=0.5)
  
  if(colorbar){
    ## make a color bar
    index <- order(unlist(z))
    rgl::lines3d(x=rep(max(ks), length(z)), y=rep(1, length(z)),
                 z=z[index], col=colors[index], lwd=6)
  }
  offset <- "      "
  rgl::box3d()
  rgl::mtext3d(text= "K", edge="x--", line=2)
  rgl::mtext3d(text= "Variables                   ", edge="y-+", line=3)
  rgl::mtext3d(text= paste(paste(rep(offset, 6), collapse=""), "Z-measures"),
               edge="z+-", line=1)
  rgl::axis3d(edge="x--", at=ks, labels=ks);
  rgl::axis3d(edge="y-+", labels=paste(geneRange, offset), at=1:length(geneRange));
  zlabs <- signif(c(min(z), 0, max(z)), 2)
  rgl::axis3d(edge="z+-", labels= zlabs, at=zlabs); 
  if (!is.null(filename)) {
    rgl::rgl.snapshot(filename)
  }
}

plotMeans <- function(averageMtx, algs, metrics, ...){
  plot(apply(averageMtx, 1, mean), main=paste("Mean across all geneRange \n", 
                                              paste(algs, collapse=", "), "\n", 
                                              paste(metrics, collapse=", ")), 
        xlab="K", xaxt='n', ylab="Combined metric Z-score", type="b",  ...)
  axis(1, at=1:dim(averageMtx)[1], labels=rownames(averageMtx) )
}

#Wrapper
plotRange3D <- function(clusRange, ks=NULL, goodAlgs=NULL, goodMeasures=NULL, 
                        filename=NULL, colorbar=T, minSize=3, plot3D=T, ...){
  meas <- clusRange[[1]][[1]]$measures
  if(is.null(ks)) ks <- as.numeric(dimnames(meas)[[2]])
  if(is.null(goodAlgs)) goodAlgs <- dimnames(meas)[[3]]
  if(is.null(goodMeasures)) goodMeasures <- dimnames(meas)[[1]]
    
  average <- getScaledAverages(data_range=clusRange[[1]], geneRange=clusRange[[2]], 
                               algs=goodAlgs, metrics=goodMeasures, ks=ks)
  
  clusterCounts <- t(sapply(clusRange[[1]], function(clusVar){  
    clusters <- count_clusters(clusVar, goodAlgs, minSize=minSize, 
                               verbose=F, returnClusters=T)
    minClusters <- sapply(clusters, function(tmp) colSums(tmp<minSize))
    rowSums(t(t(minClusters)))
  }))
  print("Count of clusters < minSize for the algorithms tested (over all K):\n")
  print(clusterCounts)
  
  if ( (ncol(average) > 1) & plot3D) {
    plotRange(averageMtx=average, clusRange[[2]], ks, filename, colorbar=colorbar)
  }
  plotMeans(averageMtx=average, algs=goodAlgs, metrics=goodMeasures, ...)
  average
}



###########  analysis functions  #####################
count_clusters <- function(COMMUNAL_out, goodAlgs, minSize=3, verbose=T, returnClusters=F){
  out <- lapply(COMMUNAL_out$ks, function(k) {
    clusters <- COMMUNAL_out$getClustering(k)
    clusters <- clusters[, goodAlgs,drop=F]
    tmp <- matrix(0, nrow=k, ncol=ncol(clusters))
    colnames(tmp) <- colnames(clusters)
    rownames(tmp) <- as.character(1:k)
    for(i in 1:ncol(clusters)){
      counts <- table(clusters[,i])
      tmp[names(counts),i] <- counts
    }
    
    if(any(tmp < minSize) & verbose){
      cat(sprintf("\n Warning: at k= %s, some clusters have fewer than %s members (counts shown below):\n",
                  k, minSize))
      print(tmp)
    } 
    
    return(tmp)
  })
  if(returnClusters) return(out)
}


clusCounts <- function(cluster.table){
  ## count the number of samples in each cluster
  ## assumes input in a list from clusterRange
  k <- max(unlist(cluster.table))
  clusTables <- lapply(cluster.table, function(x){
    y <- rep(0, k) 
    y[as.numeric(unique(x))] <- table(x)
    y
  })
  out <- Reduce(cbind, clusTables)
  colnames(out) <- names(cluster.table)
  out
}


monotoneClusterRange <- function(clusRange, goodAlgs=NULL){
  ## test a given measure for monotonicity over a range of K
  is.monotonic <- function(x) !is.unsorted(x, na.rm=T) | !is.unsorted(rev(x), na.rm=T)
  
  if(is.null(goodAlgs)) goodAlgs <- dimnames(clusRange$all.results[[1]]$measures)[[3]]
  pct.monotone <- sapply(clusRange$all.results, function(out){
    apply(out$measures, 1, function(measure){
      measure <- measure[, goodAlgs]
      sum(apply(measure, 2, is.monotonic))/dim(measure)[[2]]
    })
  })
  rowMeans(pct.monotone)
}

meas.flip <- function(meas.alg){
  ## some measures lower is better; rotate these around their mean, so that they can be compared to 
  ## standard measures; input is a matrix with measures in rows, K in columns
  lowerIsBetter <- c("widestgap", "wb.ratio", "max.diameter", "within.cluster.ss",
                     "average.within", "Connectivity", "entropy", "g3")
  lowerIsBetter <- rownames(meas.alg)[rownames(meas.alg) %in% lowerIsBetter]
  meas.alg[lowerIsBetter, ] <- -1*meas.alg[lowerIsBetter, ] + 2*rowMeans(meas.alg[lowerIsBetter, ])
  meas.alg
}


testAlgsMinSize <- function(clusRange, algs="all", minSize=3){
  if("all" %in% algs) algs <- dimnames(clusRange[[1]][[1]]$measures)[[3]]
  
  ## for each vars in clusRange, get cluster sample assignments and count < minSize
  counts <- lapply(1:length(clusRange[[1]]), function(i){
    cat("\n Testing for cluster assignments at ", clusRange[[2]][i], "variables: ")
    tmp <- count_clusters(clusRange[[1]][[i]], goodAlgs=algs, minSize=minSize, 
                          verbose=T, returnClusters=T)
    colSums(Reduce(rbind, tmp)<minSize)
  })
  
  
  ## sum over all varRange
  sumCounts <- colSums(Reduce(rbind, counts))
  
  ks <- as.numeric(dimnames(clusRange$all.results[[1]]$measures)[[2]])
  possible <- length(clusRange[[1]]) * tail(cumsum(ks), 1)
  pctCounts <- round(sumCounts/possible, 3)
  
  cat("Final fraction of algs returning clusters < minSize for all vars and K:\n")
  print(pctCounts)
  return(pctCounts)
}

## wrapper for testAlgsMinSize that filters to get algs < mean(algs)
getGoodAlgs <- function(clusRange, algs="all", minSize=3){
  if("all" %in% algs) algs <- dimnames(clusRange[[1]][[1]]$measures)[[3]]
  algs <- testAlgsMinSize(clusRange, algs=algs, minSize=minSize)
  names(algs)[algs <= mean(algs)] ## mean, not median
}


measuresCorr <- function(clusRange, goodMeasures){
  ## get all measures across all k and vars for each measure
  allMeasures <- Reduce(rbind, lapply(clusRange$all.results, function(vars){
    Reduce(rbind, lapply(1:dim(vars$measures)[3], function(i) t(meas.flip(vars$measures[,,i])) ))
  }))
  
  ## get correlation matrix
  allMeasures <- allMeasures[complete.cases(allMeasures), goodMeasures]
  cor(allMeasures)
}

getNonCorrNonMonoMeasures <- function(clusRange, goodMeasures="all", goodAlgs=NULL, numMeasures=4){
  ## first get counts of monotonicity for all measures
  ## limit to 'good' measures if desired
  ## throw out monotonic greater than mean
  monotone <- monotoneClusterRange(clusRange, goodAlgs=goodAlgs)
  if(!("all" %in% goodMeasures)) monotone <- monotone[goodMeasures]
  monotone <- names(monotone)[monotone <= median(monotone)]
  
  ## measure correlation of measures left
  ## cluster and cut at user-defined K
  monoCorr <- measuresCorr(clusRange, monotone)
  measuresClus <- cutree(hclust(dist(monoCorr)), k=numMeasures)
  
  ## within each cluster of measures, choose the one with lowest
  ## overall correlation
  finalMeasures <- unlist(lapply(unique(measuresClus), function(km){
    measuresClustered <- names(measuresClus)[measuresClus==km]
    measuresClusCorr <- rowSums(monoCorr)[measuresClustered]
    names(measuresClusCorr)[which.min(measuresClusCorr)]
  }))
  return(finalMeasures)
}



