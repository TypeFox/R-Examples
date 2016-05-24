# $Id: validClimR.R, v1.2.3 2015/08/05 12:00:00 hsbadr EPS JHU            #
#-------------------------------------------------------------------------#
# This function is a part of HiClimR R package.                           #
#-------------------------------------------------------------------------#
#  HISTORY:                                                               #
#-------------------------------------------------------------------------#
#  Version  |  Date      |  Comment   |  Author          |  Email         #
#-------------------------------------------------------------------------#
#           |  May 1992  |  Original  |  F. Murtagh      |                #
#           |  Dec 1996  |  Modified  |  Ross Ihaka      |                #
#           |  Apr 1998  |  Modified  |  F. Leisch       |                #
#           |  Jun 2000  |  Modified  |  F. Leisch       |                #
#-------------------------------------------------------------------------#
#   1.0.0   |  03/07/14  |  HiClimR   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.1   |  03/08/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.2   |  03/09/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.3   |  03/12/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.4   |  03/14/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.5   |  03/18/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.6   |  03/25/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.0.7   |  03/30/14  |  Hybrid    |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.0.8   |  05/06/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.0.9   |  05/07/14  |  CRAN      |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.0   |  05/15/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.1   |  07/14/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.2   |  07/26/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.3   |  08/28/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.4   |  09/01/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.1.5   |  11/12/14  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.1.6   |  03/01/15  |  GitHub    |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
#   1.2.0   |  03/27/15  |  MVC       |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.1   |  05/24/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.2   |  07/21/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#   1.2.3   |  08/05/15  |  Updated   |  Hamada S. Badr  |  badr@jhu.edu  #
#-------------------------------------------------------------------------#
# COPYRIGHT(C) 2013-2015 Earth and Planetary Sciences (EPS), JHU.         #
#-------------------------------------------------------------------------#
# Function: Validation of Hierarchical Climate Regionalization            #
#-------------------------------------------------------------------------#

validClimR <- function(y = NULL, k = NULL, minSize = 1, alpha = 0.05, 
	verbose = TRUE, plot = FALSE, colPalette = NULL, pch = 15, cex = 1) {

    # Coordinates
    lon <- y$coords[, 1]
    lat <- y$coords[, 2]
    
    # Preprocessed raw or PCA-reconstructed data
    x <- y$data
    
    # Cut tree based on minimum significant inter-regional correlation
    if (is.null(k)) {
    	if (verbose) write("---> Cutting tree based on minimum significant correlation...", "")
        # Check clustering method
        if (y$method != "regional" && is.null(y$treeH)) {
            write("--->\t WARNING: objective tree cut is supported only for regional linkage method!", 
                "")
            write(paste("--->\t WARNING: ", y$method, "method requires a prespecified number of clusters!"), 
                "")
        } else {
            if (is.null(y$treeH)) {
                cutHight <- y$height
            } else {
                cutHight <- y$treeH$height
            }
            
            # multi-variate clustering
            nvars <- y$nvars
            missVal <- y$missVal
            n.missVal <- 0
            if (length(missVal) > 0) {
            	for (nvar in 1:length(missVal)) {
            		n.missVal <- n.missVal + length(missVal[[nvar]])
            	}
            }

            # Minimum significant correlation coefficient
			if (verbose) write("---> Computing minimum significant correlation coefficient...", "")
            # for sample size of n years
            #nn <- dim(x)[2] - length(y$missVal)
            nn <- (dim(x)[2] - n.missVal) / nvars
            RsMin <- minSigCor(n = nn, alpha = alpha, r = seq(0, 1, by = 1e-06))$cor
            
            k <- ifelse(length(which(1 - cutHight < RsMin)) > 0, (length(cutHight) - 
                min(which(1 - cutHight < RsMin)) + 1), 2)
            
            # Set minimum k = 2, for objective tree cutting
            k <- ifelse(k < 2, 2, k)
        }
    } else {
        alpha <- NULL
        RsMin <- NULL
    }
    
    index <- list()
    if (!is.null(k)) {
        # Tree cut
        if (is.null(y$treeH)) {
            cutTree <- cutree(y, k = k)
        } else {
            # The reconstructed upper part tree
            if (verbose) write("---> Retrieving the reconstructed upper-part tree...", "")
            yH <- y$treeH
            
            cutTreeH <- cutree(yH, k = k)
            cutTree0 <- cutree(y, k = length(cutTreeH))
            
            cutTree <- cutTree0 + NA
            for (ik in 1:k) {
                cutTree[which(cutTree0 %in% as.integer(names(which(cutTreeH == 
                  ik))))] <- ik
            }
        }
        # table(cutTree)
        
        # Check cluster size
        clustFlag <- rep(1, k)
        if (minSize > 1) {
            clustFlag[which(table(cutTree) < minSize)] <- NA
        }
        
        # Region Means
        if (verbose) write("---> Computing cluster means...", "")
        RM <- t(apply(x, 2, function(r) tapply(r, cutTree, mean)))
        
        # Correlation between Region Means
		if (verbose) write("---> Computing inter-cluster correlations...", "")
        RMcor <- t(fastCor(RM, upperTri = TRUE, verbose = verbose) * clustFlag) * clustFlag
        #RMcor[lower.tri(RMcor, diag = TRUE)] <- NA
        
        # Correlation between Region Means and Region Members
		if (verbose) write("---> Computing intra-cluster correlations...", "")
        Rcor <- cor(RM, t(x))
        
        # Average Correlation between Region Means and Members
        RcorAvg <- t(apply(Rcor, 1, function(r) tapply(r, cutTree, mean)) * 
            clustFlag)
        
        clustFlag[is.na(clustFlag)] <- 0

        if (verbose) write("---> Computing summary statistics...", "")
        index$cutLevel <- c(alpha, RsMin)
        index$clustMean <- RM
        index$clustSize <- table(cutTree)
        
        index$clustFlag <- clustFlag
        names(index$clustFlag) <- 1:length(clustFlag)
        
        index$interCor <- RMcor[!is.na(RMcor)]
        i1.interCor <- which(RMcor %in% index$interCor)
        i2.interCor <- grid2D(c(1:nrow(RMcor)), c(1:ncol(RMcor)))
        names(index$interCor) <- paste(i2.interCor$lat[c(i1.interCor)], 
            i2.interCor$lon[c(i1.interCor)], sep = " & ")
        
        index$intraCor <- diag(RcorAvg)[!is.na(diag(RcorAvg))]
        
        index$diffCor <- index$intraCor - max(index$interCor)
        
        index$statSum <- cbind(summary(index$interCor, digits = 7), summary(index$intraCor, 
            digits = 7), summary(index$diffCor, digits = 7))
        colnames(index$statSum) <- c("interCor", "intraCor", "diffCor")
        
        # Correct average intra-cluster correlation by cluster size
        index$statSum[4, 2] <- sum(index$intraCor * index$clustSize[which(clustFlag == 
            1)])/sum(index$clustSize[which(clustFlag == 1)])
        
        # Oredered regions vector for the selected regions
        ks <- sum(index$clustFlag)
        Regions <- rep(NA, length(lon))
        # Regions[gMask] <- NA
        if (is.null(y$mask)) {
            Regions <- cutTree
        } else {
            Regions[-y$mask] <- cutTree
        }
        Regions[which(Regions %in% which(index$clustFlag != 1))] <- NA
        for (i in 1:ks) {
            Regions[which(Regions == which(index$clustFlag == 1)[i])] <- i
        }
        index$region <- Regions
        index$regionID <- (1:sum(index$clustFlag)) * index$clustFlag[-which(index$clustFlag != 
            1)]
        
        if (plot) {
            if (verbose) write("Generating region map...", "")
            if (is.null(colPalette)) {
                # colPalette <- colorRampPalette(c('#00007F', 'blue', '#007FFF',
                # 'cyan', '#7FFF7F', 'yellow', '#FF7F00', 'red', '#7F0000'))
                colPalette <- colorRampPalette(c("#8DD3C7", "#FFFFB3", 
                  "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", 
                  "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"))
            }
            dev.new()
            Longitude <- y$coords[, 1]
            Latitude <- y$coords[, 2]
            plot(Longitude, Latitude, col = colPalette(max(Regions, na.rm = TRUE))[Regions], 
                pch = pch, cex = cex)
        }
        
        class(index) <- "HiClimR"
    }
    
    #gc()
    return(index)
}
