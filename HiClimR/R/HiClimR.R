# $Id: HiClimR.R, v1.2.3 2015/08/05 12:00:00 hsbadr EPS JHU               #
#-------------------------------------------------------------------------#
# This is the main function of                                            #
# HiClimR (Hierarchical Climate Regionalization) R package                #
#-------------------------------------------------------------------------#
#  References:                                                            #
#                                                                         #
#  Badr, H. S., Zaitchik, B. F. and Dezfuli, A. K. (2015).                #
#  A Tool for Hierarchical Climate Regionalization.                       #
#  Earth Science Informatics, In Review.                                  #
#                                                                         #
#  Badr, H. S., Zaitchik, B. F. and Dezfuli, A. K. (2015).                #
#  Hierarchical Climate Regionalization. CRAN,                            #
#  http://cran.r-project.org/package=HiClimR.                             #
#                                                                         #
#  Source Code: https://github.com/hsbadr/HiClimR                         #
#-------------------------------------------------------------------------#
#  Clustering Methods:                                                    #
#                                                                         #
#  0. REGIONAL linakage or minimum inter-regional correlation             #
#  1. WARD's minimum variance or error sum of squares method              #
#  2. SINGLE linkage or nearest neighbor method                           #
#  3. COMPLETE linkage or diameter                                        #
#  4. AVERAGE linkage, group average, or UPGMA method                     #
#  5. MCQUITTY's or WPGMA method                                          #
#  6. MEDIAN, Gower's or WPGMC method                                     #
#  7. CENTROID or UPGMC method                                            #
#-------------------------------------------------------------------------#
# This code is modified by Hamada S. Badr <badr@jhu.edu> from:            #
# File src/library/stats/R/hclust.R                                       #
# Part of the R package, http://www.R-project.org                         #
#                                                                         #
# Copyright(C)  1995-2015  The R Core Team                                #
#                                                                         #
# This program is free software; you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation; either version 2 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# This program is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            #
# GNU General Public License for more details.                            #
#                                                                         #
# A copy of the GNU General Public License is available at                #
# http://www.r-project.org/Licenses                                       #
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
# Function: Hierarchical Climate Regionalization                          #
#-------------------------------------------------------------------------#

HiClimR <- function(x = list(), 
	lon = NULL, lat = NULL, lonStep = 1, latStep = 1, 
    geogMask = FALSE, gMask = NULL, continent = NULL, region = NULL, country = NULL, 
	meanThresh = if(class(x) == "list") vector("list", length(x)) else list(NULL), 
	varThresh = if(class(x) == "list") as.list(rep(0, length(x))) else list(0), 
	detrend = if(class(x) == "list") as.list(rep(FALSE, length(x))) else list(FALSE), 
	standardize = if(class(x) == "list") as.list(rep(FALSE, length(x))) else list(FALSE), 
	weightedVar = if(class(x) == "list") as.list(rep(1, length(x))) else list(1),
	nPC = NULL, method = "ward", hybrid = FALSE, kH = NULL, members = NULL, 
	nSplit = 1, upperTri = TRUE, verbose = TRUE, 
	validClimR = TRUE, rawStats = TRUE, k = NULL, minSize = 1, alpha = 0.05, 
    plot = TRUE, dendrogram = TRUE, colPalette = NULL, hang = -1, labels = FALSE, 
    pch = 15, cex = 1) {

    start.time.proc <- proc.time()
    start.time.sys <- Sys.time()
    if (verbose) write("\nPROCESSING STARTED\n", "")
    
    # Check number of variable
    if (verbose) write("Checking Multi-Variate Clustering (MVC)...", "")
    nvars <- 1
    if (class(x) == "list" ) {
		if (verbose) write("---> x is a list", "")

    	nvars <- length(x)
		mm <- vector(mode="numeric", length=nvars)
    	if (nvars == 0) {
 	       stop("empty data list")
		}
		
		xx <- x[[1]]
		mm[1] <- dim(x[[1]])[2]
    	if (nvars > 1) {
    		if (length(meanThresh) != nvars || 
    			length(varThresh) != nvars || 
    			length(detrend) != nvars ||
    			length(standardize) != nvars ||
    			length(weightedVar) != nvars) {
    			stop("meanThresh, varThresh, detrend, standardize and weightedVar length mismatch")
    		}
    		n1 <- dim(x[[1]])[1]
    		for (nvar in 2:nvars ) {
		    	# Check dimensions for multi-variate analysis
    			if (dim(x[[nvar]])[1] != n1) {
    				stop("matrices in x list should have the same number of rows!")
    			}
				xx <- cbind(xx, x[[nvar]])
				mm[nvar] <- dim(x[[nvar]])[2]
    		}
	    }
		x <- xx
		rm(xx)
    } else if (class(x) == "matrix") {
    	if (verbose) write("---> x is a matrix", "")
    	mm <- dim(x)[2]
    } else {
    	if (verbose) write("---> x should a matrix or a list of matrices", "")
    }
		
	if (verbose) {
	    if (nvars == 1) {
	    	write("---> single-variate clustering: 1 variable", "")
    	} else {
    		write(paste("---> multi-variate clustering:", nvars, "variables"), "")
    	}
    }

	# Check variable weights
	if (nvars > 1) {
	    if (verbose) write("Checking variable weights...", "")
		for (nvar in 1:nvars ) {
			if (verbose) write(paste("---> weight for variable #", nvar, ": ", 
				weightedVar[[nvar]], sep=""), "")
			if (! weightedVar[[nvar]] > 0) {
				stop("variable weights should be positive!")
			}
		}
	} else {
		weightedVar <- 1
	}

    # Coarsening spatial resolution
    if (lonStep > 1 && latStep > 1) {
        if (verbose) write("Coarsening spatial resolution...", "")
        xc <- coarseR(x = x, lon = lon, lat = lat, lonStep = lonStep, latStep = latStep, 
        	verbose = verbose)
    } else {
        xc <- coarseR(x = x, lon = lon, lat = lat, lonStep = 1, latStep = 1, verbose = verbose)
    }
    lon <- xc$lon
    lat <- xc$lat
    x <- xc$x
    rm(xc)

    if (verbose) write("Checking data...", "")

    # Check data dimensions
    if (verbose) write("---> Checking dimensions...", "")
    n <- dim(x)[1]
    m <- dim(x)[2]
    if (is.null(n)) 
        stop("invalid data size")
    # if (is.na(n) || n > 65536L) stop(' size cannot be NA nor exceed
    # 65536')
    if (is.na(n)) 
        stop("size cannot be NA")
    if (n < 2) 
        stop("must have n \u2265 2 objects to cluster")
    
    # Check row names (important if detrending is requested)
    if (verbose) write("---> Checking row names...", "")
    if (is.null(rownames(x))) {
        rownames(x) <- seq(1, n)
    }
    
    # Check column names (important if detrending is requested)
    if (verbose) write("---> Checking column names...", "")
    if (is.null(colnames(x))) {
        colnames(x) <- seq(1, m)
    }
    
    # Mask geographic region
    mask <- NULL

    if (geogMask) {
        if (verbose) write("Geographic masking...", "")

        if (is.null(gMask)) {
            gMask <- geogMask(continent = continent, region = region, country = country, 
                lon = lon, lat = lat, verbose = verbose, plot = FALSE)
        } else {
            if (verbose) write("---> Geographic mask is provided!", "")
        }

		if (length(gMask) > 0 && class(gMask) != "list") {
	        if (min(gMask) >= 1 && max(gMask) <= n) {
	            mask <- union(mask, as.integer(gMask))
	        }
	    }
    }
    
    mm0 <- 0
    xx <- x
    vv <- vector("list", nvars)
    if (verbose) write("Data filtering...", "")
    for (nvar in 1:nvars ) {
    	if (nvars > 1) {
    		if (verbose) write(paste("---> VARIABLE #", nvar,":", sep=""), "")
    	}
    	
    	m <- mm[nvar]
    	x <- xx[,(mm0 + 1):(mm0 + m)]

        if (verbose) write("---> Computing mean for each row...", "")
	    # xmean <- rowMeans(x, na.rm=TRUE)
	    xmean <- rowMeans(x)
        
        # Remove rows with observations mean bellow meanThresh
	    if (!is.null(meanThresh[[nvar]])) {
            if (verbose) write("---> Checking rows with mean bellow meanThresh...", "")
	        meanMask <- which(is.na(xmean) | xmean <= meanThresh[[nvar]])
            if (verbose) write(paste("--->", length(meanMask), "rows found, mean \u2264 ", 
                meanThresh[[nvar]]), "")
	        if (length(meanMask) > 0) {
	            mask <- union(mask, meanMask)
	        }
	    }
	    
	    # Center data (this has no effect on correlations but speedup compuations)
	    x <- x - xmean
        if (verbose) write("---> Computing variance for each row...", "")
	    v <- rowSums(x^2, na.rm = TRUE)
    
	    # Remove rows with near-zero-variance observations
	    if (is.null(varThresh[[nvar]])) {
	        varThresh[[nvar]] <- 0
	    }
	    varMask <- which(is.na(v) | v <= varThresh[[nvar]])
        if (verbose) write("---> Checking rows with near-zero-variance...", "")
        if (verbose) write(paste("--->", length(varMask), "rows found, variance \u2264 ", 
            varThresh[[nvar]]), "")
	    if (length(varMask) > 0) {
        	mask <- union(mask, varMask)
    	}
    	
    	vv[[nvar]] <- v
    	xx[,(mm0 + 1):(mm0 + m)] <- x
    	mm0 <- mm0 + m
    }

    mm0 <- 0
    xxx <- NULL
    missVal <- list()
    if (verbose) write("Data preprocessing...", "")
    for (nvar in 1:nvars ) {
    	if (nvars > 1) {
    		if (verbose) write(paste("---> VARIABLE #", nvar,":", sep=""), "")
    	}
    	
    	m <- mm[nvar]
    	v <- vv[[nvar]]
    	x <- xx[,(mm0 + 1):(mm0 + m)]
        
    	if (verbose) write("---> Applying mask...", "")
    	# Mask data
    	if (length(mask) > 0) {
    	    x <- x[-mask,]
    	    v <- v[-mask]
    	}
    
	    # Remove columns with missing values
    	if (verbose) write("---> Checking columns with missing values...", "")
    	x <- t(na.omit(t(x)))
    	if (verbose && length(attr(x, "na.action")) > 0) {
        	write(paste("--->\t WARNING:", length(attr(x, "na.action")), "columns found with missing values"), 
        	    "")
    	}
        
    	# Detrend data if requested
    	if (detrend[[nvar]]) {
        	if (verbose) write("---> Removing linear trend...", "")
        	x <- x - t(fitted(lm(t(x) ~ as.integer(colnames(x)))))
    	}
    
    	# Standardize data if requested
    	if (standardize[[nvar]]) {
    	    if (verbose) write("---> Standardizing data...", "")
        	x <- x/sqrt(v)
        	# Variance of each variable (object/station)
        	v <- rep(1, n)
        
        	# Standardized data
        	x <- x * sqrt(m - 1)
    	} else {
        	# Variance of each variable (object/station)
        	v <- v / (m - 1)
    	}    
    	# Re-adding the mean for nonstandardized data (July 26, 2014)
    	if (!standardize[[nvar]]) {
    	    if (length(mask) > 0) {
        	    x <- x + xmean[-mask]
        	} else {
            	x <- x + xmean
        	}
    	}
    	
    	xxx <- cbind(xxx, weightedVar[[nvar]] * x)
    	missVal[[nvar]] <- attr(x, "na.action")
    }
	x <- xxx
    
    # Free memory
	rm(xxx, xx, mm0)

    # Recheck data dimensions
    n <- dim(x)[1]
    m <- dim(x)[2]
	if (is.null(n)) 
        stop("invalid data size")
    # if (is.na(n) || n > 65536L) stop('size cannot be NA nor exceed
    # 65536')
    if (is.na(n)) 
        stop("size cannot be NA")
    if (n < 2) 
        stop("must have n \u2265 2 objects to cluster")

    # Update variance for multi-variate clustering
    if (nvars > 1) {
        if (verbose) write("---> Updating variance for multi-variate clustering...", "")
        v <- rowSums(x^2, na.rm = TRUE) / (m - 1)
    }

    # Reconstruct data from PCs if requested
    if (!is.null(nPC)) {
        if (verbose) write("---> Reconstructing data from PCs...", "")
        if (nPC >= 1 && nPC <= min(n, m)) {
            xSVD <- La.svd(t(x), nPC, nPC)
            eigVal <- xSVD$d
            expVar <- eigVal^2/sum(eigVal^2) * 100
            accVar <- sapply(seq(1, length(expVar)), function(r) sum(expVar[1:r]))
            x1 <- xSVD$u %*% diag(xSVD$d[1:nPC], nPC, nPC) %*% xSVD$vt
            x1 <- t(x1) - colMeans(x1)
            
        	# Cleanup memory from unnecessary variables
        	rm(xSVD)
        } else {
            stop(paste("invalid number of PCs,", 1, "\u2264 nPC \u2264", min(m, 
                n)))
        }
    }

	if (verbose) write("Agglomerative Hierarchical Clustering...", "")

    # Correlation matrix (fast calculation using BLAS library)
	if (verbose) write("---> Computing correlation/dissimilarity matrix...", "")
    if (!is.null(nPC)) {
        v1 <- rowSums(x1^2, na.rm = TRUE)
        # Correlation matrix (fast calculation using BLAS library)
        #r1 <- tcrossprod(x1/sqrt(v1))
        r1 <- fastCor(t(x1), nSplit = nSplit, upperTri = TRUE, verbose = verbose)
        # This is equivalent to upper triangular part of dissimilarity matrix
        r1 <- r1[col(r1) < row(r1)]
        # Variance of each variable (object/station)
        v1 <- v1/(m - 1)
    } else {
	    r <- fastCor(t(x), nSplit = nSplit, upperTri = TRUE, verbose = verbose)
	    # This is equivalent to upper triangular part of dissimilarity matrix
	    r <- r[col(r) < row(r)]

        x1 <- x
        v1 <- v
        r1 <- r
    }
    
    # Compute validation indices based on raw (100% of the total variance)
    # or PCA-filtered data?  Note that in both cases detrending and/or
    # standarding options are applied (before PCA)
    if (!rawStats) {
        x <- x1
        v <- v1
        r <- r1
    }
    
    # Dissimilarity matrix (correlation distance)
    d <- 1 - r1
    
    # Check dissimilarity matrix
    len <- as.integer(n * (n - 1)/2)
    if (length(d) != len) 
        (if (length(d) < len) 
            stop else warning)("data of improper length")
    
    
    # Cleanup memory from unnecessary variables
    rm(x1, r1)
    
    # Clustering method
    METHODS <- c("regional", "ward", "single", "complete", "average", "mcquitty", 
        "median", "centroid")
    method <- pmatch(method, METHODS) - 1
    if (is.na(method)) 
        stop("invalid clustering method")
    if (method == -1) 
        stop("ambiguous clustering method")
    
    # Check for restart clustering
    if (is.null(members)) 
        members <- rep(1, n) else if (length(members) != n) 
        stop("invalid length of members")
    
    if (verbose) write("---> Starting clustering process...", "")
    # Call Fortran subroutine for agglomerative hierarchical clustering
    storage.mode(d) <- "double"
    hcl <- .Fortran("HiClimR", n = n, len = len, method = as.integer(method), 
        ia = integer(n), ib = integer(n), crit = double(n), members = as.double(members), 
        var = v, diss = d, PACKAGE = "HiClimR")
    
    # interpret the output from previous step (such as merge, height, and
    # order lists)
    hcass <- .Fortran("hcass2", n = n, ia = hcl$ia, ib = hcl$ib, order = integer(n), 
        iia = integer(n), iib = integer(n), PACKAGE = "HiClimR")
    
    if (verbose) write("---> Constructing dendrogram tree...", "")
    # Construct 'hclust'/'HiClimR' dendogram tree
    tree <- list(merge = cbind(hcass$iia[1L:(n - 1)], hcass$iib[1L:(n - 
        1)]), height = hcl$crit[1L:(n - 1)], order = hcass$order, labels = rownames(x), 
        method = METHODS[method + 1], call = match.call(), dist.method = "correlation")
    class(tree) <- "hclust"
    
    tree$skip <- c(lonStep, latStep)
    names(tree$skip) <- c("lonStep", "latStep")
    
    if (!is.null(nPC)) {
        tree$PCA = cbind(eigVal = eigVal, expVar = expVar, accVar = accVar)
    }
    
    # return coordinates
    if (!is.null(lon) && !is.null(lat)) {
        tree$coords <- cbind(lon, lat)
        colnames(tree$coords) <- c("Lon", "Lat")
    }
    
    # return preprocessed raw or PCA-reconstructed data
    if (rawStats) {
        tree$data <- x
    } else {
        tree$data <- x1
    }
    
    # return number of variables and equivqlent number of columns
    tree$nvars <- nvars
    tree$ncols <- mm
    
    # Return mask vector
    if (length(mask) > 0) {
        tree$mask <- mask
    }
    
    # Return locations of missing values
    #if (length(attr(x, "na.action")) > 0) {
        tree$missVal <- missVal
    #}
    
    if (hybrid) {
        if (verbose) write("Hybrid Hierarchical Clustering...", "")

        if (verbose && tree$method == "regional") {
            write("---> WARNING: hybrid option is redundant when using regional linkage method!", 
                "")
        } else {
            if (is.null(kH)) {
                kH <- length(tree$height) - min(which(diff(tree$height) > 
                  mean(diff(tree$height)))) + 1
            }
            kH <- as.integer(kH)
            if (kH < 2) 
                stop("must have kH \u2265 2 objects to cluster")
            
            lenH <- as.integer(kH * (kH - 1)/2)
            
            methodH <- pmatch("regional", METHODS) - 1
            
            # Update variances dissimilarities of the upper part of the tree
		    if (verbose) write("---> Updating correlation/dissimilarity matrix...", "")
            cutTreeH <- cutree(tree, k = kH)
            RMH <- t(apply(tree$data, 2, function(r) tapply(r, cutTreeH, 
                mean)))
            vH <- apply(RMH, 2, var)
            rH <- fastCor(RMH, nSplit = nSplit, upperTri = TRUE, verbose = verbose)
            rH <- rH[col(rH) < row(rH)]
            dH <- 1 - rH
            
            # Call Fortran subroutine for agglomerative hierarchical clustering
            if (verbose) write("---> Reonstructing the upper part of the tree...", "")
            storage.mode(d) <- "double"
            hclH <- .Fortran("HiClimR", n = kH, len = lenH, method = as.integer(methodH), 
                ia = integer(kH), ib = integer(kH), crit = double(kH), 
                members = as.double(table(cutTreeH)), var = vH, diss = dH, 
                PACKAGE = "HiClimR")
            
            # interpret the output from previous step (such as merge, height, and
            # order lists)
            hcassH <- .Fortran("hcass2", n = kH, ia = hclH$ia, ib = hclH$ib, 
                order = integer(kH), iia = integer(kH), iib = integer(kH), 
                PACKAGE = "HiClimR")
            
            # Construct 'hclust'/'HiClimR' dendogram tree for the upper part
            treeH <- list(merge = cbind(hcassH$iia[1L:(kH - 1)], hcassH$iib[1L:(kH - 
                1)]), height = hclH$crit[1L:(kH - 1)], order = hcassH$order, 
                labels = names(table(cutree(tree, kH))), method = METHODS[methodH + 
                  1], call = match.call(), dist.method = "correlation")
            class(treeH) <- "hclust"
            
            # Add the new upper part of the tree to the original tree
            if (verbose) write("---> Merging upper and lower parts of the tree...", "")
            tree$treeH <- treeH
        }
    }
    
    # Plot dendrogram tree
    if (plot && dendrogram) {
        dev.new()
        if (hybrid && tree$method != "regional") {
            opar <- par(mfrow = c(1, 2))
            plot(tree, hang = hang, labels = labels, main = paste(METHODS[method + 
                1], "Method | Original Tree"))
            # rect.hclust(tree, k=kH, border='blue', cluster=cutTreeH)
            plot(treeH, hang = hang, labels = FALSE, main = paste("Regional Linkage |", 
                kH, "clusters"), axes = FALSE, ylab = "Maximum Inter-Regional Correlation")
            axis(2, at = pretty(treeH$height), labels = 1 - pretty(treeH$height))
            par(opar)
        } else {
            if (tree$method == "regional") {
                plot(tree, hang = hang, labels = labels, axes = FALSE, 
                  ylab = "Maximum Inter-Regional Correlation")
                axis(2, at = pretty(tree$height), labels = 1 - pretty(tree$height))
            } else {
                plot(tree, hang = hang, labels = labels)
            }
        }
    }
    
    # Cluster validation
    if (validClimR) {
        if (verbose) write("Calling cluster validation...", "")
        z <- validClimR(y = tree, k = k, minSize = minSize, alpha = alpha, verbose = verbose, 
            plot = plot, colPalette = colPalette, pch = pch, cex = cex)
        
        tree <- c(tree, z)
    }

    class(tree) <- c("hclust", "HiClimR")
    if (verbose) write("\nPROCESSING COMPLETED", "")
    
    # Print running time
    write("\nRunning Time:", "")
    print(proc.time()- start.time.proc) ; print(Sys.time()-start.time.sys)
    
    gc()
    
    # Output Tree
    tree
}
