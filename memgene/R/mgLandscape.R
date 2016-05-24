
mgLandscape <- function(resistance, genD, coords, euclid=TRUE, forwardPerm=100, forwardAlpha=0.05, finalPerm=1000, verbose=TRUE) {
 
   if (!(class(resistance) %in% c("RasterLayer", "RasterStack", "RasterBrick"))) {
    stop("memgene: resistance must be a resistance surface provided in a RasterLayer format, or if multiple surfaces are being compared use a RasterStack or RasterBrick")
  }
  
  if (grepl("longlat", projection(resistance))) {
        stop("memgene: raster must be in a planar projection that is the same as coordinates.  Use raster::projectRaster()", call.=FALSE)
  }
  
    
  if (sum(sum(resistance <= 0)[], na.rm=TRUE) > 0) {
        stop("memgene: all values on resistance surface must be either missing (NA) or >=1", call.=FALSE)
  }
  
  if (euclid) {
    results <- vector("list", nlayers(resistance)+1)
  }
  else {
    results <- vector("list", nlayers(resistance))
  }
  
  for (iMap in 1:length(results)) {
    
    results[[iMap]] <- NA
    
    if (iMap==1 && euclid) {
      if (verbose) cat("\nAnalyzing Euclidean surface (landscape model ", iMap, " of ", length(results), ")\n", sep="")
      if (verbose) cat("\tExtracting Moran's eigenvectors from Euclidean distance matrix\n")
      mem <- tryCatch(mgMEM(dist(coords)), error=function(e) return(NA))
    }
    else {
      if (verbose) cat("\nAnalyzing resistance surface (landscape model ", iMap, " of ", length(results), ") ", "[", names(resistance)[[ifelse(euclid, iMap-1, iMap)]], "]\n", sep="")
      cat("\tCalculating least-cost path distance matrix\n")
      trans <- transition(resistance[[ifelse(euclid, iMap-1, iMap)]], function(x) 1/mean(x), 8)
      lcp <- as.matrix(costDistance(trans, SpatialPoints(coords)))
      
      if (verbose) cat("\tExtracting Moran's eigenvectors from least-cost path distance matrix\n")
      mem <- tryCatch(mgMEM(lcp), error=function(e) return(NA))
    }

    if (any(is.na(mem))) {
      if (verbose) cat("\tFailure to extract Moran's eigenvectors.  Skipping this surface.\n")
      if (iMap==1 && euclid) {
        results[[iMap]] <- data.frame(hypothesis="Euclidean", matrix(NA, 1, 9))
      }
      else {
        results[[iMap]] <- data.frame(hypothesis=names(resistance[[ifelse(euclid, iMap-1, iMap)]]), matrix(NA, 1, 9))
      }
    }
    else {
      pos <- mem$vectorsMEM[, mem$valuesMEM > 0]
      neg <- mem$vectorsMEM[, mem$valuesMEM < 0]
  
      if (verbose) cat("\tForward selections of positive Moran's eigenvectors\n")
      selectedPos <- tryCatch(na.omit(mgForward(genD, pos, perm=forwardPerm,
                      alpha=forwardAlpha)$selectedMEM),
                        error=function(e) return(0))
      if (verbose) cat("\t----Selected:", ifelse(length(selectedPos)==0, "None", paste(sort(selectedPos), collapse=", ")), "\n")
      
      if (verbose) cat("\tForward selections of negative Moran's eigenvectors\n")
      selectedNeg <- tryCatch(na.omit(mgForward(genD, neg, perm=forwardPerm,
                          alpha=forwardAlpha)$selectedMEM),
                              error=function(e) return(0))
      if (verbose) cat("\t----Selected:", ifelse(length(selectedNeg)==0, "None", paste(sort(selectedNeg), collapse=", ")), "\n")

      memSelected <- cbind(pos[, selectedPos], neg[, selectedNeg])

      if (ncol(memSelected)==0) {
          if (verbose) cat("\tThere are no significant Moran's eigenvectors.\n\tSpatial genetic signal is possibly too weak for memgene analysis.\n\tSkipping this surface.\n")
          results[[iMap]] <- data.frame(hypothesis=names(resistance[[ifelse(euclid, iMap-1, iMap)]]), matrix(NA, 1, 9))
      }
      else {
        if (verbose) cat("\tPartitioning spatial genetic variation\n")
        results[[iMap]] <- mgVarPart(genD, memSelected, coords, perm=finalPerm)
        if (iMap==1 && euclid) {
          results[[iMap]] <- data.frame(hypothesis="Euclidean", matrix(t(results[[iMap]]), 1, 14, byrow=FALSE), forwardPerm, finalPerm)
        }
        else {
          results[[iMap]] <- data.frame(hypothesis=names(resistance[[iMap-1]]), matrix(t(results[[iMap]]), 1, 14, byrow=FALSE), forwardPerm, finalPerm)
        }
        results[[iMap]] <- results[[iMap]][, -c(3,4,5,6,13,15)]
      }
    }
  }
  allResults <- list()
  allResults$summary <- do.call(rbind, results)
  names(allResults$summary) <- c("model", "[abc]", "P[abc]", "[a]", "P[a]", "[c]", "P[c]", "[b]", "[d]", "forwardPerm", "finalPerm")
  class(allResults) <- "mgLandscape"
  return(allResults)

}