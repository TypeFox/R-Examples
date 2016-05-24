
#
#' Split markers into chromosomes
#' 
#' This function splits markers into linkage groups (LG), which ideally represent chromosomes.
#' The split is based on hierarchical clustering with a single linkage distance.
#'  
#' @param rf Matrix of pairwise recombination frequencies.
#' @param height Threshold value for grouping the markers.
#' @param nchr Expected number of chromosomes.
#' @param method Default is "single", which is used for the hierarchical clustering.
#' @param filter Logical, if the result should be filtered or not. Default is FALSE. Creates zeros for the markers below the threshold.
#' @param thresh Threshold for filtering. Default is 0.05, i.e. linkage groups with less than 5\% of markers, are filtered out.
#' @param rm.dup Logical, if the duplicated markers should be filtered out.
#' TRUE is highly recommended because the markers have no added value for the linkage map.
#' @return Vector of cluster relationship. Same length and order as the matrix of recombination frequencies.
#' @examples
#' data(simTetra)
#' simTetrageno<-bases2genotypes(simTetra, 4)
#' rfMat<-calcRec(simTetrageno, 4)
#' splitChr(rfMat, nchr = 7)
#' @export
splitChr <- function(rf, height = 0.4, nchr = NULL, method = "single", filter = FALSE, 
                     thresh = 0.05, rm.dup = TRUE){
  split <- data.frame(names = rownames(rf), split = 1, dup = 0)
  rownames(split) <- split$names
  if(rm.dup == TRUE){
    zeroes <- which(rf == 0, arr.ind = TRUE)
    offdiag <- which(zeroes[, 1] < zeroes[, 2])
    for(j in offdiag){
      split$split[zeroes[j, 2]] <- 0
      split$dup[zeroes[j, 2]] <- zeroes[j, 1]
    }
  } 
  rfsub <- rf[split$split > 0, split$split > 0]
  tree <- stats::hclust(as.dist(rfsub), method = method)
  minleaves <- thresh * ncol(rfsub)
  if(!is.null(nchr)){
    output <- cutree(tree = tree, k = nchr)
    while(filter){
      filtClus <- which(table(output) < minleaves)
      if(length(filtClus) > 0){
        tooFilt <- which(output %in% filtClus)
        if(length(tooFilt) < ncol(rfsub)){
          split$split[split$split > 0][tooFilt] <- 0
          rfsub <- rfsub[-tooFilt, -tooFilt]
          tree <- stats::hclust(as.dist(rfsub), method = method)
          output <- stats::cutree(tree = tree, k = nchr)
        }else{
          stop("Could not split data into ", nchr, " clusters." )
        }
      }else{
        filter <- FALSE
      }
    }
  }else{
    output <- stats::cutree(tree = tree, h = height)  
    if(filter){
      filtClust <- which(table(output) < minleaves)
      output[output %in% filtClust] <- 0
    }
  }
  split$split[split$split > 0] <- output
  return(split)   
}
