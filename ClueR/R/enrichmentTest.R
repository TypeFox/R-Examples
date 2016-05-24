#' Fisher's exact test-based enrichment test
#' 
#' Takes a vector of names representing phosphorylation sites that are partitioned in the same cluster 
#' and an kinase-substrate annotation. Test for enrichment of the kinase based on the name vector. 
#' 
#' 
#' @param clust a vector of names representing phosphorylation sites that are partitioned in the same cluster
#' @param annotation a list with names correspond to kinases and elements correspond to substrates belong to each kinase
#' @param universe the universe of names to compare against
#' @param alter indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less"
#' @return a matrix that contains enrichment of each kinase based on the input name vector.
#' 
#' @export
#' 
enrichmentTest <- function(clust, annotation, universe, alter="greater") {
  
  fisherTest.mat <- matrix("NA", ncol=4, nrow=length(annotation))
  colnames(fisherTest.mat) <- c("kinase", "pvalue", "# of substrates", "substrates")
  for (i in 1:length(annotation)) {
    di <- length(intersect(clust, annotation[[i]]))
    dn <- length(setdiff(clust, annotation[[i]]))
    ndi <- length(setdiff(annotation[[i]], clust))
    ndn <- length(setdiff(universe, union(clust, annotation[[i]])))
    
    p <- fisher.test(rbind(c(di, ndi), c(dn, ndn)), alternative=alter)$p.value
    substrates <- paste(intersect(clust, annotation[[i]]), collapse="|")
    fisherTest.mat[i,] <- c(names(annotation)[i], p, di, substrates)
  }
  
  # sort by pvalues and return the result
  fisherTest.mat <- fisherTest.mat[order(as.numeric(fisherTest.mat[, 2])),]
  return(fisherTest.mat)
}
