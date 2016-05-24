#' Run CLUster Evaluation
#' 
#' Takes in a time-course matrix and test for enrichment of the clustering using cmeans or kmeans clustering algorithm with a reference annotation.
#' @param Tc a numeric matrix to be clustered. The columns correspond to the time-course and the rows correspond to phosphorylation sites.
#' @param annotation a list with names correspond to kinases and elements correspond to substrates belong to each kinase.
#' @param rep number of times the clustering is to be applied. This is to account for variability in the clustering algorithm.
#' @param kRange the range of k to be tested for clustering.
#' @param clustAlg the clustering algorithm to be used. The default is cmeans clustering.
#' @param effectiveSize the size of annotation groups to be considered for calculating enrichment. Groups that are too small
#' or too large will be removed from calculating overall enrichment of the clustering.
#' @param pvalueCutoff a pvalue cutoff for determining which kinase-substrate groups to be included in calculating overall enrichment of the clustering.
#' @param alpha a penalty factor for penalizing large number of clusters.
#' @return a clue output that contains the input parameters used for evaluation and the evaluation results. Use ls(x) to see details of output. 'x' be the output here.
#' @export
#' @examples
#' # load the human ES phosphoprotoemics data (Rigbolt et al. Sci Signal. 4(164):rs3, 2011)
#' data(hES)
#' # load the PhosphoSitePlus annotations (Hornbeck et al. Nucleic Acids Res. 40:D261-70, 2012)
#' data(PhosphoSite)
#' 
#' # make a subset of hES dataset for demonstrating the example in a short time frame
#' ids <- c("CK2A1", "ERK1", "ERK2", "CDK7", 
#' "p90RSK", "p70S6K", "PKACA", "CDK1", "DNAPK", "ATM", "CDK2")
#' hESs <- hES[rownames(hES) %in% unlist(PhosphoSite.human[ids]),]
#' 
#' # run CLUE with a repeat of 3 times and a range from 2 to 13
#' set.seed(2)
#' clueObj <- runClue(Tc=hESs, annotation=PhosphoSite.human, rep=2, kRange=13)
#' 
#' # visualize the evaluation outcome
#' Ms <- apply(clueObj$evlMat, 2, mean, na.rm=TRUE)
#' Ss <- apply(clueObj$evlMat, 2, sd, na.rm=TRUE)
#' library(Hmisc)
#' errbar(1:length(Ms), Ms, Ms+Ss, Ms-Ss, cex=1.2, type="b", xaxt="n", xlab="k", ylab="E")
#' axis(1, at=1:12, labels=paste("k=", 2:13, sep=""))
#' 
#' # generate the optimal clustering results
#' best <- clustOptimal(clueObj, rep=10, mfrow=c(3, 4))
#' 
#' # list enriched clusters
#' best$enrichList
#' 
#' # obtain the optimal clustering object (not run)
#' # best$clustObj
#' 
runClue <- function(Tc, annotation, rep=10, kRange, clustAlg="cmeans", effectiveSize=c(5, 100), pvalueCutoff=0.05, alpha=0.5) {
  
  # standardize the matrix by row
  means <- apply(Tc, 1, mean)
  stds <- apply(Tc, 1, sd)
  tmp <- sweep(Tc, 1, means, FUN="-")
  Tc <- sweep(tmp, 1, stds, FUN="/")
  
  ## filter the annotation groups that has no entry from the Tc
  annotation.intersect <- lapply(annotation, intersect, rownames(Tc))
  annotation.filtered <- annotation.intersect[lapply(annotation.intersect, length) > 0]
  
  # apply CLUE
  repeat.list <- list()
  for(rp in 1:rep) {
    cat("repeat", rp, "\n");
    enrichment <- c()
    for (k in 2:kRange) {
      clustered <- c()
      if (clustAlg == "cmeans") {
        clustered <- cmeans(Tc, centers=k, iter.max=50, m=1.25)
      } else if (clustAlg == "kmeans"){
        clustered <- kmeans(Tc, centers=k, iter.max=50)
      } else {
        print("Unknown clustering algorithm specified. Using cmeans clustering instead")
        clustered <- cmeans(Tc, centers=k, iter.max=50, m=1.25)
      }
      
      # compute clustering enrichment
      evaluate <- clustEnrichment(clustered, annotation.filtered, effectiveSize, pvalueCutoff)
      fisher.pvalue <- evaluate$fisher.pvalue
      escore <- -log10(fisher.pvalue) - alpha * nrow(clustered$centers)
      
      enrichment <- c(enrichment, escore)
    }
    repeat.list[[rp]] <- enrichment
  }
  
  # combine the multiple testing results
  x <- do.call(rbind, repeat.list)
  # transform the pvalue 
  #x.transform <- -log10(x)
  # scale the values into [0, 1]
  x.normalize <- (x - min(x)) / (max(x) - min(x))
  rownames(x.normalize) <- paste("repeat", 1:rep, sep="")
  colnames(x.normalize) <- paste("k", 2:kRange, sep="=")
  
  # identify the k that maximize the enrichment
  maxK <- which.max(apply(x.normalize, 2, mean)) + 1
  
  # return the evaluation results
  result <- list()
  # input parameters
  result$Tc <- Tc
  result$annotation <- annotation.filtered
  result$clustAlg <- clustAlg
  result$effectiveSize <- effectiveSize
  result$pvalueCutoff <- pvalueCutoff
  # clue parameters
  result$evlMat <- x.normalize
  result$maxK <- maxK
  
  return(result)
}
