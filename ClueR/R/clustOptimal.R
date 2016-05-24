#' Generate optimal clustering
#' 
#' Takes a clue output and generate the optimal clustering of the time-course data.
#' @param clueObj the output from runClue.
#' @param rep number of times the clustering is to be applied to find the best clustering result.
#' @param user.maxK user defined optimal k value for generating optimal clustering. If not provided, the optimal k that is identified by clue will be used.
#' @param visualize a boolean parameter indicating whether to visualize the clustering results.
#' @param ... pass additional parameter for controlling the plot if visualize is TRUE.
#' @return return a list containing optimal clustering object and enriched kinases or gene sets. 
#'
#' @export
#' @examples
#' # simulate a time-series data with 4 distinctive profile groups and each group with
#' # a size of 50 phosphorylation sites.
#' simuData <- temporalSimu(seed=1, groupSize=50, sdd=1, numGroups=4)
#' 
#' # create an artificial annotation database. Generate 20 kinase-substrate groups each
#' # comprising 10 substrates assigned to a kinase.
#' # among them, create 4 groups each contains phosphorylation sites defined to have the
#' # same temporal profile.
#' kinaseAnno <- list()
#' groupSize <- 50
#' for (i in 1:4) {
#'  kinaseAnno[[i]] <- paste("p", (groupSize*(i-1)+1):(groupSize*(i-1)+10), sep="_")
#' }
#'
#' for (i in 5:20) {
#'  set.seed(i)
#'  kinaseAnno[[i]] <- paste("p", sample.int(nrow(simuData), size = 10), sep="_")
#' }
#' names(kinaseAnno) <- paste("KS", 1:20, sep="_")
#'
#' # run CLUE with a repeat of 2 times and a range from 2 to 7
#' set.seed(1)
#' clueObj <- runClue(Tc=simuData, annotation=kinaseAnno, rep=5, kRange=7)
#' 
#' # visualize the evaluation outcome
#' Ms <- apply(clueObj$evlMat, 2, mean, na.rm=TRUE)
#' Ss <- apply(clueObj$evlMat, 2, sd, na.rm=TRUE)
#' library(Hmisc)
#' errbar(1:length(Ms), Ms, Ms+Ss, Ms-Ss, cex=1.2, type="b", xaxt="n", xlab="k", ylab="E")
#' axis(1, at=1:6, labels=paste("k=", 2:7, sep=""))
#' 
#' # generate optimal clustering results using the optimal k determined by CLUE
#' best <- clustOptimal(clueObj, rep=3, mfrow=c(2, 3))
#' 
#' # list enriched clusters
#' best$enrichList
#' 
#' # obtain the optimal clustering object (not run)
#' # best$clustObj
#' 
clustOptimal <- function(clueObj, rep, user.maxK=NULL, visualize=TRUE, ...) {
  
  bstPvalue <- 1
  bst.clustObj <- c()
  bst.evaluation <- c()
  for(i in 1:rep) {
    clustObj <- c()
    if (clueObj$clustAlg == "cmeans") {
      if (is.null(user.maxK)) {
        # use clue determined max k value
        clustObj <- cmeans(clueObj$Tc, centers=clueObj$maxK, iter.max=50, m=1.25)
      } else {
        # use user specific k value
        clustObj <- cmeans(clueObj$Tc, centers=user.maxK, iter.max=50, m=1.25)
      }
    } else {
      if (is.null(user.maxK)) {
        # use clue determined max k value
        clustObj <- kmeans(clueObj$Tc, centers=clueObj$maxK, iter.max=50)
        # set the membership as 1 for all partitions
        clustObj$membership <- matrix(1, nrow=nrow(clueObj$Tc), ncol=nrow(clustObj$centers))
        rownames(clustObj$membership) <- names(clustObj$cluster)
        
      } else {
        # use user specific k value
        clustObj <- kmeans(clueObj$Tc, centers=user.maxK, iter.max=50)
        # set the membership as 1 for all partitions
        clustObj$membership <- matrix(1, nrow=nrow(clueObj$Tc), ncol=nrow(clustObj$centers))
        rownames(clustObj$membership) <- names(clustObj$cluster)
      }
    }
    
    evaluation <- clustEnrichment(clustObj, clueObj$annotation, clueObj$effectiveSize, clueObj$pvalueCutoff)
    currentPvalue <- evaluation$fisher.pvalue
    
    if (currentPvalue < bstPvalue) {
      bstPvalue <- currentPvalue
      bst.clustObj <- clustObj
      bst.enrichList <- evaluation$enrich.list
    }
  }
  
  if(visualize) {
    fuzzPlot(clueObj$Tc, clustObj = bst.clustObj, ...)
  }
  
  results <- list()
  results$clustObj <- bst.clustObj
  results$enrichList <- bst.enrichList
  return(results)
}
