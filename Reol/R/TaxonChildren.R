TaxonChildren <- function (MyHiers){
  MyHiers <- RemoveNAFiles(MyHiers)
  TC <- matrix(ncol=5, nrow=0)
  colnames(TC) <- c("Taxon", "TaxonChild", "TaxonRank", "eolID", "hierID")  
  for(tax in sequence(length(MyHiers))){
    oneFile <- OneFileHierarchy(MyHiers[tax])
    Taxon <- oneFile[which(oneFile[,6] == GetHierID(MyHiers[tax])), 1]
    whichRows <- which(oneFile[,5] == GetHierID(MyHiers[tax]))
    if(length(whichRows) > 0){
      oneFile <- matrix(oneFile[whichRows,], nrow=length(whichRows))  #find which are to parent
      #if(any(is.na(oneFile[,2])))  
      #  oneFile <- oneFile[-which(is.na(oneFile[,2])),] #remove any NAs
      if(dim(oneFile)[1] > 0){
        for(i in sequence(dim(oneFile)[1])) {
          TC <- rbind(TC, c(Taxon, oneFile[i,c(1,2,4,6)]))
        }
      }
    }
  }  
  TC <- data.frame(TC, stringsAsFactors=FALSE)
  return(TC)
}


TaxonParents <- function(MyHier) {
  hierID <- GetHierID(MyHier)
  parents <- subsetDataForHierTrees(OneFileHierarchy(MyHier), hierID)
  parents <- data.frame(paste(parents[,2], parents[,1]))
  names(parents) <- ""
  return(parents)
}