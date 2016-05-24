GatherSynonyms <- function(MyHiers, output=c("detail", "counts")) {
  MyHiers <- RemoveNAFiles(MyHiers)
  output <- match.arg(output)
  syns <- matrix(ncol=3, nrow=0)
  colnames(syns) <- c("Taxon", "hierID", "Synonym")
  SynCounts <- matrix(nrow=length(MyHiers), ncol=3)
  colnames(SynCounts) <- c("Taxon", "hierID", "NumberOfSynonyms")
  for(i in sequence(length(MyHiers))){
    resOneFile <- OneFileHierarchy(MyHiers[i])
    Taxon <- resOneFile[which(resOneFile[,6] == GetHierID(MyHiers[i])), 1]
    whichSyn <- c(which(resOneFile[,3] == "Synonym"), which(resOneFile[,3] == "synonym"))
    SynCounts[i,] <- c(Taxon, GetHierID(MyHiers[i]), length(whichSyn))
    if(length(whichSyn) > 0) {
      for(j in sequence(length(whichSyn))) {
        syns <- data.frame(rbind(syns, c(Taxon, GetHierID(MyHiers[i]), resOneFile[whichSyn[j], 1])), stringsAsFactors=FALSE) 
      }
    }
  }
  SynCounts <- data.frame(SynCounts, stringsAsFactors=FALSE)
  SynCounts[,3] <- as.numeric(SynCounts[,3])
  if(output == "detail")
    return(syns)
  if(output == "counts")
    return(SynCounts)
}
