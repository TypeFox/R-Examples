MatchHierPageToEOLdata <- function(MyHiers, EOLdata){
  MyHiers <- RemoveNAFiles(MyHiers)
  if(class(MyHiers) == "list")
    fileNames <- names(MyHiers)
  else
    fileNames <- MyHiers
  matchedData <- matrix(nrow=length(MyHiers), ncol=3)
  matchedData[,2] <- as.character(sapply(fileNames, GetHierID))
  for(i in sequence(length(MyHiers))) {
    resOneFile<-OneFileHierarchy(MyHiers[i])
    matchedData[i,1] <- resOneFile[which(matchedData[i,2]==resOneFile[,6]), 1] #taxonName
    matchedData[i,3] <- resOneFile[which(matchedData[i,2]==resOneFile[,6]), 4] #HierID
  }
  matches <- match(matchedData[,3], EOLdata[,2])
  matchedData <- cbind(matchedData, EOLdata[matches, 3:dim(EOLdata)[2]])
  colnames(matchedData) <- c( "HierTaxon", "HierID", "eolID", colnames(EOLdata[3:dim(EOLdata)[2]]))
  matchedData <- data.frame(matchedData, row.names=1, stringsAsFactors=FALSE)
  return(matchedData)
}


MatchDataToTreeTips <- function(Tree, Data){
  if(length(grep("HierID", colnames(Data), ignore.case=T)) == 0) 
    stop("EOLdata must be matched to HierarchyID first, use MatchHierPageToEOLdata(MyHiers, EOLdata)")
  if(any(is.na(Data)))
    Data <- Data[-unique(which(is.na(Data), arr.ind=TRUE)[,1]),]
  if(length(grep("Taxon", colnames(Data), ignore.case=T)) > 0)  #if taxon names are a column, then make them rownames
    rownames(Data) <- Data[,grep("Taxon", colnames(Data), ignore.case=T)]
  NewOrderMatchedData <- data.frame(matrix(ncol=dim(Data)[2], nrow=length(Tree$tip.label)))
  colnames(NewOrderMatchedData) <- colnames(Data)
  for(i in sequence(length(Tree$tip.label))){
    if(any(rownames(Data) %in% Tree$tip.label[i])) {
      NewOrderMatchedData[i,] <- Data[which(rownames(Data) %in% Tree$tip.label[i]),]
      rownames(NewOrderMatchedData)[i] <- rownames(Data[which(rownames(Data) %in% Tree$tip.label[i]),])
    }
    else {
      NewOrderMatchedData[i,] <- rep(NA, dim(NewOrderMatchedData)[2])
      rownames(NewOrderMatchedData)[i] <- Tree$tip.label[i]
    }
  }
  if(all(rownames(NewOrderMatchedData) == Tree$tip.label))
    return(NewOrderMatchedData)
  else
    stop("something wonky in data")
}









