subsetDataForHierTrees <- function(oneFileHier, HierID) {
  oneFileHier <- matrix(oneFileHier[1:which(oneFileHier[,6] == HierID),], ncol=7)  #stop at HierID to avoid taking children
  oneFileHier  <- matrix(oneFileHier[which(!duplicated(oneFileHier[,2])),], ncol=7) #delete repeats
  oneFileHier  <- matrix(oneFileHier[!is.na(oneFileHier[,2]),], ncol=7) #delete NAs; ncol is hard coded to be 7 columns so that R doesn't convert to a vector when there is a single row.  
  if(any(oneFileHier[,2] == "unranked clade"))
    oneFileHier  <- oneFileHier[-which(oneFileHier[,2] == "unranked clade"),] #delete unranked
  return(oneFileHier)
}

MergeTaxonomies <- function(i, j) {
  combined <- i
  `%ni%` = Negate(`%in%`) 
  outlier.pos <- which(j %ni% i)
  for (outlier.index in sequence(length(outlier.pos))) {
    if (outlier.pos[outlier.index] == 1) {
      paste("Haven't dealt with this yet")
      return(combined)
    }
    previous.pos <- which(grepl(j[outlier.pos[outlier.index] - 1], combined))[1]
    if ((previous.pos+1) <= length(combined))
      combined <- c(combined[1:previous.pos], j[outlier.pos[outlier.index]], combined[(previous.pos+1): length(combined)]) 
    else
      combined <- c(combined[1:previous.pos], j[outlier.pos[outlier.index]])
  }
  return(combined)
}

CombineHierarchyInfo <- function(MyHiers) {
  CombFiles <- matrix(nrow=0, ncol=7)
  longestHierTaxon <- 0  #start at 0, so it accepts the first file as the longest
  MergedTax <- NULL
  for(i in sequence(length(MyHiers))) {
    OFH <- OneFileHierarchy(MyHiers[i])
    if(dim(OFH)[1] <= 1)
      cat(paste(OFH[1,1], "has no higher taxonomic rankings, so it must be dropped\n"))
    if(dim(OFH)[1] > 1){
      oneFile <- subsetDataForHierTrees(OFH, GetHierID(MyHiers[i]))
      if(dim(oneFile)[1] <= 1)
        cat(paste(OFH[1,1], "has no higher taxonomic rankings, so it must be dropped\n"))
      if(dim(oneFile)[1] > 1) {
        Tax <- oneFile[,2]
        MergedTax <- MergeTaxonomies(Tax, MergedTax)
        if(length(oneFile[,2]) > longestHierTaxon) {
          longestHierTaxon <- max(longestHierTaxon, length(oneFile[,2]))    
          CombFiles <- rbind(oneFile, CombFiles)  #puts longest hierarchies first in combined files
          CombFiles <- as.data.frame(CombFiles, stringsAsFactors=FALSE)
        }
      }
      else
        CombFiles <- rbind(CombFiles, oneFile)  #puts shorter hierarchies after
    }
  }
  return(list(CombFiles, MergedTax))
}

MakeTreeData <- function(MyHiers) {
  MyHiers <- RemoveNAFiles(MyHiers)
  CombFiles <- CombineHierarchyInfo(MyHiers)
  whichColumns <- CombFiles[[2]]
  TreeData <- data.frame(matrix(nrow=length(MyHiers), ncol=length(whichColumns)))
  colnames(TreeData) <- whichColumns
  for(i in sequence(length(MyHiers))) {
    #here go one at a time and add each row
    oneFile <- subsetDataForHierTrees(OneFileHierarchy(MyHiers[i]), GetHierID(MyHiers[i]))
    for(j in sequence(dim(oneFile)[1])){
      colPosition <- which(colnames(TreeData) == oneFile[j,2])
      TreeData[i,colPosition] <- oneFile[j,1]
      TreeData <- as.data.frame(TreeData, stringsAsFactors=T)
    }
  }
  return(TreeData)
}  

RepeatDataToDrop <- function(TreeData) {
  #Drop any columns that are the same or NA
  if(length(unique(TreeData)) == 1)
    return(TRUE)
  if(any(is.na(TreeData)))
    return(TRUE)
  else
    return(FALSE)
}

DropADim <- function(TreeData) {
  while(any(is.na(TreeData))){
    rowPercent <- apply(is.na(TreeData), 1, sum)/dim(TreeData)[2]
    colPercent <- apply(is.na(TreeData), 2, sum)/dim(TreeData)[1]
    maxPercent <- max(c(rowPercent, colPercent))
    if(maxPercent %in% rowPercent)
      TreeData  <- TreeData[-which(rowPercent == maxPercent), ]
    if(maxPercent %in% colPercent)
      TreeData  <- TreeData[,-which(colPercent == maxPercent)]
  }
  return(TreeData)
}

AutofillTaxonNames <- function(TreeData){
  #autofill in missing data with child taxon names
  for(i in sequence(dim(TreeData)[1])){
    columnNAs <- which(is.na(TreeData[i,]))
    columnInfo <- which(!is.na(TreeData[i,]))
    if(any(columnNAs)){
      if(max(columnNAs) > max(columnInfo)){
        columnNAs <- columnNAs[-which(columnNAs > max(columnInfo))]
      }
      for (j in rev(sequence(length(columnNAs)))){
        ChildTaxonPlace <- which(columnInfo > columnNAs[j])[1]
        #TreeData[i, columnNAs[j]] <- paste(colnames(TreeData)[columnNAs[j]], TreeData[i, columnInfo[ChildTaxonPlace]], sep="")
        TreeData[i, columnNAs[j]] <- paste(TreeData[i, columnInfo[ChildTaxonPlace]], sep="")
      }
    }
  }
  return(TreeData)
}
  

MakeHierarchyTree <- function(MyHiers, missingData=NULL, includeNodeLabels=TRUE, userRanks=NULL) {
  TreeData <- MakeTreeData(MyHiers)
  pattern <- paste("~", paste(colnames(TreeData), sep="", collapse="/"), sep="")
  if(!is.null(userRanks)){
    TreeData <- TreeData[,which(colnames(TreeData) %in% userRanks)]
    TreeData <- AutofillTaxonNames(TreeData)   
    #TreeData <- DropADim(TreeData)
    pattern <- paste("~", paste(colnames(TreeData), sep="", collapse="/"), sep="")   
  }
  if(any(apply(TreeData, 2, RepeatDataToDrop))) {
    if(any(is.na(TreeData))){
      TreeData <- AutofillTaxonNames(TreeData)
      if(colnames(TreeData[dim(TreeData)[2]]) != strsplit(pattern, "/")[[1]][length(strsplit(pattern, "/")[[1]])]){
        if(is.null(missingData))
          stop("Tip taxa are differing hierarchical ranks, you need to choose an option in missingData whether to drop taxa or ranks.")
        if(missingData == "pruneTaxa")
          TreeData <- TreeData[-which(is.na(TreeData[,dim(TreeData)[2]])),]
      }
    }
    DataToDrop <- which(apply(TreeData, 2, RepeatDataToDrop))
    if(any(DataToDrop))
      pattern <- paste("~", paste(colnames(TreeData)[-DataToDrop], sep="", collapse="/"), sep="")
  }
  if(colnames(TreeData[dim(TreeData)[2]]) != strsplit(pattern, "/")[[1]][length(strsplit(pattern, "/")[[1]])]){
    paste("Your hierarchy files contain information to the", colnames(TreeData[dim(TreeData)[2]]), "level, however not all taxa have this information. In order to make a tree the tips must align, so information was pruned to", strsplit(pattern, "/")[[1]][length(strsplit(pattern, "/")[[1]])])
    paste("If you want to make a", colnames(TreeData[dim(TreeData)[2]]), "tree, change the missingData argument to pruneTaxa, which will remove taxa without", colnames(TreeData[dim(TreeData)[2]]))
  }
  if(pattern == "~")
    stop("Error in Tree Building: try MakeTreeData(MyHiers) to see if there is hierarchical data associated with your files")
  fo <- as.formula(pattern)
  TreeData <- as.data.frame(apply(TreeData, 2, factor))
  tree <- ladderize(as.phylo.formula(fo, data=TreeData))  
  if(includeNodeLabels)
    tree <- makeNodeLabel(tree, method="u", nodeList=NodeLabelList(MyHiers, "all", missingData=missingData), fixed=TRUE)  #maybe change this later when other options
  return(tree)
}

ReturnTaxSet <- function(Taxon, TreeData) {
	whichRows <- which(TreeData == Taxon, TRUE)[,1]
	return(TreeData[whichRows, dim(TreeData)[2]])
}

NodeLabelList <- function(MyHiers, label="all", missingData) {  #also make an option to just label genus, etc. 
  TreeData <- MakeTreeData(MyHiers)
  if(any(is.na(TreeData))){
    TreeData <- AutofillTaxonNames(TreeData)
    if(any(is.na(TreeData)) && !is.null(missingData)){
      if(missingData == "pruneTaxa")
        TreeData <- TreeData[-which(is.na(TreeData[,dim(TreeData)[2]])),]
    }
  }
  DataToDrop <- which(apply(TreeData, 2, RepeatDataToDrop))
  if(any(DataToDrop))
    TreeData <- TreeData[,-DataToDrop]
  if(label == "all")
    uniqueTaxSets <- c(apply(TreeData[,-dim(TreeData)[2]], 2, unique), recursive=T)
  ListOfSpeciesPerNode <- lapply(uniqueTaxSets, ReturnTaxSet, TreeData= TreeData)
  names(ListOfSpeciesPerNode) <- uniqueTaxSets
  ListOfSpeciesPerNode <- ListOfSpeciesPerNode[which(lapply(ListOfSpeciesPerNode, length) > 1)]
  if(any(c(lapply(ListOfSpeciesPerNode, duplicated), recursive=T))){
    for(i in rev(sequence(length(ListOfSpeciesPerNode)))){
      if(any(duplicated(ListOfSpeciesPerNode[[i]])))
        ListOfSpeciesPerNode[[i]] <- ListOfSpeciesPerNode[[i]][-which(duplicated(ListOfSpeciesPerNode[[i]]))]
    }
  }
  if(any(c(lapply(ListOfSpeciesPerNode, length), recursive=T) == 1))
    ListOfSpeciesPerNode <- ListOfSpeciesPerNode[-which(c(lapply(ListOfSpeciesPerNode, length), recursive=T) == 1)]
  if(any(is.na(names(ListOfSpeciesPerNode))) || any(names(ListOfSpeciesPerNode) == "Not assigned")){
    whichMissing <- c(which(is.na(names(ListOfSpeciesPerNode))), which(names(ListOfSpeciesPerNode) == "Not assigned"))
    ListOfSpeciesPerNode <- ListOfSpeciesPerNode[-whichMissing]
  }  
  return(ListOfSpeciesPerNode)
}

