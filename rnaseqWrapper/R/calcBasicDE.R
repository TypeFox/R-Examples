calcBasicDE <-
function(data,colID="_mean_FPKM_log2",
         whichIndex="allPairwise",matchEnd=TRUE,appendName="diff"){
  
  ## data should be matrix or data.frame
  ## colID tells which measure to use for DE analysis, can be a vector for several things
  ##   either character (to be grepped) or numeric to indicate which cols to analyze
  ##   accepts regexp arguments, so be careful
  ## whichIndex tells which groupID (e.g. "Control") the others should compare against.
  ##    Can be numeric (which of the colID to use) or characther for grep (can have regexp)
  ##    whichIndex="allPairwise" calculates all pairwise comparisons
  ## appendName adds a bit to the end of the name, to prevent it from matching in future things
  
  colIndex <- NULL # set a variable to hold the column index
  if(is.character(colID)){
    for(k in 1:length(colID)){
      if(matchEnd){
        colIndex <- c(colIndex,grep(paste(colID[k],"$",sep=""),colnames(data),value=TRUE))  
      } else {
        colIndex <- c(colIndex,grep(colID[k],colnames(data),value=TRUE))  
      }
      
    } 
  } else if(is.numeric(colID)) {
    colIndex <- colID
  } else {
    stop("ERROR: colID must be either numeric or character")
  }
  
  
  tempData <- data[,colIndex]
  diffsFrame <- data.frame(row.names=rownames(tempData))
  
  ## Set the index value
  if(whichIndex!="allPairwise"){
    if(whichIndex==0|whichIndex==""|is.null(whichIndex)){
      whichIndex <- "allPairwise"
    } else if(is.character(whichIndex)){
      whichIndex <- grep(whichIndex,colnames(tempData))
      if(length(whichIndex) == 0){
        whichIndex <- "allPairwise"  
        warning("Warning: no matches for index column. Using all pairwise instead.")
      } else if(length(whichIndex) > 1){
        whichMatch <- colnames(tempData)[whichIndex]
        warning("Warning: whichIndex returned more than one index column: ",paste(whichMatch,collapse="; "),". \n\t\tOnly ", whichMatch[1]," used.")
        whichIndex <- whichIndex[1]
      }
    } else if(is.numeric(whichIndex)){
      if(length(whichIndex) > 1){
        whichMatch <- colnames(tempData)[whichIndex]
        warning("Warning: whichIndex returned more than one index column: ",paste(whichMatch,collapse="; "),". \n\t\tOnly ", whichMatch[1]," used.")
        whichIndex <- whichIndex[1]
      } else if(whichIndex>length(names(tempData))){
        warning("Warning: no matches for index column. Using all pairwise instead.")
      }
    } else{
      whichIndex <- "allPairwise"
      warning("Warning: whichIndex must be numeric or character. Using all pairwise instead.")
    }
    
  }
  
  
  for(k in 1:length(names(tempData))){
    firstSet <- tempData[,k]
    firstName <- names(tempData)[k]
    
    if(whichIndex=="allPairwise"){
      for(j in 2:length(names(tempData))){
        if(j>k){
          secondSet <- tempData[,j]
          secondName <- names(tempData)[j]
          diffsFrame[,paste(firstName,"minus",secondName,appendName,sep="_")] <- firstSet-secondSet
        }
      }
    } else{
      if(k!=whichIndex){
        secondSet <- tempData[,whichIndex]
        secondName <- names(tempData)[whichIndex]
        diffsFrame[,paste(firstName,"minus",secondName,appendName,sep="_")] <- firstSet-secondSet
      }
    }
  }

  return(diffsFrame)
}
