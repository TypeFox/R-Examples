calcCombVals <-
function(data,groupID,colID="FPKM",combineCols=TRUE,
                         matchEnd=TRUE,FUN="mean",functionName=NULL){
  ## groupID should either me vector of characters defining group, 
  ##   or a list of the column's for each group (e.g. list(Vfpkm=c(1,2,3),Cfpkm=(4,5,6)))
  ##   If using characters make sure the names are unique enough to match only what you want them to
  ## data should be matrix or data.frame
  ## colID tells which measure you want the mean of, can be a vector for several things
  ## combineCols tells whether or not to put the colums into one data.frame
  ## FUN should either be a function or a character to match to function
  ## If not a character, set functionName to assign column names for the output
  
  ## Name the function Set which function to use:
  if(is.null(functionName)){
    functionName <- as.character(FUN  )
  }
  FUN<-match.fun(FUN)
  
  colIndex <- list() # set a variable to hold the column index
  
  if(is.character(groupID)){
    for(j in 1:length(colID)){
      ## Set colID to "*" for searching, if it is blank
      if(colID[j]==""|colID[j]=="*"|colID[j]=="all"){
        colID[j]<-"*"
        colName <- "all"
      } else{
        colName <- colID[j]
      }
      
      indCols <- character()
      for(k in 1:length(groupID)){
        ## set groupID to "*" if blank, to match all instead of none
        if(groupID[k]==""|groupID[k]=="*"){
          groupID[k]<-"*"
          groupName <- "all"
        }else{
          groupName <- groupID[k]
        }
        if(matchEnd){
          indCols <- c(indCols,grep(paste(colID[j],"$",sep=""),colnames(data),value=TRUE) )  
        } else {
          indCols <- grep(colID[j],colnames(data),value=TRUE)
        }

        colIndex[[paste(groupName,colName,functionName,sep="_")]] <- grep(groupID[k],indCols,value=TRUE)
      } 
    } 
  } else if(is.list(groupID)) {
      colIndex <- groupID
  } else {
      stop("ERROR: colID must be either numeric or character")
  }
      
  temp <- lapply(colIndex,function(x) {apply(data[,x],1,FUN)})
    
  if(combineCols){
    tempOut <- as.data.frame(do.call(cbind,temp) )  
  } else{
    tempOut <- temp
  }
  
  return(tempOut)
}
