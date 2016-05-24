calcLogVal <-
function(data,colID="FPKM",offset=0.1,setBase=2,matchEnd=TRUE){
  ## If using colID make sure the names are unique enough to match only what you want them to
  ##    colID assumes that the ID is at the end of the line, unless matchEnd=FALSE
  ## data should be matrix or data.frame
  ## offset tells how much to add before taking log to avoid undefined
  ## setBase tells which base value to use
  colIndex <- NULL # set a variable to hold the column index
  if(is.character(colID)){
    for(k in 1:length(colID)){
      if(matchEnd){
        colIndex <- c(colIndex,grep(paste(colID[k],"$",sep=""),colnames(data)))  
      } else {
        colIndex <- c(colIndex,grep(colID[k],colnames(data)))  
      }
      
    } 
  } else if(is.numeric(colID)) {
    colIndex <- colID
  } else {
    stop("ERROR: colID must be either numeric or character")
  }
  
  if(length(colIndex)>1){ ## only run apply if you actually have a matrix, not just one vector
    temp <- as.data.frame(apply(data[,colIndex],c(1,2),function(x) {log(x+offset,base=setBase)}))  
    names(temp)<- paste(names(temp),"_log",setBase,sep="")
  } else {
    temp <- as.data.frame(log(data[,colIndex]+offset,base=setBase))
    names(temp) <- paste(names(data)[colIndex],"_log",setBase,sep="")
  }
  

  return(temp)
}
