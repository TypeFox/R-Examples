
##' Sort a list alphabetically or given names
##'
##' NA
##' @title Sort a list alphabetically or given names
##' @param x 
##' @param xnames 
##' @return A sorted list
##' @author Xiaobei Zhao
sortList <- function(x, xnames=NULL){
  xnames0 <- names(x)
  if (is.null(xnames)){
    xnames <- sort(xnames0)
  } 
  xnames <- unique(xnames)
  if(! all(xnames %in% xnames0)){
    stop(sprintf('xnames should be within: ',paste(xnames,sep="",collapse=",")))
  }

  ret <- c()
  for (a in xnames){
    ret <- c(ret,x[xnames0==a])
  }
  return(ret)
} 


##' Covert a list to a name-content table
##'
##' NA
##' @title Covert a list to a name-content table
##' @param x 
##' @param sep 
##' @return A string
##' @author Xiaobei Zhao
list2table.string <- function(x,sep="\t"){
  ##print(x)
  ret=""
  for (i in 1:length(x)){
    ret <- paste(ret,sprintf("%s%s%s\n",names(x)[i],sep,x[[i]]),sep="")
  }
  return(ret)
}



##' Covert a matrix to a human-readable table
##'
##' NA
##' @title Covert a matrix to a human-readable table
##' @param x 
##' @param sep 
##' @param row.names 
##' @param col.names 
##' @return A string
##' @author Xiaobei Zhao
matrix2table.string <- function(x,sep=" ",row.names=TRUE,col.names=FALSE){
  x <- cbind(rownames(x),"[",x,"]")
  x <- apply(x,2,format,justify="right")
  paste(apply(x,1,paste,sep="",collapse=sep),sep="",collapse="\n")
}



##' Write a string to file
##'
##' NA
##' @title Write a string to file
##' @param outFpath 
##' @param entry 
##' @param mode 
##' @return NULL
##' @author Xiaobei Zhao
writeString <- function(outFpath,entry,mode="w"){
  input <- file(outFpath,open=mode)
  writeLines(entry, input, sep="")
  close(input)
}
