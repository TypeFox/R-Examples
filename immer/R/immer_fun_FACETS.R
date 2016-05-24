# Functions for immer_FACETS

# 1. Check if no name is longer than 8 bits (and 3 bits for its extension)

grepInput <- function(name){
  name <- strsplit(name,c("=|;"))[[1]]
  name <- gsub("\\s","",name)
  name[ grep("\\w\\.\\w",name) ]
  
}

# rename names if they are longer than 8 Bit (without extensions)
bit8 <- function(x){
  splitTmp <- strsplit(x,"\\.")[[1]][1]
  tmp <- nchar(splitTmp)
  if(tmp>8){
    splitTmp <- strsplit(splitTmp,"|")[[1]][1:8]
    newName <- paste0(paste0(splitTmp,collapse = ""),".",strsplit(x,"\\.")[[1]][2])
    oldName <- x
    return(c(newName,oldName))
  }else{
    return(x)
  }
}

