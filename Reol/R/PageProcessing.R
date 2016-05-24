PageProcessing <- function(MyEOL, ...) {
  if(class(MyEOL) == "character" || class(MyEOL) == "vector")
    res <- xmlToList(xmlRoot(xmlParse(MyEOL, getDTD=FALSE), ...), simplify=FALSE)
  if(class(MyEOL) == "list")
    res <- xmlToList(xmlRoot(xmlParse(MyEOL[[1]], getDTD=FALSE), ...), simplify=FALSE)
  if(is.na(GetHierID(MyEOL)))
    return(paste("Filenames contain NAs"))
  if(!is.null(res$error)){
    paste("Bad file", MyEOL, "has an error:", res$error)
    return(NULL)
  }
  if(is.null(res$error))
    return(res)
}


RemoveNAFiles <- function(MyFiles){
  if(any(is.na(GetHierID(MyFiles)))) {
    whichNAs <- which(is.na(GetHierID(MyFiles)))
    MyFiles <- MyFiles[-whichNAs]
  }
  return(MyFiles)
}