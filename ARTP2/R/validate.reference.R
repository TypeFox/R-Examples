
validate.reference <- function(reference){
  
  # validate reference
  tmp <- (c("data.frame", "matrix") %in% class(reference))
  if(!any(tmp)){
    msg <- "pathway should be either an external file name or a data.frame"
    stop(msg)
  }else{
    if("matrix" %in% class(reference)){
      reference <- as.data.frame(reference)
    }
  }
  
  if(reference.type(reference) == 'ref.geno'){
    return(NULL)
  }
  
  header <- c("bed", "bim", "fam")
  tmp <- (header %in% colnames(reference))
  if(!any(tmp)){
    msg <- paste("Columns below were not found in reference:\n", paste(header[!tmp], collapse = " "))
    stop(msg)
  }
  
  reference <- reformat.reference.path(reference)
  
  tmp <- !file.exists(reference$bed)
  if(any(tmp)){
    msg <- paste(c("Files below were not found: ", reference$bed[tmp]), collapse = "\n")
    stop(msg)
  }
  
  tmp <- !file.exists(reference$bim)
  if(any(tmp)){
    msg <- paste(c("Files below were not found: ", reference$bim[tmp]), collapse = "\n")
    stop(msg)
  }
  
  tmp <- !file.exists(reference$fam)
  if(any(tmp)){
    msg <- paste(c("Files below were not found: ", reference$fam[tmp]), collapse = "\n")
    stop(msg)
  }
  
}
