
validate.pathway.definition <- function(pathway){
  
  # validate pathway
  if(is.character(pathway)){
    if(length(pathway) > 1){
      msg <- "pathway should be either an external file name or a data.frame"
      stop(msg)
    }
    
    if(!file.exists(pathway)){
      msg <- paste0("File below was not found: \n", pathway)
      stop(msg)
    }
  }else{
    tmp <- (c("data.frame", "matrix") %in% class(pathway))
    if(!any(tmp)){
      msg <- "pathway should be either an external file name or a data.frame"
      stop(msg)
    }
  }
  
}
