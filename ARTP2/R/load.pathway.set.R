
load.pathway.set <- function(pathway, options){
  
  if(!is.list(pathway) && !is.character(pathway)){
    msg <- 'pathway should be a list of data frame, or a character vector'
    stop(msg)
  }
  
  if(is.character(pathway)){
    pathway <- as.list(pathway)
  }
  
  npath <- length(pathway)
  if(npath <= 1){
    msg <- 'at least two pathways should be specified in argument pathway'
    stop(msg)
  }
  
  path <- list()
  for(i in 1:npath){
    validate.pathway.definition(pathway[[i]])
    path[[i]] <- load.pathway.definition(pathway[[i]], options)
  }
  
  path
  
}
