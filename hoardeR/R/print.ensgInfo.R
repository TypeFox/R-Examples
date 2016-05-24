`print.ensgInfo` <- function(x, full=FALSE, ...){
  if(full){
    X <- data.frame(ENSG.ID=x[[1]], GeneName=x[[2]], GeneFull=x[[3]], GeneType=x[[4]], Summary=x[[5]], PubMedHits=x[[6]])
  } else {
    X <- data.frame(ENSG.ID=x[[1]], GeneName=x[[2]])    
  }

  print(X,...)
} 
