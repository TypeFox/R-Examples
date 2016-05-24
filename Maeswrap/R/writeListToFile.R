
writeListsToFile <- function(lis, fn){
  
  
  txt <- list()
  nml <- names(lis)
  for(i in 1:length(lis)){
    
    txt[[i]] <- formatNameList(nml[i], lis[[i]])
    txt[[i]] <- c(txt[[i]], "")
  }
  txt <- unlist(txt)
  
  writeLines(txt, fn)
}
