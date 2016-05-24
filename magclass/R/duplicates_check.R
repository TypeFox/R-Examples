.duplicates_check <- function(coord) {
  coord <- as.data.frame(coord)
  duplicates <- duplicated(coord)
  if(any(duplicates)) {
    warning("Duplicate entries found, only the last entry will be used (duplicate entries: ",paste(apply(rbind(NULL,unique(coord[duplicates,])),1,paste,collapse="|"),collapse=", "),")!")    
  } 
}