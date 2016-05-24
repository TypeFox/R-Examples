FirstTwo <- function(name) {
  if(is.null(name))
    name <- NA
  if(!is.na(name)) {
  	#trim2 <- function (x) gsub("^\\s+|\\s+$", "", x)
    name <- gsub("^\\s+|\\s+$", "", name)
    if(length(strsplit(name, " ")[[1]]) > 2)
      name <- paste(strsplit(name," ")[[1]][1:2],sep=" ",collapse=" ")
  }
  return(name)
}