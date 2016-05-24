getAssemblies <- function(xml){
  xmlOne <- xml[[1]]
  for(i in 1:length(xml)){
    xmlOne <- rbind(xmlOne,xml[[i]])
  }
  xmlOI <- xmlOne[,1]
  assemblies <- sapply(strsplit(xmlOI,", "),"[",2)
  assemblies
}