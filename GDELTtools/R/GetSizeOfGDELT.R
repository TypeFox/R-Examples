GetSizeOfGDELT <- function(filesize.url="http://data.gdeltproject.org/events/filesizes") {
  
  # Returns the size of the complete data set, compressed, in GB
  
  fs <- read.delim(filesize.url, header=FALSE, sep=" ")
  names(fs) <- c("size.bytes", "file.name")
  fs$size.bytes <- as.numeric(fs$size.bytes)
  fs <- fs[fs$file.name %in% ListAllGDELTFiles()$compressed,]
  gb <- sum(fs$size.bytes) / (1024^3)
  
  return(gb)
}