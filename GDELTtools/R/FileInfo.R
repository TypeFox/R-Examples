FileInfo <- function(path, include.dirs=FALSE) {
  # Returns info on files in path in a data.frame
  info.on.files <- ldply(paste(path, "/", dir(path), sep=""), file.info)
  info.on.files <- data.frame(name=dir(path), info.on.files, stringsAsFactors=FALSE)
  if(!include.dirs) info.on.files <- info.on.files[!info.on.files$isdir,]
  return(info.on.files)
}
