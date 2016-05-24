exifTagNames <- function(inDir, whichSubDir = 1, returnMetadata = FALSE){

  stopifnot(file.exists(inDir))
  stopifnot(is.logical(returnMetadata))
  if(Sys.which("exiftool") == "") stop("cannot find ExifTool")
  dirs.tmp <- list.dirs(inDir, recursive = FALSE, full.names = TRUE)

  if(file.exists(dirs.tmp[whichSubDir]) == FALSE) stop("the specified subdirectory does not exist")

  file.tmp <- list.files(dirs.tmp[whichSubDir],
                         full.names = TRUE,
                         pattern    = ".JPG$|.jpg$",
                         recursive  = TRUE)[1]
  if(length(file.tmp) == 0) stop(paste("found no jpg in ", dirs.tmp[whichSubDir], sep = "\n"))

  if(returnMetadata == FALSE){
    command.tmp  <- paste('exiftool -csv "', file.tmp, '"', sep = "")
    metadata.tmp <- system(command.tmp, intern=TRUE)
    tagnames     <- sort(unlist(strsplit(metadata.tmp[[1]], split = ",")))
    return (tagnames)
  } else {
    command.tmp  <- paste('exiftool  "', file.tmp, '"', sep = "")
    metadata.tmp <- system(command.tmp, intern=TRUE)
    return (metadata.tmp)
  }
}