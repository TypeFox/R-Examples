DownloadGdelt <- function(f,
                          local.folder,
                          max.local.mb,
                          data.url.root="http://data.gdeltproject.org/events/",
                          verbose=TRUE) {
  # Dowloads a single file, then removes files if necessary to get under max.local.mb
  # Returns TRUE if file downloaded successfully, FALSE otherwise
  # DOES NOT GIVE A WARNING if non-gdelt files are in the local.folder
  
  # Guardians
  if(!missing(max.local.mb)) stopifnot(max.local.mb >= 0)
  # Add guardian to ensure URLs end with a slash
  # Add guardian to ensure local.folder does NOT end with a slash or backslash
  
  if(missing(local.folder)) local.folder <- tempdir()
  # Coerce ending slashes as needed
  local.folder <- StripTrailingSlashes(path.expand(local.folder))
  data.url.root <- paste(StripTrailingSlashes(data.url.root), "/", sep="")
  
  # Download the file
  op <- options()
  options(HTTPUserAgent=paste("GDELTtools v", packageVersion("GDELTtools"),
                              " in ", getOption("HTTPUserAgent"),
                              sep=""))
  result <- download.file(url=paste(data.url.root, f, sep=""),
                          destfile=paste(local.folder, "/", f, sep=""),
                          quiet=!verbose)
  if(0 != result) return(FALSE)
  options(op)
  
  # Clean up if necessary
  if(!missing(max.local.mb)) {
    info.on.files <- FileInfo(local.folder)
    mb.currently.stored <- sum(info.on.files$size, na.rm=TRUE) / 2^20
    #browser()
    while(mb.currently.stored > max.local.mb) {
      # delete file in folder accessed longest ago, BUT NOT CURRENT FILE
      info.on.files <- info.on.files[-which(dir(local.folder, include.dirs=FALSE)==f),]  # remove current file from consideration for deletion
      info.on.files <- info.on.files[info.on.files$size > 0,]  # remove size-zero files
      if(0 == nrow(info.on.files)) {
        # exit, because current file is the only file
        mb.currently.stored <- 0
      } else {
        old.file.ids <- which(min(info.on.files$atime, na.rm=TRUE)==info.on.files$atime)
        if(length(old.file.ids) < 1) stop("No local files to delete.")
        else del.file.id <- old.file.ids[1]
        
        file.remove(paste(local.folder, "/", dir(local.folder, include.dirs=FALSE)[del.file.id], sep=""))
        
        # update
        info.on.files <- FileInfo(local.folder)
        mb.currently.stored <- sum(info.on.files$size, na.rm=TRUE) / 2^20
      }
    }
  }
  return(TRUE)
}
