# Verify the integrity of a GDELT data file
#
# Compares the MD5 hash of a downloaded file to the known hash provided
# on the server.
# 
IsValidGDELT <- function(f,
                         local.folder) {
  
  md5.url <- "http://data.gdeltproject.org/events/md5sums"
  
  md5.df <- tryCatch(read.delim(md5.url, sep=" ", header=FALSE, stringsAsFactors=FALSE), 
                     error=function(e) stop(simpleError(paste("unable to read MD5 file at", md5.url), 
                                                        "IsValidGDELT")))
  
  this.md5 <- md5.df[ md5.df[,ncol(md5.df)]==f ,1]
  if(length(this.md5) != 1) {
    warning("Unable to find MD5 for ", f)
    return(FALSE)
  }
  
  observed.md5 <- tryCatch(md5sum(paste(StripTrailingSlashes(local.folder), "/", f, sep="")),
                           error=function(e) stop(simpleError("unable to calculate MD5 for downloaded file",
                                                              "IsValidGDELT")))
  
  return( observed.md5 == this.md5 )
}
