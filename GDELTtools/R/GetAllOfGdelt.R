#' Download all the GDELT files to a local folder
#'
#' Downloads all (missing) GDELT files. ** This takes a long time and a lot of space. **
#' 
#' @aliases GetAllOfGDELT
#' @param local.folder character, path to the file to be validated.
#' @param data.url.root character, URL for the folder with GDELT data files.
#' @param force logical, if TRUE then the download is carried out without further prompting the user.
#' @return logical, TRUE if all files were downloaded successfully.
#' @export
#' @references
#' GDELT: Global Data on Events, Location and Tone, 1979-2012.  
#' Presented at the 2013 meeting of the International Studies Association
#' in San Francisco, CA.
#' \url{http://www.gdeltproject.org/}
#' @author 
#' \tabular{ll}{
#'   Stephen R. Haptonstahl \tab \email{srh@@haptonstahl.org}\cr
#' }
#' @examples
#' \dontrun{
#' GetAllOfGDELT("~/gdeltdata")} 
GetAllOfGDELT <- function(local.folder,
                          data.url.root="http://data.gdeltproject.org/events/",
                          force=FALSE) {
  
  if(FALSE == force) {
    # ask the user if they are sure
    size.of.GDELT <- GetSizeOfGDELT()
    w <- strwrap(paste("The compressed GDELT data set is currently ",
                       round(size.of.GDELT, 1),
                       "GB. It will take a long time to download and requires a lot of room (",
                       round(size.of.GDELT, 1),
                       "GB) where you store it. Please verify that you have sufficient free space on the drive where you intend to store it.",
                       sep=""))
    writeLines(w)
    response <- readline("Are you ready to proceed? (y/n) ")
    
    if(FALSE == grepl("[yY]", response)) {
      return(FALSE)
    }
  }
  
  # Coerce ending slashes as needed
  local.folder <- StripTrailingSlashes(path.expand(local.folder))
  data.url.root <- paste(StripTrailingSlashes(data.url.root), "/", sep="")
  # create the local.folder if is doesn't exist
  dir.create(local.folder, showWarnings=FALSE, recursive = TRUE)
  
  source.files <- ListAllGDELTFiles()$compressed
  
  res <- sapply(source.files, function(this.file) {
    this.res <- FALSE
    
    # validate if already downloaded
    if( this.file %in% LocalVersusRemote(filelist=source.files, local.folder=local.folder)$local ) {
      if(FALSE == IsValidGDELT(f=this.file, local.folder=local.folder)) {
        # remove the offending file; it'll be downloaded in the later if/then
        file.remove(paste(local.folder, "/", this.file, sep=""))
        Sys.sleep(1)
      } else {
        this.res <- TRUE
      }
    }
    # download if not already downloaded
    if( this.file %in% LocalVersusRemote(filelist=source.files, local.folder=local.folder)$remote ) {
      download.result <- DownloadGdelt(f=this.file,
                                       local.folder=local.folder,
                                       max.local.mb=Inf,
                                       data.url.root=data.url.root,
                                       verbose=TRUE)
      if(TRUE == download.result) {
        if(FALSE == IsValidGDELT(f=this.file, local.folder=local.folder)) {
          # try again
          download.result <- DownloadGdelt(f=this.file,
                                           local.folder=local.folder,
                                           max.local.mb=Inf,
                                           data.url.root=data.url.root,
                                           verbose=TRUE)
          if(TRUE == IsValidGDELT(f=this.file, local.folder=local.folder)) {
            this.res <- TRUE
          }
        } else {
          this.res <- TRUE
        }
      }
    }
    
    # return results for this.file
    if(TRUE == this.res) {
      cat("Downloading or verifying", this.file, "succeeded.\n")
    } else {
      cat("Downloading or verifying", this.file, "FAILED.\n")
    }
    return(this.res)
  })
  
  return( all(res) )
}
