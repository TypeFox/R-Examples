#' Download and subset GDELT data
#'
#' Download the GDELT files necessary for a data set, import them,
#' filter on various crieteria, and return a data.frame. 
#' 
#' @aliases GetGDELT
#' @param start.date character, just about any human-readable form of the earliest date to include.
#' @param end.date character, just about any human-readable form of the latest date to include.
#' @param filter list, named list encoding the values to include for specified fields. See Details.
#' @param local.folder character, if specified, where downloaded files will be saved.
#' @param max.local.mb numeric, the maximum size in MB of the downloaded files that will be retained.
#' @param allow.wildcards logical, must be TRUE to use * in \code{filter} to specify 'any character(s)'.
#' @param use.regex logical, if TRUE then \code{filter} will be processed as a \code{\link{regular expression}}.
#' @param data.url.root character, URL for the folder with GDELT data files.
#' @param verbose logical, if TRUE then indications of progress will be displayed.
#' @return data.frame
#' @export
#' @details
#' 
#' If \code{local.folder} is not specified then downloaded files are stored in
#' \code{tempdir()}. If a needed file has already been downloaded to \code{local.folder}
#' then this file is used instead of being downloaded. This can greatly speed up future
#' 
#' 
#' Dates are parsed with \code{dateParse} in the TimeWarp package. 
#' Years must be given with four digits.
#' 
#' @section Filtering Results:
#' 
#' This is how you write the \code{filter}.
#' 
#' @references
#' GDELT: Global Data on Events, Location and Tone, 1979-2012.  
#' Presented at the 2013 meeting of the International Studies Association
#' in San Francisco, CA.
#' \url{http://www.gdeltproject.org/}
#' @author 
#' \tabular{ll}{
#'   Stephen R. Haptonstahl \tab \email{srh@@haptonstahl.org}\cr
#'   Thomas Scherer \tab \email{tscherer@@princeton.edu}\cr
#'   John Beieler \tab \email{jub270@@psu.edu}\cr
#' }
#' @examples
#' \dontrun{
#' test.filter <- list(ActionGeo_ADM1Code=c("NI", "US"), ActionGeo_CountryCode="US")
#' test.results <- GetGDELT(start.date="1979-01-01", end.date="1979-12-31", 
#'   filter=test.filter)
#' table(test.results$ActionGeo_ADM1Code)
#' table(test.results$ActionGeo_CountryCode)}
#' 
#' # Specify a local folder to store the downloaded files
#' \dontrun{
#' test.results <- GetGDELT(start.date="1979-01-01", end.date="1979-12-31", 
#'                          filter=test.filter,
#'                          local.folder="~/gdeltdata",
#'                          max.local.mb=500)}
GetGDELT <- function(start.date,
                     end.date=start.date,
                     filter,
                     local.folder=tempdir(),
                     max.local.mb=Inf,
                     allow.wildcards=FALSE, 
                     use.regex=FALSE,
                     data.url.root="http://data.gdeltproject.org/events/",
                     verbose=TRUE) {
  
  # Coerce ending slashes as needed
  local.folder <- StripTrailingSlashes(path.expand(local.folder))
  data.url.root <- paste(StripTrailingSlashes(data.url.root), "/", sep="")
  # create the local.folder if is doesn't exist
  dir.create(local.folder, showWarnings=FALSE, recursive = TRUE)
  
  start.date <- strftime(dateParse(start.date), format="%Y-%m-%d")
  end.date <- strftime(dateParse(end.date), format="%Y-%m-%d")
  
  if(start.date <= "2014-01-25" & end.date >= "2014-01-23") warning("Your date range includes some or all of Feb 23-25, 2014. GDELT data is not available for these dates.")
  
  out.initialized <- FALSE
  
  # Determine file list based on dates
  source.files <- FileListFromDates(start.date=start.date, end.date=end.date)
  if(0 == length(source.files)) stop("No GDELT files available for the dates specified")
  
  # Ingest and filter local files
  for(this.file in LocalVersusRemote(filelist=source.files, local.folder=local.folder)$local) {
    if(FALSE == IsValidGDELT(f=this.file, local.folder=local.folder)) {
      # remove the offending file; it'll be downloaded in the later loop
      file.remove(paste(local.folder, "/", this.file, sep=""))
      next
    }
    new.data <- GdeltZipToDataframe(f=paste(local.folder, "/", this.file, sep=""),
                                    daily=grepl("export.CSV", this.file, fixed=TRUE),
                                    verbose=verbose)
    if(!missing(filter)) new.data <- FilterGdeltDataframe(x=new.data,
                                                          filter=filter,
                                                          allow.wildcards=allow.wildcards,
                                                          use.regex=use.regex,
                                                          verbose=verbose)
    if(out.initialized) {
      out <- rbind(out, new.data)
    } else {
      out <- new.data
      out.initialized <- TRUE
    }
  }
  
  # Download, ingest, and filter remote files
  for(this.file in LocalVersusRemote(filelist=source.files, local.folder=local.folder)$remote) {
    download.result <- DownloadGdelt(f=this.file,
                                     local.folder=local.folder,
                                     max.local.mb=max.local.mb,
                                     data.url.root=data.url.root,
                                     verbose=verbose)
    if(FALSE == download.result) {
      stop("Unable to download file ", this.file, ". Please try again. If you get this result again, the file might not be available on the server.")
    }
    if(FALSE == IsValidGDELT(f=this.file, local.folder=local.folder)) {
      # try again
      download.result <- DownloadGdelt(f=this.file,
                                       local.folder=local.folder,
                                       max.local.mb=max.local.mb,
                                       data.url.root=data.url.root,
                                       verbose=verbose)
      if(FALSE == IsValidGDELT(f=this.file, local.folder=local.folder)) {
        stop("Unable to verify the integrity of ", this.file)
      }
    }
    
    new.data <- GdeltZipToDataframe(f=paste(local.folder, "/", this.file, sep=""),
                                    daily=grepl("export.CSV", this.file, fixed=TRUE),
                                    verbose=verbose)
    if(!missing(filter)) new.data <- FilterGdeltDataframe(x=new.data,
                                                          filter=filter,
                                                          allow.wildcards=allow.wildcards,
                                                          use.regex=use.regex,
                                                          verbose=verbose)
    if(out.initialized) {
      out <- rbind(out, new.data)
    } else {
      out <- new.data
      out.initialized <- TRUE
    }
  }
  
  # Filter one more time on dates
  start.date.numeric <- as.numeric(strftime(dateParse(start.date), format="%Y%m%d"))
  end.date.numeric <- as.numeric(strftime(dateParse(end.date), format="%Y%m%d"))
  out <- out[out$SQLDATE >= start.date.numeric & out$SQLDATE <= end.date.numeric,]
  
  return(out)
}
