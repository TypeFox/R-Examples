#' Download a file, using http, https, or ftp
#' 
#' This is a wrapper for \code{\link{download.file}} and takes all the same 
#' arguments. The only difference is that, if the protocol is https, it changes 
#' some settings to make it work. How exactly the settings are changed differs 
#' among platforms.
#' 
#' This function also should follow http redirects on all platforms, which is 
#' something that does not happen by default when \code{curl} is used, as on Mac
#' OS X.
#' 
#' With Windows, it either uses the \code{"wininet"} method (for R 3.2) or uses 
#' the \code{"internal"} method after first ensuring that \code{setInternet2}, 
#' is active (which tells R to use the \code{internet2.dll}).
#' 
#' On other platforms, it will try to use \code{libcurl}, \code{wget}, then 
#' \code{curl}, and then \code{lynx} to download the file. R 3.2 will typically 
#' have the \code{libcurl} method and for previous versions of R Linux platforms
#' will have \code{wget} installed, and Mac OS X will have \code{curl}.
#' 
#' Note that for many (perhaps most) types of files, you will want to use 
#' \code{mode="wb"} so that the file is downloaded in binary mode.
#'
#' @param url The URL to download.
#' @param ... Other arguments that are passed to \code{\link{download.file}}.
#'
#' @seealso \code{\link{download.file}} for more information on the arguments
#'   that can be used with this function.
#'
#' @export
#' @examples
#' \dontrun{
#' # Download the downloader source, in binary mode
#' download("https://github.com/wch/downloader/zipball/master",
#'          "downloader.zip", mode = "wb")
#' }
#'
#' @importFrom utils download.file
download <- function(url, ...) {
  # First, check protocol. If http or https, check platform:
  if (grepl('^https?://', url)) {
    
    # Check whether we are running R 3.2
    isR32 <- getRversion() >= "3.2"
    
    # Windows
    if (.Platform$OS.type == "windows") {
      
      if (isR32) {
        method <- "wininet"
      } else {
        
        # If we directly use setInternet2, R CMD CHECK gives a Note on Mac/Linux
        seti2 <- `::`(utils, 'setInternet2')
        
        # Check whether we are already using internet2 for internal
        internet2_start <- seti2(NA)
        
        # If not then temporarily set it
        if (!internet2_start) {
          # Store initial settings, and restore on exit
          on.exit(suppressWarnings(seti2(internet2_start)))
          
          # Needed for https. Will get warning if setInternet2(FALSE) already run
          # and internet routines are used. But the warnings don't seem to matter.
          suppressWarnings(seti2(TRUE))
        }
        
        method <- "internal"
      }
      
      # download.file will complain about file size with something like:
      #       Warning message:
      #         In download.file(url, ...) : downloaded length 19457 != reported length 200
      # because apparently it compares the length with the status code returned (?)
      # so we supress that
      suppressWarnings(download.file(url, method = method, ...))
      
    } else {
      # If non-Windows, check for libcurl/curl/wget/lynx, then call download.file with
      # appropriate method.
      
      if (isR32 && capabilities("libcurl")) {
        method <- "libcurl"
      } else if (nzchar(Sys.which("wget")[1])) {
        method <- "wget"
      } else if (nzchar(Sys.which("curl")[1])) {
        method <- "curl"

        # curl needs to add a -L option to follow redirects.
        # Save the original options and restore when we exit.
        orig_extra_options <- getOption("download.file.extra")
        on.exit(options(download.file.extra = orig_extra_options))

        options(download.file.extra = paste("-L", orig_extra_options))

      } else if (nzchar(Sys.which("lynx")[1])) {
        method <- "lynx"
      } else {
        stop("no download method found")
      }

      download.file(url, method = method, ...)
    }

  } else {
    download.file(url, ...)
  }
}
