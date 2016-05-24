#' @rdname is_r
#' @export
is_architect <- function()
{
  if(!"package:rj" %in% search() ||
    is.null(device_name <- formals(getOption("device"))$name) ||
    device_name != "rj.gd")
  {
    return(false("You are not running Architect/StatET."))
  }
  TRUE
}

# This is to avoid a NOTE by R CMD check.  See 
# http://stackoverflow.com/q/9439256/134830 and
# http://stackoverflow.com/q/8096313/134830
if(getRversion() >= "2.15.1")  utils::globalVariables("Revo.version")

#' @rdname is_r
#' @export
is_revo_r <- function()
{
  if(!exists("Revo.version", "package:base", inherits = FALSE) || 
    !is.list(Revo.version))
  {
    return(false("You are not running Revolution R."))
  }
  TRUE
}

#' @rdname is_r
#' @export
is_rstudio <- function()
{
  # Can also check 
  # gui <- .Platform$GUI
  # if(is.null(gui) || gui != "RStudio")
  # but this does not work when called from .Rprofile.site
  env <- Sys.getenv("RSTUDIO", "0")
  if(env != "1")
  {
    return(false("You are not running RStudio."))
  }
  TRUE
}

#' Is RStudio the current version?
#' 
#' Checks to see if the running version of RStudio is the current version.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_rstudio_current} returns \code{TRUE} or \code{FALSE}, and
#' \code{assert_is_rstudio_current} throws an error in the event of an out of
#' date RStudio.  Non-RStudio IDEs throw an error.
#' @references This function is engineered from the \code{downloadUpdateInfo} 
#' function from 
#' \url{https://github.com/rstudio/rstudio/blob/master/src/cpp/session/modules/SessionUpdates.R}
#' where the string for the OS is described in \code{beginUpdateCheck} from
#' \url{https://github.com/rstudio/rstudio/blob/master/src/cpp/session/modules/SessionUpdates.cpp}
#' @seealso \code{\link{is_rstudio}}, \code{\link{is_rstudio_desktop}}
#' @export
is_rstudio_current <- function()
{
  assert_is_rstudio()
  os <- if(is_windows()) "windows" else if(is_mac()) "mac" else "linux"
  current_version <- rstudio_version_info()$version
  if(is.na(current_version))
  {
    return(
      false(
        "RStudio is out of date; it is old enough that the version check API has changed."
      )
    ) 
  }
  update_url <- sprintf(
    "http://www.rstudio.org/links/check_for_update?version=%s&os=%s&format=kvp",
    current_version,
    os
  )
  lines <- readLines(update_url, warn = FALSE)
  update_version <- substring(
    grep("update-version=([0-9.]*)", strsplit(lines, "&")[[1]], value = TRUE),
    16
  )
  if(nzchar(update_version))
  {
    return(
      false(
        "RStudio is out of date; you are running %s but version %s is available.",
        current_version,
        update_version
      )
    )
  }
  TRUE
}

#' Is RStudio running in desktop or server mode?
#' 
#' Checks for RStudio desktop or server version.
#' @references The values that RStudio uses for its mode are defined in
#' \url{https://github.com/rstudio/rstudio/blob/master/src/cpp/session/include/session/SessionConstants.hpp}
#' via the constants \code{kSessionProgramModeDesktop} and 
#' \code{kSessionProgramModeServer}.
#' @seealso \code{\link{is_rstudio}}, \code{\link{is_rstudio_current}}
#' @examples 
#' is_rstudio_desktop()
#' is_rstudio_server()
#' @export
is_rstudio_desktop <- function()
{
  if(!(ok <- is_rstudio()))
  {
    return(false(cause(ok)))
  }
  info <- rstudio_version_info()
  if(!is.list(info) && is.na(info)) # very old RStudio
  {
    return(info)
  }
  if(info$mode != "desktop")
  {
    return(false(gettext("You are running the server version of RStudio.")))
  }
}

#' @rdname is_rstudio_desktop
#' @export
is_rstudio_server <- function()
{
  if(!(ok <- is_rstudio()))
  {
    return(false(cause(ok)))
  }
  info <- rstudio_version_info()
  if(!is.list(info) && is.na(info)) # very old RStudio
  {
    return(info)
  }
  if(info$mode != "server")
  {
    return(false(gettext("You are running the desktop version of RStudio.")))
  }
}

#' Get RStudio's version information
#' 
#' Wrapper to .rs.api.versionInfo.
rstudio_version_info <- function()
{
  assert_is_rstudio()
  tools_rstudio <- "tools:rstudio"
  if(!tools_rstudio %in% search())
  {
    return(
      na(
        gettext(
          "'tools:rstudio' is not loaded, so the RStudio API is not available."
        )
      )
    )
  }
  e <- as.environment(tools_rstudio)
  if(!".rs.api.versionInfo" %in% ls(e, all.names = TRUE))
  {
    return(
      na(
        gettext(
          "You are using an old version of RStudio, which does not tell you version information."
        )
      )
    )
  }
  e$.rs.api.versionInfo()
}

