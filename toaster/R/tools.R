#' Give a deprecation error, warning, or messsage, depending on version number.
#'
#' Version numbers have the format <major>.<minor>.<subminor>, like 0.9.2.
#' This function compares the current version number of toaster against the
#' specified \code{version}, which is the most recent version before the
#' function (or other object) was deprecated.
#'
#' \code{toa_dep} will give an error, warning, or message, depending on the
#' difference between the current toaster version and the specified
#' \code{version}.
#'
#' If the current major number is greater than \code{version}'s major number,
#' or if the current minor number is more than 1 greater than \code{version}'s
#' minor number, give an error.
#'
#' If the current minor number differs from \code{version}'s minor number by
#' one, give a warning.
#'
#' If the current subminor number differs from \code{version}'s subminor
#' number, print a message.
#'
#' @param version The last version of toaster where this function was good
#'   (in other words, the last version where it was not deprecated).
#' @param msg The message to print.
#' @seealso \link[ggplot2]{gg_dep}
#' @export
toa_dep <- function(version, msg) {
  v <- as.package_version(version)
  cv <- utils::packageVersion("toaster")
  
  # If current major number is greater than last-good major number, or if
  #  current minor number is more than 1 greater than last-good minor number,
  #  give error.
  if (cv[[1,1]] > v[[1,1]]  ||  cv[[1,2]] > v[[1,2]] + 1) {
    stop(msg, " (Defunct; last used in version ", version, ")",
         call. = FALSE)
    
    # If minor number differs by one, give warning
  } else if (cv[[1,2]] > v[[1,2]]) {
    warning(msg, " (Deprecated; last used in version ", version, ")",
            call. = FALSE)
    
    # If only subminor number is greater, give message
  } else if (cv[[1,3]] > v[[1,3]]) {
    message(msg, " (Deprecated; last used in version ", version, ")")
  }
  
  invisible()
}


toaSqlQuery <- function(channel, sql, errors = FALSE, ..., closeOnError=FALSE) {
  
  rs = sqlQuery(channel, sql, errors=errors, ...)
  
  if (length(rs) == 1 && is.numeric(rs) && rs == -1) {
    msg = odbcGetErrMsg(channel)
    if (closeOnError) 
      close(channel)
    stop(msg)
  }else
    return(rs)
}