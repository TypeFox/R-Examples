#' Append a command to the Log tab dealing with namespaces
#'
#' Time-stamp: <2015-11-15 10:03:29 gjw>
#'
#' @param comment      A message to include as a comment.
#' @param ...          The command(s) to report in the log.
#' @param include.libs Include any required library() calls.
#'
#' Report a command to the rattle Log tab textview. We check the
#' commands for any namespace usage and then include an appropriate
#' library() call for each and remove them from the commands
#' themselves.
#'
#' Each command will be printed on a new line.
#' 
appendLibLog <- function(comment, ..., include.libs=TRUE)
{
  # 150828 This started as the old appendLog but with a simplified
  # parameter list and added in the extraction of namespaces and then
  # rewrite the commands to not include the namespace.

  # Only continue if this is called from inside Rattle.
  
  if (is.null(crv$rattleGUI)) return()

  # Identify namespace string and namespace string with function.
  
  ns  <- '([a-zA-Z0-9_\\.]+)::'
  nsf <- stringr::str_c(ns, '([a-zA-Z0-9_\\.]+)')
  
  cmds <-
    list(...) %>%
    unlist() %>%
    stringr::str_c(collapse="\n")

  libs <-
    cmds %>%
    stringr::str_extract_all(nsf) %>%
    unlist() %>%
    unique() %>%
    stringr::str_split('::') %>%
    unlist()

  # 150917 Keep make check quiet....
  pkg <- fun <- funs <- "." <- NULL
  
  if (is.null(libs))
    include.libs <- FALSE
  else
    libs %<>%
      matrix(, ncol=2, byrow=TRUE) %>%
      data.frame(stringsAsFactors=FALSE) %>%
      magrittr::set_names(c("pkg", "fun")) %>%
      dplyr::group_by(pkg) %>%
      dplyr::summarise(funs=paste(fun, collapse="(), ")) %>%
      dplyr::group_by(pkg) %>%
      dplyr::summarise(cmd=sprintf("library(%s) # Provides %s().", pkg, funs)) %>%
      magrittr::extract2(2) %>%
      stringr::str_c(collapse="\n")

  cmds %<>%
    stringr::str_replace_all(ns, "")

  msg <-
    (if (include.libs) libs %s+% "\n\n" else "") %s+%
    cmds %>%
    paste(sep="", crv$start.log.comment, comment, crv$end.log.comment, .)

  # Always place the text at the end of the Log tab textview
  # irrespective of where the cursor is.

  log.buf <-
    theWidget("log_textview") $
    getBuffer()
  
  location <-
    log.buf $
    getEndIter() $
    iter
  
  log.buf $ insert(location, msg)
}

