list_files = function(from) 
  list.files(from, full.names = TRUE, pattern = "*.RData")


get_empty_results = function(from) {
  files = list_files(from)
  empty_files = rep(NA, length(files))
  for(i in seq_along(files))
    empty_files[i] = is.null(readRDS(files[i])$results)
  empty_files
}

#' Functions for manipulating uploaded results
#' 
#' Functions used for moving and creating the \code{past_results} data set from 
#' uploaded data. The \code{move_files} function is used to moved files from the server
#' to another location, whilst removing any empty data sets. 
#' @note One of the unit tests uploads an empty results file. 
#' Files where the results are \code{NULL} are moved to a sub-directory (called)
#' \code{empty} in the \code{to} directory. If the \code{empty} directory doesn't 
#' exist, it is created.
#' 
#' Currently these functions are specfic to my set-up.
#' @param from A directory containing the uploaded results.
#' @param to Destination directory
#' @export
move_files = function(from, to) {
  empty_files = get_empty_results(from)
  dir.create(paste0(to, "/empty/"), showWarnings = FALSE)
  to = ifelse(empty_files, paste0(to, "/empty/"), to)
  
  files = list_files(from)
  for(i in seq_along(files)) {
    cmd = paste("mv -v", files[i],  to[i])
    system(cmd)
  }
  invisible(files)
}
