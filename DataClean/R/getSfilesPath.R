#' If you want to collect all files under certain folder, this function should be the perfect one. It will collect all files with certain name. Then this function will return a list will all paths of those files so that further import or read is feasible.
#'
#' @title Collecting paths of some specified files that you want to import or read.
#'
#' @param root.path is the root path including all folders and files that you would like to search.
#' @param filename is the name of files that you want to collect.
#'
#' @return The whole paths of all files that meet the criteria were saved as a list.
#' @export
#' @examples
#' getSfilesPath(root.path = R.home(), filename = "?.exe")
getSfilesPath <- function(root.path, filename){
  if (missing(root.path)) {
    stop("Two 'files' must be a character string of name if it is saved in working directory, or it should include saving path of file.")
  }
  whole.path <- list.files(path = root.path, pattern = filename, full.names = TRUE, recursive = TRUE)
  return(whole.path)
}
