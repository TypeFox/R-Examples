#' Load a bunch of packages
#' 
#' @param ... List of packages
#' @export
#' @examples
#' \dontrun{
#' packages('glmnet', 'caret)
#' packages('glmnet caret')
#' packages(list('glmnet', 'caret'), c('e1071', 'parallel multicore'), 'stringr')
#' }
packages <- function(...) {
  split_packages <- function(string) strsplit(string, '[^a-zA-Z.$0-9_]+')[[1]]
  pkgs <- unlist(lapply(unlist(list(...)), split_packages))
  sapply(pkgs, load_package)
  TRUE
}
