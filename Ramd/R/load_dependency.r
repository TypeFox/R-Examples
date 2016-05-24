#' Load a bunch of dependencies by filename
#'
#' @name load_dependency
#' @param dep Name of dependency, e.g., relative filename (without .r)
#' @examples
#' \dontrun{
#'   helper <- load_dependency('path/to/helper')
#' }
load_dependency <- function(dep) {
  path <- suppressWarnings(base::normalizePath(file.path(current_directory(), dep)))
  if (!file.exists(path)) {
    new_path <- paste(path, '.r', sep = '')
    if (!file.exists(new_path)) new_path <- paste(path, '.R', sep = '')
    path <- new_path
  }
  if (!file.exists(path))
    stop(paste("Unable to load dependency '", dep, "'", sep = ''))
  fileinfo <- file.info(path)
  mtime <- fileinfo$mtime

  value <- NULL
  
  if (path %in% get_src_cache_names() &&
      mtime == (cache_hit = get_src_cache(path))$mtime)
    value <- cache_hit$value
  else {
    # We fetch "source" from the global environment to allow other packages
    # to inject around sourcing files and be compatible with Ramd.
    value <- get('source', globalenv())(path)$value
    set_src_cache(list(value = value, mtime = mtime), path)
  }
  invisible(value)
}
