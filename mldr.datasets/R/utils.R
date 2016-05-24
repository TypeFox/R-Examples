#' @title Check if an mldr object is locally available and download it if needed
#' @description This function checks if the mldr object whose name is given as input is locally available, loading it in memory.
#' If necessary, the dataset will be downloaded from the GitHub repository and saved locally.
#' @param mldr.name Name of the dataset to load
#' @examples
#'\dontrun{
#' library(mldr.datasets)
#' check_n_load.mldr("bibtex")
#' bibtex$measures
#' }
#' @export
check_n_load.mldr <- function(mldr.name) {
  if(exists(mldr.name, .GlobalEnv, mode = "list"))
    cat('The ', mldr.name, ' dataset is already loaded.')
  else {

    fpath <- paste(find.package('mldr.datasets'), '/extdata/', mldr.name, '.rda', sep = '')

    if (!file.exists(fpath)) {
      url <- availableMlds[availableMlds$Name == mldr.name, "URL"]

      if (length(url) > 0)
        download.file(url, fpath)
    }

    if (file.exists(fpath))
      load(fpath, .GlobalEnv)
    else
      stop('The ', mldr.name, ' dataset is not available. Try calling mldrs() to update the list of available datasets')
  }
}

#' @title Obtain and show a list of additional datasets available to download
#' @description The function downloads from GitHub the most up to date list of additional datasets. Those datasets are not
#' included into the package, but can be downloaded and saved locally.
#' @examples
#'\dontrun{
#' library(mldr.datasets)
#' mldrs()
#' }
#' @export
mldrs <- function() {
  availableMlds <- read.csv("https://fcharte.github.io/mldr.datasets/availableMlds.csv", stringsAsFactors = FALSE)
  save(availableMlds, file = paste0(find.package('mldr.datasets'), "/R/sysdata.rda"), compress = "gzip")

  View(availableMlds[-(length(availableMlds))], 'List of additional datasets available at the mldr.datasets repository')
}
