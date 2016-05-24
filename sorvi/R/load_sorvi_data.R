#' load_sorvi_data
#'
#' Arguments:
#' @param data.id data ID to download (suffix before .rda). Investigate the contents of the url path to check data.ids
#' @param verbose verbose 
#'
#' Return:
#' @return translations 
#'
#' @examples # translations <- load_sorvi_data("translations")
#'
#' @export
#' @references
#' See citation("sorvi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @keywords utilities

load_sorvi_data <- function(data.id, verbose = TRUE) {

  # Circumvent warnings
  fi.en.maakunnat <- NULL
  kuntarajat.maa.shp <- NULL

  url <- ropengov_storage_path()
  filepath <- paste(url, "/louhos/", data.id, ".rda", sep = "")
  if (verbose) { message(paste("Loading ", filepath, sep = "")) }
  load(url(filepath))  

  if ( data.id == "translations" ) { return(fi.en.maakunnat) }
  if ( data.id == "kuntarajat.maa.shp" ) { return(kuntarajat.maa.shp) }

}

