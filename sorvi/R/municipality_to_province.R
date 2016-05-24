#' List province for each municipality in Finland.
#' @aliases municipality2province
#' @param municipalities NULL 
#' @param municipality.info NULL 
#' @return Mapping vector listing the province for each municipality in Finland.
#' @export 
#' @references
#' See citation("sorvi") 
#' @author Leo Lahti \email{louhos@@googlegroups.com}
#' @examples 
#' # Info table for municipalities:
#' # municipality.info <- get_municipality_info_mml()
#' # List all municipalities: 
#' # all.municipalities <- as.character(municipality.info$Kunta) 
#' # Pick province for given municipalities:
#' # mapping between municipalities (kunta) and provinces (maakunta)
#' # m2p <- municipality_to_province(c("Helsinki", "Tampere", "Turku")) 
#' # Speed up by providing predefined table of municipality info:
#' # m2p <- municipality_to_province(c("Helsinki", "Tampere", "Turku"), municipality.info)
#' @keywords utilities

municipality_to_province <- function (municipalities = NULL, municipality.info = NULL) {

  if (is.null(municipality.info)) { 
    municipality.info <- get_municipality_info_mml()
  }

  m2p <- as.character(municipality.info$Maakunta.FI)
  names(m2p) <- as.character(municipality.info$Kunta.FI)

  if (!is.null(municipalities)) {
    m2p <- m2p[as.character(municipalities)]
  }

  m2p

}


