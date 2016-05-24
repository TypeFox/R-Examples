#' Definition of list of wavebands used in plant biology
#'
#' Defined according to different authors.
#'
#' @param std a character string "sensory", "sensory10", "sensory20" or ""
#'
#' @return a list of wavebands
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' Plant_bands()
#' Plant_bands("sensory")
#' Plant_bands("sensory10")
#' Plant_bands("sensory20")
#' Plant_bands("")
#'
#' @family lists of unweighted wavebands
#'
Plant_bands <- function(std = "sensory20") {
  if (std %in% c("sensory", "sensory10", "sensory20")) {
    if (std == "sensory10") {
      RFRstd <- "Smith10"
    } else {
      RFRstd <- "Smith20"
    }
    return(list(UVB(), UVA(), Blue("Sellaro"), Green("Sellaro"),
                Red(RFRstd), Far_red(RFRstd)))
  } else if (std %in% c("energy", "")) {
    return(list(UVB(), UVA(), PAR()))
  } else {
    stop("Unrecognized std: '", std, "'")
  }
}

