#' @title Frequency of species from DAS file
#' @description Returns number of sightings of each species code in a DAS file.
#'
#' @param x filename of a DAS file or a \code{data.frame} from \code{\link{das.read}}.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
das.spp.freq <- function(x) {
  df <- if(is.data.frame(x)) x else das.read(x)
  df <- df[df$Event == "A" & df$OnEffort, ]
  spp <- unlist(df[, c("Spp1", "Spp2", "Spp3")])
  table(spp)
}
