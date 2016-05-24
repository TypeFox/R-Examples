#' Check invasive species status for a set of species from GISD database
#'
#' @export
#'
#' @param x character; a vector of scientific species names in the form of
#'    c("Genus species").
#' @param simplify logical; returns a data.frame with the species name and the
#'    values "Invasive", "Not in GISD". I recomend to check first the not
#'    simplified version (default), which contains raw information about the
#'    level of invasiveness.
#' @param verbose logical; If TRUE (default), informative messages printed.
#'
#' @return A list with species names, native range countries, and invasive
#' range countries
#'
#' @description This function check which species (both plants and animals) are
#' considered "invaders" somewhere in the world.
#'
#' For that end, it checks GISD (http://www.iucngisd.org/gisd) and
#' returns a value, either "Not in GISD" or the brief description presented in
#' GISD.
#'
#' Note that the webpage contains more information. Also note that the function
#' won't tell you if it's exotic in your area, a lot of exotic species are not
#' considered invaders (yet).
#'
#' As expected, the function is as good as the database is, which I find quite
#' reliable and well maintained.
#' The database is also able to recognize a lot (but not all) of the species
#' synonyms.
#'
#' Note that \code{eol_invasive} with source of gisd or gisd100 may end up with
#' different results as this function goes directly to the GISD website, whereas
#' eol_invasive only updates their GISD data occassionally. See notes in
#' \code{eol_invasive}.
#'
#' @author Ignasi Bartomeus \email{nacho.bartomeus@@gmail.com}
#' @examples \dontrun{
#' sp <- c("Carpobrotus edulis", "Rosmarinus officinalis")
#' ## first species is invasive, second one is not.
#' gisd(sp)
#' gisd(sp, simplify = TRUE)
#'
#' sp <- c("Carpobrotus edulis", "Rosmarinus officinalis", "Acacia mangium",
#' "Archontophoenix cunninghamiana", "Antigonon leptopus")
#' gisd(sp)
#' gisd(sp, simplify = TRUE)
#' }
gisd <- function(x, simplify = FALSE, verbose = TRUE) {
  outlist <- list()
  for (i in seq_along(x)) {
    mssg(verbose, paste("Checking", x[i]))
    out <- gbif_find(x[i])
    if (length(out) == 0) {
      outlist[[i]] <- list(species = x[i], status = "Not in GISD")
    } else{
      #Parse url and extract table
      doc <- xml2::read_html(paste0(gisd_base(), out$taxonID))
      if (!simplify) {
        alien <- gsub("^\\s+|\\s+$", "", gsub("\\[|\\]|[[:digit:]]", "", xml_text(xml_find_all(doc, '//div[@id="ar-col"]//ul/li'))))
        native <- gsub("^\\s+|\\s+$", "", xml_text(xml_find_all(doc, '//div[@id="nr-col"]//ul/li')))
        outlist[[i]] <- list(species = x[i], alien_range = alien, native_range = native)
      } else {
        outlist[[i]] <- list(species = x[i], status = "Invasive")
      }
    }
  }
  names(outlist) <- x
  mssg(verbose, "Done")
  return(outlist)
}

gbif_find <- function(x) {
  args <- list(datasetKey = 'b351a324-77c4-41c9-a909-f30f77268bc4', name = x)
  out <- GET('http://api.gbif.org/v1/species', query = args)
  stop_for_status(out)
  fromJSON(content(out, "text", encoding = "UTF-8"))$results
}

gisd_base <- function() "http://www.iucngisd.org/gisd/species.php?sc="
