#' Suggest a valid name from a misspelled one
#' 
#' This function tries to suggest a valid name according to the Brazilian Flora
#' Checklist using a possibly incorrect one as a starting point.
#' 
#' @param taxon a character vector containing a single name
#' @param max.distance a numeric value indicating how conservative the function 
#'   should be when searching for suggestions. Values close to 1 are very 
#'   conservative
#' @param return.na a logical indicating whether to return a \code{NA} or the original 
#'   input when no suggestion is found
#' @param ignore.words \code{NULL} or a character vector with words to be ignored by the function. 
#'   Useful if you are automatizing a workflow and wants the function to ignore
#'   words or phrases such as "not found", "dead", "undetermined", and so on
#' @export
#' @return A character vector or \code{NA}
#' @examples
#' \dontrun{
#' suggest.names("Cofea arabyca")
#' suggest.names("Myrcia bela")
#' }
suggest.names <-
  function(taxon, max.distance = 0.75, return.na = TRUE, ignore.words = NULL) {
    #taxon <- iconv(taxon, to = "ASCII//TRANSLIT")
    taxon <- fixCase(taxon)
    taxon.orig <- taxon
    uncertain <- regmatches(taxon, regexpr("[a|c]f+\\.", taxon))
    taxon <- gsub("^\\s+|\\s+$", "", taxon)
    if (length(uncertain) != 0L) taxon <- gsub("[a|c]f+\\.", "", taxon)
    ident <- regmatches(taxon, regexpr("\\s+sp\\.+\\w*", taxon))
    if (length(ident) != 0L) taxon <- unlist(strsplit(taxon, " "))[1]
    if (!nzchar(taxon)) return(NA)
    first.letter <- strsplit(taxon, "")[[1]][1]
    species.first.letter <- all.taxa$search.str[grep(paste("^", first.letter, sep = ""), all.taxa$search.str)]
    l1 <- length(taxon)
    l2 <- length(species.first.letter)
    out <- adist(taxon, species.first.letter)
    distance <- 1 - (out/pmax(nchar(taxon), 
                                  nchar(species.first.letter)))
    max.dist <- max(distance, na.rm = TRUE)
    if (max.dist >= max.distance) {
      if (length(ident) == 0L) {
        res <- species.first.letter[distance == max(distance, na.rm = TRUE)][1]
        if (length(uncertain) == 0L) {
          return(res)
        } else {
          res <- unlist(strsplit(res, " "))
          return(paste(res[1], uncertain, res[2:length(res)]))
        }
      } else {
        paste(species.first.letter[distance == max(distance, na.rm = TRUE)][1], ident, sep = "")
      }
    } else {
      if (return.na) {
        NA
      } else {
        taxon.orig
      }
    }
  }
