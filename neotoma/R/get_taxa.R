
#' Get taxon information from Neotoma.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @param taxonid Numeric taxon identifier used in Neotoma
#' @param taxonname A character string representing the full or partial name of taxa of interest.
#' @param status The current status of the taxon, one of 'extinct', 'extant', 'all'.
#' @param taxagroup The taxonomic grouping for the taxa. See \url{http://api.neotomadb.org/doc/resources/taxa} for the list of approved groupings.
#' @param ecolgroup The ecological group of the taxa. More detailed than \code{taxagroup}, can be obtained using \code{get_table("EcolGroupTypes")}.
#'
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return Returns a data frame with the following components:
#'
#'  \item{ \code{TaxonID} }{Unique database record identifier for a taxon}
#'  \item{ \code{TaxonCode} }{Shorthand notation for a taxon identification}
#'  \item{ \code{TaxonName} }{Name of the taxon}
#'  \item{ \code{Author} }{Author(s) of the name. Used almost exclusively with beetle taxa}
#'  \item{ \code{Extinct} }{True if extinct; false if extant}
#'  \item{ \code{TaxaGroup} }{Code for taxa group to which taxon belongs}
#'  \item{ \code{EcolGroups} }{Array of ecological group codes to which the taxon belongs}
#'  \item{ \code{HigherTaxonID} }{TaxonID of the next higher taxonomic rank}
#'  \item{ \code{PublicationID} }{Publication identification number}
#'  \item{ \code{Notes} }{Free-form notes or comments about the taxon}

#' @examples
#' \dontrun{
#' ## Return all species taxa with "Abies" in name - note wildcard
#' taxa <- get_taxa(taxonname = "Abies*")
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/contacts
#' @keywords IO connection
#' @export
get_taxa <- function(taxonid, taxonname, status, taxagroup, ecolgroup) {

  base.uri <- 'http://api.neotomadb.org/v1/data/taxa'

  cl <- as.list(match.call())
  cl[[1]] <- NULL
  cl <- lapply(cl, eval, envir = parent.frame())

  # Parameter check on taxagroup:
  if ('taxagroup' %in% names(cl)) {
    taxon.codes <- c('AVE', 'BIM', 'BRY',
                     'BTL', 'FSH', 'HRP',
                     'LAB', 'MAM', 'MOL',
                     'PHY', 'TES', 'VPL')

    if (!cl$taxagroup %in% taxon.codes) {
      stop(paste0('taxonGroup is not an accepted code. ',
                  'Use get_table(\'TaxaGroupTypes\') ',
                  'to obtain acceptible classes'))
    }
  }

  # Parameter check on taxonname and taxonids, I'm allowing
  # only one, but I think it can accept two.
  if (any(c('taxonids', 'taxonname') %in% names(cl))) {

    if (all(c('taxonids', 'taxonname') %in% names(cl))) {
      stop('Can only accept either taxonids OR taxonname, not both.')
    }
    if ('taxonids' %in% names(cl) & !is.numeric(cl$taxonids)) {
      stop(paste0('The variable taxonids must be numeric. ',
                  'To obtain a list of taxon IDs use the get_table command.'))
    }
    if ('taxonname' %in% names(cl) & !is.character(cl$taxonname)) {
      stop(paste0('The variable taxonname must be a character string. ',
                  'To obtain a list of taxon names use the get_table command.'))
    }
  }

  if ('status' %in% names(cl)) {
    if (!cl$status %in% c('extinct', 'extant', 'all')) {
      stop('Status must be one of: \'extinct\', \'extant\', or \'all\'')
    }
  }

  # Call the API:
  neotoma_content <- httr::content(httr::GET(base.uri, query = cl), as = "text")
  if (identical(neotoma_content, "")) stop("")
  aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)
  
  if (aa[[1]] == 0) {
    stop(paste('Server returned an error message:\n', aa[[2]]), call. = FALSE)
  }
  if (aa[[1]] == 1) {
    output <- aa[[2]]
    
    rep_NULL <- function(x) { 
      if (is.null(x)) {NA}
      else{
        if (class(x) == 'list') {
          # Recursive function to go through the list & clear NULL values.
          lapply(x, rep_NULL)
        } else {
          return(x)
        }
      }
    }
    
    # Clear NULLs from the output object & replace with NA values.
    output <- lapply(output, function(x)rep_NULL(x))
    
    cat('The API call was successful, you have returned ',
        length(output), 'records.\n')
  }

  if (class(aa) == 'try-error') {
    output <- aa
  } else {

      # Don't need anaonymous function here, call `[[()` with
      # argument "TaxonName". Equivalent of output[[i]][, "TaxonName"]
      names(output) <- sapply(output, `[[`, "TaxonName")

      # There are some values in here that are empty lists:
      output <- lapply(output, function(x) {
          len <- sapply(x, length) == 0
          if (any(len)) {
              x[[which(len)]] <- NA
          }
          x
      })

      # Bind each list into a single data frame
      output <- do.call(rbind.data.frame, output)
      rownames(output) <- NULL
  }

  output
}
