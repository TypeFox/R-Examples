#' Get plant taxonomical and distribution data
#' 
#' This function collects taxonomic information and distribution from the
#' Brazilian Flora Checklist. Synonyms and misspelled names are resolved 
#' automatically. Results can be combined with life form, habitat, vernacular
#' name, and occurrence data.
#' 
#' @param taxa a character vector containing one or more taxa, without authors 
#'   see \code{\link{remove.authors}} if you have a list with authorities
#' @param replace.synonyms should the function automatically replace synonyms?
#' @param suggest.names should the function try to correct misspelled names?
#' @param life.form include the life form of the taxon?
#' @param habitat include the habitat of the taxon?
#' @param vernacular include vernacular names and localities?
#' @param states include occurrence data?
#' @param establishment include the establishment type (native, cultivated or 
#'   naturalized)?
#' @param drop NULL or character vector with names of columns with taxonomic
#'   information to be removed from the returned data frame. Available names: 
#'   "id", "scientific.name", "accepted.name", "family", "genus",
#'   "specific.epiteth", "infra.epiteth", "taxon.rank", "authorship",
#'   "taxon.status", "name.status", "threat.status", and "search.str".
#' @param suggestion.distance a value between 0 and 1 indicanting how conservative the
#'   name suggestion algorithm should be. Values closer to 1 are very
#'   conservative. Be very careful, lower values can give wrong suggestions.
#' @param parse Parse names through the GBIF parser to remove authors?
#' @details The returned data frame will contain a variable number of rows and 
#'   columns depending on how the function was called. For instance, since there
#'   might be more than one vernacular name for each taxon, some rows
#'   will be duplicated if \code{vernacular} is set to \code{TRUE}. All misspelled taxa
#'   are automatically corrected if the function can come up with a reasonable
#'   guess for the name. Conservation status follows the IUCN nomenclature.
#' @return a data frame
#' @export
#' @examples 
#' \dontrun{
#' data(plants)
#' get.taxa(plants)
#' get.taxa(plants, life.form = TRUE, establishment = TRUE)
#' }
get.taxa <- function (taxa, replace.synonyms = TRUE, suggest.names = TRUE, 
                      life.form = FALSE, habitat = FALSE, vernacular = FALSE, states = FALSE, 
                      establishment = FALSE, drop = c("authorship", "genus", "specific.epiteth", 
                                                      "infra.epiteth", "name.status"), 
                      suggestion.distance = 0.9, parse = FALSE) 
{
  taxa <- trim(taxa)
  taxa <- taxa[nzchar(taxa)]
  if (length(taxa) == 0L) 
    stop("No valid names provided.")
  original.search <- taxa
  ncol.taxa <- ncol(all.taxa)
  res <- data.frame(matrix(vector(), length(taxa), ncol.taxa + 
                             1, dimnames = list(c(), c(names(all.taxa), "notes"))), 
                    stringsAsFactors = FALSE)
  minus.notes <- seq_len(ncol.taxa)
  index <- 0
  for (taxon in taxa) {
    notes <- NULL
    index <- index + 1
    if (parse) {
      url <- "http://api.gbif.org/v1/parser/name"
      request <- try(POST(url, body = list(taxon), encode = "json"))
      if (inherits(request, "try-error")) {
        warning("Couldn't connect with the GBIF data servers. Check your internet connection or try again later.")
      } else {
        warn_for_status(request)
        taxon <- content(request)[[1]]$canonicalName
      }
    }
    taxon <- fixCase(taxon)
    uncertain <- regmatches(taxon, regexpr("[a|c]f+\\.", 
                                           taxon))
    if (length(uncertain) != 0L) {
      taxon <- gsub("\\s[a|c]f+\\.", "", taxon)
    }
    ident <- regmatches(taxon, regexpr("\\s+sp\\.+\\w*", 
                                       taxon))
    if (length(ident) != 0L) {
      split.name <- unlist(strsplit(taxon, " "))
      taxon <- split.name[1]
      infra <- split.name[2]
    }
    found <- length(with(all.taxa, {
      which(search.str == taxon)
    })) > 0L
    if (!found) {
      if (suggest.names) {
        taxon <- suggest.names(taxon, max.distance = suggestion.distance)
      }
      else {
        res[index, "notes"] <- "not found"
        next
      }
      if (is.na(taxon)) {
        res[index, "notes"] <- "not found"
        next
      }
      else {
        notes <- "was misspelled"
      }
    }
    accepted <- all.taxa[with(all.taxa, {
      which(search.str == taxon & taxon.status == "accepted")
    }), ]
    if (nrow(accepted) > 0) {
      if (nrow(accepted) == 1L) {
        res[index, minus.notes] <- accepted
      }
      else {
        notes <- c(notes, "check +1 accepted")
      }
      res[index, "notes"] <- paste(notes, collapse = "|")
      if (length(ident) != 0L) res[index, "search.str"] <- paste(res[index, "search.str"], infra)
      next
    }
    synonym <- all.taxa[with(all.taxa, {
      which(search.str == taxon & taxon.status == "synonym")
    }), ]
    nrow.synonym <- nrow(synonym)
    if (nrow.synonym > 0L) {
      if (replace.synonyms) {
        related <- relationships[with(relationships, {which(related.id %in% synonym$id)}), ]   
        accepted <- all.taxa[with(all.taxa, {
          which(id %in% related$id & taxon.status == "accepted")
        }), ]
        nrow.accepted <- nrow(accepted)
        if (nrow.accepted == 0L) {
          if (nrow.synonym == 1L) {
            notes <- c(notes, "check no accepted name")
            res[index, minus.notes] <- synonym
          }
          if (nrow.synonym > 1L) {
            notes <- c(notes, "check no accepted +1 synonyms")
          }
        }
        if (nrow.accepted == 1L) {
          notes <- c(notes, "replaced synonym")
          res[index, minus.notes] <- accepted
        }
        if (nrow.accepted > 1L) {
          notes <- c(notes, "check +1 accepted")
          if (nrow.synonym == 1L) {
            res[index, minus.notes] <- synonym
          }
        }
      }
      else {
        if (nrow(synonym) == 1L) {
          res[index, minus.notes] <- synonym
        }
        else {
          notes <- c(notes, "check +1 entries")
        }
      }
      res[index, "notes"] <- paste(notes, collapse = "|")
      if (length(ident) != 0L) res[index, "search.str"] <- paste(res[index, "search.str"], infra)
      next
    }
    
    undefined <- all.taxa[with(all.taxa, {
      which(search.str == taxon & is.na(taxon.status))
    }), ]
    
    nrow.undefined <- nrow(undefined)
    
    if (nrow.undefined == 0L) {
      notes <- c(notes, "check undefined status")
    }
    
    if (nrow.undefined == 1L) {
      notes <- c(notes, "check undefined status")
      res[index, minus.notes] <- undefined
    }
    
    if (nrow.undefined > 1L) {
      notes <- c(notes, "check undefined status")
    }
    
    res[index, "notes"] <- paste(notes, collapse = "|")
    if (length(ident) != 0L) res[index, "search.str"] <- paste(taxa, infra)
  }
  if (is.null(drop)) {
    res <- data.frame(res, original.search, stringsAsFactors = FALSE)
  }
  else {
    res <- data.frame(res[, !names(res) %in% drop], original.search, 
                      stringsAsFactors = FALSE)
  }
  if (life.form) {
    res <- left_join(res, species.profiles[, c("id", "life.form")], 
                     by = "id")
  }
  if (habitat) {
    res <- left_join(res, species.profiles[, c("id", "habitat")], 
                     by = "id")
  }
  if (vernacular) {
    res <- left_join(res, vernacular.names[, c("id", "vernacular.name")],
                     by = "id")  
  }
  if (states) {
    res <- left_join(res, distribution[, c("id", "occurrence")], 
                     by = "id")
  }
  if (establishment) {
    res <- left_join(res, distribution[, c("id", "establishment")], 
                     by = "id")
  }
  res
}