#' Search metadata for search terms using regex (more powerful than searching
#'    online without regex).
#' @import RCurl XML stringr
#' @param input Dryad metadata list, from e.g., getalldryad_metadata function,
#'    or load xml from directory (in which case, provide directory)
#' @param terms search terms e.g., 'plants', 'Whickam'
#' @param fuzzy (logical) do fuzzy search, TRUE (uses agrep) or FALSE (uses grep)
#' @param ignorecase (logical) if FALSE, pattern matching is case sensitive, and if
#'    TRUE, case is ignored during matching
#' @param value (logical) if FALSE, a vector containing integer (row) indices of the
#'    matches returned, and if TRUE, a vector containing the matching elements
#'    themselves is returned
#' @param maxdistance maximum distance allowed for a match. As integer, OR fraction
#'    of the pattern length, OR a list with possible entries:
#'    all (max. overall distance), insertions (max. number/fraction of
#'    insertions), deletions (max. number/fraction of deletions),
#'    and substitutions (max. number/fraction of substitutions)
#' @param loc where you want to search, any of title, creator, description, date, type,
#'    identifier, relation, OR 'all' for search over all metadata fields
#' @details Input is a Dryad metadata data frame from function getalldryad_metadata,
#'    or from directory (if latter, give path with arg 'input').
#' @return A numeric vector of OAI identifier's for datasets that match search.
#' @export
#' @examples \dontrun{
#' # Search data.frame in R
#' mymetdata <- getalldryad_metadata(T, progress = 'text', T, '/Mac/R_stuff/Blog_etc/Dryad/')
#' search_dryad(mymetdata, 'map', fuzzy=F, loc='type', maxdistance='all')
#' search_dryad(mymetdata, 'asddddf', fuzzy=T, loc='all')
#' search_dryad(mymetdata, 'clustal', fuzzy=F, ignorecase=T, value=F, loc='all')
#'
#' # Or search from a saved data.frame on file
#' search_dryad('/Mac/R_stuff/Blog_etc/Dryad/dryadmetadata.csv', 'me', fuzzy=T)
#' }
search_dryad <-
function(input, terms, fuzzy = "FALSE", ignorecase = "TRUE",
    value = "FALSE", maxdistance = 0.1, loc = "all") {
    searchdf <- function(x, terms) {
        identifier.1 <- "NA"
        rm(identifier.1)
        # The above two non-sense lines are to allow check() to ignore the global var warning.
        if (loc == "all") {
            if (fuzzy == "TRUE") {
                rowslist <- apply(x, 2, agrep, pattern = terms, ignore.case = ignorecase,
                  value = value, max.distance = maxdistance)
            } else if (fuzzy == "FALSE") {
                rowslist <- apply(x, 2, grep, pattern = terms, ignore.case = ignorecase,
                  value = value)
            }
        } else if (!loc == "all") {
            dat <- as.data.frame(str_match(names(x), "[a-z]+"))
            dat[, 1] <- as.character(dat[, 1])
            dat$rows <- rownames(dat)
            cols_ <- as.numeric(dat[dat[, 1] %in% loc, 2])
            if (fuzzy == "TRUE") {
                rowslist <- as.list(apply(data.frame(x[, cols_], x[,
                  cols_]), 2, agrep, pattern = terms, ignore.case = ignorecase,
                  value = value, max.distance = maxdistance))
            } else if (fuzzy == "FALSE") {
                rowslist <- as.list(apply(data.frame(x[, cols_], x[,
                  cols_]), 2, grep, pattern = terms, ignore.case = ignorecase,
                  value = value))
            }
        }

        if (class(try(do.call(c, rowslist), silent = T)) %in% "try-error") {
            stop("Awwwww snap. No datasets contain your search results")
        } else {
            rows <- do.call(c, rowslist)
        }
        rowsus <- sort(unique(rows))
        ids <- subset(x, rownames(x) %in% rowsus, identifier.1)
        oais <- as.numeric(apply(ids, 1, function(x) str_split(as.character(x),
            "dryad.")[[1]][2]))
        oais_ <- oais
    }

    if (class(input) == "data.frame") {
        oais_ <- searchdf(input, terms)
        oais_ <- oais_[!is.na(oais_)]
    } else if (class(input) == "character") {
        xdf <- read.csv(input)
        oais_ <- searchdf(xdf, terms)
        oais_ <- oais_[!is.na(oais_)]
    } else {
        stop("Error: input must be one of class data.frame or directory-file\nlocation\n or file name if in directory already")
    }
    oais_
}
