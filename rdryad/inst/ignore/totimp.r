#' Get Total Impact (http://totalimpact.org/) metrics or bibliography on a paper.
#' @import stringr RJSONIO plyr
#' @param id The DOI, URL, or PubMed ID - see examples.
#' @param fields One of 'metrics' or 'biblio'.
#' @param sleep Time (in seconds) before function sends API call - defaults to
#'    zero.  Set to higher number if you are using this function in a loop with
#'    many API calls.  Although, with TotalImpact, you can also just put in
#'    multiple id's in one call (much faster) instead of using a loop (slower).
#' @param url The base URL (do not change from default).
#' @return A list.
#' @export
#' @examples \dontrun{
#' totimp(id = '10.5061/dryad.8671')
#' totimp(id = '10.5061/dryad.8671', 'biblio')
#' totimp(id = list('10.5061/dryad.8671','18428094','http://opensciencesummit.com/program/'), 'biblio')
#' totimp(id = '18428094', 'biblio')
#' }
totimp <- function(id = NA, fields = "metrics", sleep = 0, 
    url = "http://total-impact.org/api/v1/items/") {
    Sys.sleep(sleep)
    id_ <- paste(laply(id, str_replace_all, pattern = "/", replacement = "%252F"), 
        collapse = ",")
    url2 <- paste(url, id_, ".json", "?fields=", fields, sep = "")
    fromJSON(url2)
}

# URL
#
#   http://total-impact.org/api/v1/items/18428094,10.1371%252Fjournal.pmed.0020124,http:%252F%252Fopensciencesummit.com%252Fprogram%252F,10.5061%252Fdryad.8048.json?fields=metrics
# Multiple
#
#   fromJSON('http://total-impact.org/api/v1/items/18428094,10.1371%252Fjournal.pmed.0020124,http:%252F%252Fopensciencesummit.com%252Fprogram%252F,10.5061%252Fdryad.8048.json?fields=metrics')
# Dryad
#
#   fromJSON('http://total-impact.org/api/v1/items/10.5061%252Fdryad.8671.json?fields=biblio')
# PubMed
#
#   fromJSON('http://total-impact.org/api/v1/items/10.1371%252Fjournal.pmed.0020124.json?fields=metrics')
# Other
#
#   fromJSON('http://total-impact.org/api/v1/items/18428094.json?fields=biblio') 
