#' Get the date when the article was published.
#'
#' @import httr
#' @importFrom stringr str_split
#' @export
#' @param doi Digital object identifier for an article in PLoS Journals (character)
#' @param pmid PubMed object identifier (numeric)
#' @param pmcid PubMed Central object identifier (numeric)
#' @param wos Web of Science identifier (character)
#' @param scp Scopus identifier (character)
#' @param url Canonical URL (character)
#' @param get Get year, month, or day; if unspecified, whole date returned.
#' @param key your PLoS API key, either enter, or loads from .Rprofile (character)
#' @param api_url API endpoint, defaults to http://alm.plos.org/api/v3/articles (character)
#' @param ... optional additional curl options (debugging tools mostly)
#' @references See a tutorial/vignette for alm at
#' \url{http://ropensci.org/tutorials/alm_tutorial.html}
#' @return Date when article was published.
#' @examples \dontrun{
#' alm_datepub(doi='10.1371/journal.pone.0026871')
#' alm_datepub('10.1371/journal.pone.0026871', get='year')
#'
#' # Provide more than one DOI
#' dois <- c('10.1371/journal.pone.0026871','10.1371/journal.pone.0048868',
#' 		'10.1371/journal.pone.0048705','10.1371/journal.pone.0048731')
#' alm_datepub(doi=dois, get="month")
#' }

alm_datepub <- function(doi = NULL, pmid = NULL, pmcid = NULL, wos = NULL, scp = NULL, url = NULL,
  get = NULL, key = NULL, api_url = 'http://alm.plos.org/api/v5/articles', ...)
{
  temp <- alm_ids(doi = doi, pmid = pmid, pmcid = pmcid, wos = wos, scp = scp, url = url, api_url = api_url,
                  info = "summary", key = key, ...)
  if(length(doi) == 1) getdate(temp$data, get) else lapply(temp$data, getdate, get=get)
}

getdate <- function(x, get) {
  date <- x$info$issued
  if(is.null(get)) { date } else {
    switch(get,
           year = as.numeric(strsplit(date,"-")[[1]][1]),
           month = as.numeric(strsplit(date,"-")[[1]][2]),
           day = as.numeric(strsplit(date,"-")[[1]][3])
    )
  }
}
