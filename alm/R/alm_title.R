#' Get title of article by inputting the doi for the article.
#'
#' @import httr
#' @importFrom stringr str_replace_all
#' @export
#' @param doi Digital object identifier for an article in PLoS Journals (character)
#' @param pmid PubMed object identifier (numeric)
#' @param pmcid PubMed Central object identifier (numeric)
#' @param wos Web of Science identifier (character)
#' @param scp Scopus identifier (character)
#' @param url Canonical URL (character)
#' @param key your PLoS API key, either enter, or loads from .Rprofile (character)
#' @param api_url API endpoint, defaults to http://alm.plos.org/api/v3/articles (character)
#' @param ... optional additional curl options (debugging tools mostly)
#' @return Title of article, in xml format.
#' @references See a tutorial/vignette for alm at
#' \url{http://ropensci.org/tutorials/alm_tutorial.html}
#' @examples \dontrun{
#' alm_title(doi='10.1371/journal.pbio.0000012')
#' dois <- c('10.1371/journal.pone.0026871','10.1371/journal.pone.0048868',
#'   	'10.1371/journal.pone.0048705','10.1371/journal.pone.0048731')
#' alm_title(doi=dois)
#' }

alm_title <- function(doi = NULL, pmid = NULL, pmcid = NULL, wos = NULL, scp = NULL, url = NULL,
  key = NULL, api_url = 'http://alm.plos.org/api/v5/articles', ...)
{
  .Deprecated(msg = "This function will be deprecated soon")
  temp <- alm_ids(doi = doi, pmid = pmid, pmcid = pmcid, wos = wos, scp = scp, url = url, api_url = api_url,
                  info = "summary", key = key, ...)
  if(length(doi) == 1) temp$data$info$title else lapply(temp$data, function(x) x$info$title)
}
