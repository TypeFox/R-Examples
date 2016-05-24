#' Easy http
#'
#' @importFrom methods is
#' @importFrom stats setNames
#' @importFrom utils head modifyList packageVersion
#' @importFrom R6 R6Class
#' @importFrom curl curl_version
#' @import httr
#' @name request-package
#' @aliases request
#' @docType package
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#' @keywords package
#'
#' @examples \dontrun{
#' ## Build API routes
#' ### Works with full or partial URLs
#' api('https://api.github.com/')
#' api('http://api.gbif.org/v1')
#' api('api.gbif.org/v1')
#'
#' ### Works with ports, full or partial
#' api('http://localhost:9200')
#' api('localhost:9200')
#' api(':9200')
#' api('9200')
#'
#' ## The above are not passed through a pipe, so simply define a URL, but don't
#' ## do a request. To make an http request, you can either pipe a url or
#' ## partial url to e.g., \code{\link{api}}, or call \code{\link{http}}
#' 'https://api.github.com/' %>% api()
#' ### Or
#' api('https://api.github.com/') %>% http()
#'
#' # Non-standard evaluation (NSE)
#' api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, issues) %>%
#'   peep
#'
#' # Standard evaluation (SE)
#' api('https://api.github.com/') %>%
#'   api_path_('repos', 'ropensci', 'rgbif', 'issues') %>%
#'   peep
#'
#' ## Templating
#' repo_info <- list(username = 'craigcitro', repo = 'r-travis')
#' api('https://api.github.com/') %>%
#'   api_template(template = 'repos/{{username}}/{{repo}}/issues', data = repo_info) %>%
#'   peep
#' }
NULL
