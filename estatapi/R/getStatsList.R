#' getStatsList API
#'
#' Search for statistical datasets via e-Stat API.
#'
#' @param appId Application ID
#' @param searchWord Keyword for searching. You can use \code{OR} and \code{AND}. (e.g.: \code{apple AND orrange}).
#' @param ... Other parameters.
#' @seealso
#' \url{http://www.e-stat.go.jp/api/e-stat-manual/#api_2_1}
#' \url{http://www.e-stat.go.jp/api/e-stat-manual/#api_3_2}
#'
#' @section Other parameters:
#' For every detailed information, please visit the URL in See Also.
#' \itemize{
#'  \item \code{surveyYears}:
#'    Year and month when the survey was conducted. The format is either \code{YYYY}, \code{YYYYMM}, or \code{YYYYMM-YYYYMM}
#'  \item \code{openYears}:
#'    Year and month when the survey result was opened. The format is the same as \code{surveyYears}
#'  \item \code{statsField}:
#'    Field of statistics. The format is either two digits (large classification) or
#'    four digits (small classification). For the detail of the classification, see
#'    \url{http://www.soumu.go.jp/toukei_toukatsu/index/seido/sangyo/26index.htm}
#'  \item \code{statsCode}:
#'     Code assigned for each statistical agency and statistics. The format can be
#'     five digits (agency), and eight digits (statistics). For the detail, see
#'     \url{http://www.stat.go.jp/info/guide/public/code/code.htm}.
#'  \item \code{searchKind}:
#'     Type of statistics. \code{1}: summary, \code{2}: regional mesh, \code{3}: Sensus.
#'  \item \code{startPosition}:
#'    integer. The the first record to get.
#'  \item \code{limit}:
#'    integer. Max number of records to get.
#'  \item \code{updatedDate}:
#'    Last updated date. The format is either \code{YYYY}, \code{YYYYMM}, \code{YYYYMMDD}, \code{YYYYMMDD-YYYYMMDD}
#' }
#'
#' @examples
#' \dontrun{
#' estat_getStatsList(
#'   appId = "XXXX",
#'   searchWord = "CD",
#'   limit = 3
#' )
#' }
#' @export
estat_getStatsList <- function(appId, searchWord, ...) {
  j <- estat_api("rest/2.0/app/json/getStatsList", appId = appId, searchWord = searchWord, ...)

  j$GET_STATS_LIST$DATALIST_INF$TABLE_INF %>%
    purrr::map(as_flattened_character) %>%
    dplyr::bind_rows()
}
