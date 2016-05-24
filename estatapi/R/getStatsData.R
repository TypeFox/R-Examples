#' getStatsData API
#'
#' Get some statistical data from e-stat API.
#'
#' @param appId application ID
#' @param statsDataId ID of the statistical dataset
#' @param ... Other parameters.
#' @seealso
#' \url{http://www.e-stat.go.jp/api/e-stat-manual/#api_2_3}
#' \url{http://www.e-stat.go.jp/api/e-stat-manual/#api_3_4}
#' @section Other parameters:
#' For every detailed information, please visit the URL in See Also.
#' \itemize{
#'  \item \code{lvTab}:
#'    Level of the meta-information. The format can be \code{X} or \code{X-Y}, \code{-X} and \code{X-}.
#'  \item \code{cdTab}:
#'    Code(s) of the meta-infomation items to select. The format can be a character vector (\code{c("001", "002")}) or
#'    a chraracter of values and commas (\code{"001,002"}).
#'  \item \code{cdTabFrom}:
#'    The code of the first meta-information item to select.
#'  \item \code{cdTabTo}:
#'    The code of the last meta-information item to select.
#'  \item \code{lvTime}:
#'     Level of the time to select. The format is the same as \code{lvTab}
#'  \item \code{cdTime}
#'     Time(s) to select. The format is the same way like \code{cdTab}
#'  \item \code{cdTimeFrom}:
#'    The first time to select. The format is the same way like \code{cdTabFrom}
#'  \item \code{cdTimeTo}:
#'    The last time to select. The format is the same way like \code{cdTabTo}
#'  \item \code{lvArea}:
#'     Level of the area to select. The format is the same as \code{lvTab}
#'  \item \code{cdArea}
#'     Code(s) of the Area to select. The format is the same way like \code{cdTab}
#'  \item \code{cdAreaFrom}:
#'    The code of the first area to select. The format is the same way like \code{cdTabFrom}
#'  \item \code{cdAreaTo}:
#'    The code of the last area to select. The format is the same way like \code{cdTabTo}
#'  \item \code{lvCat01}, \code{cdCat01}, \code{cdCat01From}, \code{cdCat01To}, ...:
#'    The same way like above.
#'  \item \code{startPosition}:
#'    integer. The the first record to get.
#'  \item \code{limit}:
#'    integer. Max number of records to get.
#' }
#'
#' @examples
#' \dontrun{
#' estat_getStatsData(
#'   appId = "XXXX",
#'   statsDataId = "0003065345",
#'   cdCat01 = c("008", "009", "010"),
#'   limit = 3
#' )
#' }
#'
#' @export
estat_getStatsData <- function(appId, statsDataId, ...)
{
  j <- estat_api("rest/2.0/app/json/getStatsData", appId = appId, statsDataId = statsDataId, ...)

  # TODO: rerun with startPosition automatically
  next_key <- j$GET_STATS_DATA$STATISTICAL_DATA$RESULT_INF$NEXT_KEY
  if(!is.null(next_key))
    message(sprintf("There are more records; please rerun with startPosition=%s", next_key))

  class_info <- get_class_info(j$GET_STATS_DATA$STATISTICAL_DATA$CLASS_INF$CLASS_OBJ)

  value_df <- j$GET_STATS_DATA$STATISTICAL_DATA$DATA_INF$VALUE %>%
    dplyr::bind_rows()

  suppressWarnings(
    value_df <- value_df %>%
      dplyr::mutate(value = readr::parse_number(`$`))
  )

  for (info_name in names(class_info)) {
    value_df <- merge_class_info(value_df, class_info, info_name)
  }

  value_df
}
