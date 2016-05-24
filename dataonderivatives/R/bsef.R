#' Get Bloomberg SEF data
#'
#' The Bloomberg Swap Execution Facility (SEF) offers customers the ability to
#' execute derivative instruments across a number of different asset classes.
#' It is required to make publicly available price, trading volume and other
#' trading data. It publishes this data on its website. I have reverse
#' engineered the JavaScript libraries used by its website to call the
#' Bloomberg Application Service using \code{POST} requests to a target URL.
#'
#' @param date the date for which data is required as Date or DateTime object.
#'   Only the year, month and day elements of the object are used. Must be of
#'   length one.
#' @param asset_class the asset class for which you would like to download
#'   trade data. Valid inputs are \code{"CR"} (credit), \code{"IR"} (rates),
#'   \code{"EQ"} (equities), \code{"FX"} (foreign exchange), \code{"CO"}
#'   (commodities). Can be a vector of these. Defaults to \code{NULL} which
#'   corresponds to all asset classes.
#' @param curate a logical flag indicating whether raw data should be returned
#'   or whether the raw data should be processed (default). The latter involves
#'   selecting particular fields and formatting these as seemed appropriate
#'   based on data and API reviews at the time the formatting was coded.
#' @return a data frame containing the requested data, or an empty data frame
#'   if data is unavailable
#' @importFrom dplyr %>%
#' @references \href{http://data.bloombergsef.com}{Bloomberg SEF data}
#' @examples
#' \dontrun{
#' library (lubridate)
#' # All asset classes
#' get_bsef_data(ymd(20140528))
#' # Only IR and FX asset classes
#' get_bsef_data(ymd(20140528), c("IR", "FX"))
#' }
#' @export

get_bsef_data <- function (date, asset_class = NULL, curate = TRUE) {
  valid_asset_classes <- c('CR', 'EQ', 'FX', 'IR', 'CO')
  if (is.null(asset_class)) {
    asset_class <- valid_asset_classes
  }
  assertthat::assert_that(all(asset_class %in% valid_asset_classes),
    lubridate::is.instant(date), length(date) == 1)
  df <- download_bsef_data(date, asset_class)
  if (!curate) {
    return (df)
  } else {
    return (df %>% format_bsef_data())
  }
}

download_bsef_data <- function (date, asset_class) {
  dplyr::bind_rows(Map(download_bsef_data_single, date, asset_class))
}

download_bsef_data_single <- function (date, asset_class) {
  # BAS doesn't appear to accept an end_date different from date: the
  # response is empty.
  start_date <- format(date, '%Y-%m-%dT00:00:00.000000Z')
  end_date <- start_date
  # Build POST body
  body <- list(list(tradeDays = list(startDay = start_date, endDay = end_date)))
  names(body) <- bsef_data_requestor(asset_class)
  body <- list(Request = body)
  response <- httr::POST(url = bsef_url(), config = bsef_header(), body = body,
    encode = 'json')
  # Convert response's content to JSON from raw
  response <- jsonlite::fromJSON(rawToChar(response$content))
  # Drill down response to data set that we are interested in
  df <- response$response[[bsef_data_responder(asset_class)]]$BsefEodData
  # Create asset_class field if necesary
  if (!is.null(df)) {
    if (is.list(df)) df <- dplyr::as_data_frame(df)
    df$assetclass <- asset_class
    return(df)
  } else {
    df <- dplyr::data_frame()
    return(df)
  }
}

#' @importFrom dplyr %>%
format_bsef_data <- function (df) {
  if (all(dim(df) == c(0, 0))) {
    return(dplyr::data_frame())
  } else {
    mutations <- list(
      ~lubridate::ymd_hms(tradeDate),
      ~as.numeric(priceOpen),
      ~as.numeric(priceHigh),
      ~as.numeric(priceLow),
      ~as.numeric(priceClose),
      ~as.numeric(settlementPrice),
      ~as.numeric(totalVolume),
      ~as.numeric(blockTradeVolume),
      ~as.numeric(totalVolumeUsd),
      ~as.numeric(blockTradeVolumeUsd)
    )
    mutation_names <- c("tradeDate", "priceOpen", "priceHigh", "priceLow",
      "priceClose", "settlementPrice", "totalVolume", "blockTradeVolume",
      "totalVolumeUsd", "blocktradevolumeusd")
    df %>%
      dplyr::mutate_(.dots = stats::setNames(mutations, mutation_names)) %>%
      dplyr::select_("tradeDate", "assetclass", "security", "currency",
        "priceOpen", "priceHigh", "priceLow", "priceClose", "settlementPrice",
        "totalVolume", "blockTradeVolume", "totalVolumeUsd",
        "blockTradeVolumeUsd")
  }
}

# URL target for data request
# Source: http://data.bloombergsef.com/assets/js/ticker.js
# Date accessed: 19 Sep 2014
bsef_url <- function () {
  'http://data.bloombergsef.com/bas/blotdatasvc'
}

# Bloomberg BAS version number
# Source: http://data.bloombergsef.com/assets/js/ticker.js
# Date accessed: 19 Sep 2014
bsef_header <- function (version = '1.9') {
  httr::add_headers('bas-version' = version)
}

# Way to tell Bloomberg which market data set is wanted
# Source: http://data.bloombergsef.com/assets/js/ticker.js
# Date accessed: 19 Sep 2014
bsef_data_requestor <- function (asset_class) {
  paste0('getBsefEod', bsef_asset_class_map(asset_class), 'DataRequest')
}

# Way to drill down into the data response that is provided
# Source: http://data.bloombergsef.com/assets/js/ticker.js
# Date accessed: 19 Sep 2014
bsef_data_responder <- function (asset_class) {
  paste0('getBsefEod', bsef_asset_class_map(asset_class), 'DataResponse')
}

bsef_asset_class_map <- function (asset_class) {
  asset_classes <- c(CR = 'Cds', EQ = 'Eqt', FX = 'Fx', IR = 'Irs', CO = 'Cmd')
  unname(asset_classes[asset_class])
}
