#' Get Bloomberg SDR data
#'
#' The Bloomberg Swap Data Repository (BSDR) is a registered U.S. swap data
#' repository that allows market participants to fulfil their public disclosure
#' obligations under U.S. legislation. BSDR is required to make publicly
#' available price, trading volume and other trading data reported to its U.S.
#' repository. It publishes this data on its website in real-time and also on a
#' historical basis. I have reverse engineered the JavaScript libraries used by
#' its website to call the Bloomberg Application Service using \code{POST}
#' requests to a target URL.
#'
#' @param date the date for which data is required as Date or DateTime object.
#'   It will use all date-time elements including year, month, day, hour,
#'   minute, second (incl. fractional seconds) and  time zone information to
#'   determine the set of trades to return. It will return the set of trades
#'   for the day starting on \code{date}.
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
#' @references \href{http://www.bloombergsdr.com/search}{BSDR search}
#' @examples
#' \dontrun{
#' library (lubridate)
#' # All asset classes for day starting 6 May 2015
#' get_bsdr_data(ymd(20150506))
#' # Only IR and FX asset classes
#' get_bsdr_data(ymd(20150506), c("IR", "FX"))
#' }
#' @export

get_bsdr_data <- function (date, asset_class = NULL, curate = TRUE) {
  valid_asset_classes <- c('CR', 'EQ', 'FX', 'IR', 'CO')
  if (is.null(asset_class)) {
    asset_class <- valid_asset_classes
  }
  assertthat::assert_that(all(asset_class %in% valid_asset_classes),
    lubridate::is.instant(date), length(date) == 1)
  df <- download_bsdr_data(date, asset_class)
  if (!curate) {
    return(df)
  } else {
    return(df %>% format_bsdr_data())
  }
}

download_bsdr_data <- function (date, asset_class) {
  dplyr::bind_rows(Map(download_bsdr_data_single, date, asset_class))
}

# Reverse engineering API:
# http://www.bloombergsdr.com/assets/js/search.js
download_bsdr_data_single <- function (date_range, asset_class, currency = NULL,
  notional_range = NULL) {
  # If only one date is specified for date_range, then add one day to this to
  # define the end date (i.e. we want date range starting from the date
  # specified). This will be only allowable behaviour for time being as this
  # will ensure consistency with other data sources and minimise API hammering.
  if (length(date_range) == 1) {
    date_range <- date_range + lubridate::days(0:1)
  }
  # Format dates to format expected by BBG API
  date_range <- format(date_range, "%Y-%m-%dT%H:%M:%OS6Z")
  # Build payload to BBG API. NB: the payload build process in search.js
  # does not assume datetime_low and _high are specified. However, form always
  # has datetime_low pre-populated. Given logic at start of function, we can
  # assume that both are defined in this scope.
  body <- list(
    data_type = "FULL",
    asset_class = asset_class,
    datetime_low = date_range[1],
    datetime_high = date_range[2])
  if (!is.null(currency)) {
    body <- c(body, currency = currency)
  }
  if (!is.null(notional_range)) {
    body <- c(body, notional_low = notional_range[1],
      notional_high = notional_range[2])
  }
  body <- list(Request = list(pubRequest = body))
  response <- httr::POST(url = bsdr_url(), body = body, config = bsef_header(),
    encode = 'json')
  # Convert response's content to JSON from raw
  # Some nested data frames. So flatten these out.
  response <- response$content
  assertthat::assert_that(is.raw(response))
  response <- jsonlite::fromJSON(rawToChar(response))
  # Drill down response to data set that we are interested in
  df <- response$Response$pubResponse$public_recs
  if (!is.null(df)) {
    return(remove_invalid_dfs(df))
  } else {
    df <- dplyr::data_frame()
    return(df)
  }
}

bsdr_url <- function () {
  "http://www.bloombergsdr.com/bas/bsdrweb"
}

# Header isn't specified in search.js. However, added for sake of consistency
# with BSEF spec.
bsdr_header <- function (version = 1.3) {
  httr::add_headers('bas-version' = version)
}

# CFTC requirements on publicly disclosed data fields
# http://www.gpo.gov/fdsys/pkg/CFR-2014-title17-vol2/pdf/CFR-2014-title17-vol2-part43-appA.pdf
# Using http://www.bloombergsdr.com/api
# http://www.bloombergsdr.com/assets/img/BSDR%20API%20for%20Slice%20Files.pdf
#' @importFrom dplyr %>%
format_bsdr_data <- function (data) {
  # Select fields that are presented on the BSDR web search interface
  mutations <- list(
    ~lubridate::fast_strptime(exec_timestamp, "%Y-%m-%dT%H:%M:%OS%z"),
    ~lubridate::fast_strptime(effective_date, "%m/%d/%Y"),
    ~lubridate::fast_strptime(end_date, "%m/%d/%Y"),
    ~as.numeric(price_notation),
    ~as.numeric(additional_price_notation),
    ~as.numeric(notional_currency_amount_1),
    ~as.numeric(option_strike_price),
    ~as.numeric(option_premium))
  mutation_names <- c("exec_timestamp", "effective_date", "end_date",
    "price_notation", "additional_price_notation", "notional_currency_amount_1",
    "option_strike_price", "option_premium")
  data %>% dplyr::mutate_(.dots = stats::setNames(mutations, mutation_names))
}

remove_invalid_dfs <- function (df) {
  replace_with_na <- function (element) {
    if (is.data.frame(element) && ncol(element) == 0) {
      return(rep(NA, nrow(df)))
    } else {
      return(element)
    }
  }
  dplyr::as_data_frame(Map(replace_with_na, df))
}
