if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("DDR_ASSET_CLASSES"))
}

DDR_ASSET_CLASSES <- c("CR" = "CREDITS", 'EQ' = "EQUITIES", 'FX' = "FOREX",
  'IR' = "RATES", 'CO' = "COMMODITIES")

ddr_file_name <- function (date, asset_class) {
  paste0(paste("CUMULATIVE", DDR_ASSET_CLASSES[asset_class], format(date, "%Y"),
    format(date, "%m"), format(date, "%d"), sep="_"))
}

ddr_url <- function (date, asset_class) {
  assertthat::assert_that(assertthat::has_name(DDR_ASSET_CLASSES, asset_class))
  stump <- "https://kgc0418-tdw-data-0.s3.amazonaws.com/slices/"
  # "https://kgc0418-tdw-data-0.s3.amazonaws.com/slices/CUMULATIVE_CREDITS_2015_04_29.zip"
  paste0(stump, ddr_file_name(date, asset_class), ".zip")
}

download_ddr_zip <- function (date, asset_class) {
  zip_url <- ddr_url(date, asset_class)
  tmpfile_pattern <- ddr_file_name(date, asset_class)
  tmpdir <- file.path(tempdir(), "ddr")
  if (!dir.exists(tmpdir)) dir.create(tmpdir, recursive = TRUE)
  tmpfile <- tempfile(tmpfile_pattern, tmpdir, fileext = ".zip")
  # Make sure that URL is ok. Is possible to download a "valid" ZIP file which
  # is empty (e.g. for a future date). This is avoided by checking URL status
  # before downloading the file.
  tmpdir <- file.path(tmpdir, date, asset_class)
  if (httr::url_ok(zip_url)) {
    downloader::download(url = zip_url, destfile = tmpfile, quiet = TRUE)
    # Create date/asset_class dir as CSV file name in zip does not reflect date.
    # This makes it harder to ensure read_ddr_file picks up the right file.
    utils::unzip(tmpfile, exdir = tmpdir)
    unlink(tmpfile)
    invisible(0)
  } else {
    # If no data exists, read_ddr_zip will return an empty data_frame()
    invisible(-1)
  }
}

specify_ddr_col_types <- function () {
  list(
    DISSEMINATION_ID = readr::col_integer(),
    ORIGINAL_DISSEMINATION_ID = readr::col_integer(),
    ACTION = readr::col_character(),
    EXECUTION_TIMESTAMP = readr::col_datetime(),
    CLEARED = readr::col_character(),
    INDICATION_OF_COLLATERALIZATION = readr::col_character(),
    INDICATION_OF_END_USER_EXCEPTION = readr::col_character(),
    INDICATION_OF_OTHER_PRICE_AFFECTING_TERM = readr::col_character(),
    "BLOCK_TRADES_AND_LARGE_NOTIONAL_OFF-FACILITY_SWAPS" = readr::col_character(),
    EXECUTION_VENUE = readr::col_character(),
    EFFECTIVE_DATE = readr::col_date(),
    END_DATE = readr::col_date(),
    DAY_COUNT_CONVENTION = readr::col_character(),
    SETTLEMENT_CURRENCY = readr::col_character(),
    ASSET_CLASS = readr::col_character(),
    "SUB-ASSET_CLASS_FOR_OTHER_COMMODITY" = readr::col_character(),
    TAXONOMY = readr::col_character(),
    PRICE_FORMING_CONTINUATION_DATA = readr::col_character(),
    UNDERLYING_ASSET_1 = readr::col_character(),
    UNDERLYING_ASSET_2 = readr::col_character(),
    PRICE_NOTATION_TYPE = readr::col_character(),
    PRICE_NOTATION = readr::col_numeric(),
    ADDITIONAL_PRICE_NOTATION_TYPE = readr::col_character(),
    ADDITIONAL_PRICE_NOTATION = readr::col_numeric(),
    NOTIONAL_CURRENCY_1 = readr::col_character(),
    NOTIONAL_CURRENCY_2 = readr::col_character(),
    ROUNDED_NOTIONAL_AMOUNT_1 = readr::col_numeric(),
    ROUNDED_NOTIONAL_AMOUNT_2 = readr::col_numeric(),
    PAYMENT_FREQUENCY_1 = readr::col_character(),
    PAYMENT_FREQUENCY_2 = readr::col_character(),
    RESET_FREQUENCY_1 = readr::col_character(),
    RESET_FREQUENCY_2 = readr::col_character(),
    EMBEDED_OPTION = readr::col_character(),
    OPTION_STRIKE_PRICE = readr::col_numeric(),
    OPTION_TYPE = readr::col_character(),
    OPTION_FAMILY = readr::col_character(),
    OPTION_CURRENCY = readr::col_character(),
    OPTION_PREMIUM = readr::col_numeric(),
    OPTION_LOCK_PERIOD = readr::col_date(),
    OPTION_EXPIRATION_DATE = readr::col_date(),
    PRICE_NOTATION2_TYPE = readr::col_character(),
    PRICE_NOTATION2 = readr::col_numeric(),
    PRICE_NOTATION3_TYPE = readr::col_character(),
    PRICE_NOTATION3 = readr::col_numeric())
}

read_ddr_file <- function (date, asset_class, curate) {
  tmpdir <- file.path(tempdir(), 'ddr/', date, "/", asset_class, '/')
  ddrfile <- list.files(tmpdir, DDR_ASSET_CLASSES[asset_class], full.names = TRUE)
  if (length(ddrfile) < 1L) {
    return(dplyr::data_frame())
  } else {
    # Should only have one file per day. Use first if multiple matches
    # Use col_types() to specify each colum type as the first 100 rows may
    # not contain valid values. Reviewing col names for different asset classes
    # at 30 Apr 2014 indicates they all have same names. So specify col types
    # explicitly
    if (!curate) {
      return(readr::read_csv(ddrfile[1]), col_types = NULL)
    } else {
      return(readr::read_csv(ddrfile[1], col_types = specify_ddr_col_types()))
    }
  }
}

clean_ddr_files <- function () {
  unlink(file.path(tempdir(), 'ddr'), recursive = TRUE)
}


#' Get DDR data
#'
#' The DTCC Data Repository is a registered U.S. swap data repository that
#' allows market participants to fulfil their public disclosure obligations
#' under U.S. legislation. This function will give you the ability to download
#' trade-level data that is reported by market participants. The field names are
#' (and is assumed to be) the same for each asset class.
#'
#' @param date the date for which data is required as Date or DateTime object.
#'   Only the year, month and day elements of the object are used and it must of
#'   be length one.
#' @param asset_class the asset class for which you would like to download trade
#'   data. Valid inputs are \code{"CR"} (credit), \code{"IR"} (rates),
#'   \code{"EQ"} (equities), \code{"FX"} (foreign exchange), \code{"CO"}
#'   (commodities). Can be a vector of these. Defaults to \code{NULL} which
#'   corresponds to all asset classes.
#' @param curate a logical flag indicating whether raw data should be returned
#'   or whether the raw data should be processed (default). The latter involves
#'   selecting particular fields and formatting these as seemed appropriate
#'   based on data reviews at the time the formatting was coded.
#' @return a \code{tbl_df} that contains the requested data. If no data exists
#'   on that date, an empty data frame (zero columns and rows) is returned.
#' @examples
#' \dontrun{
#' library("lubridate")
#' get_ddr_data(ymd(20140430), "IR")
#' }
#' @references \href{https://rtdata.dtcc.com/gtr/}{DDR Real Time Dissemination
#' Platform}
#' @export

get_ddr_data <- function (date, asset_class = NULL, curate = TRUE) {
  valid_asset_classes <- c('CR', 'EQ', 'FX', 'IR', 'CO')
  if (is.null(asset_class)) {
    asset_class <- valid_asset_classes
  }
  assertthat::assert_that(lubridate::is.instant(date), length(date) == 1,
    all(asset_class %in% valid_asset_classes))
  Map(download_ddr_zip, date, asset_class)
  on.exit(clean_ddr_files())
  dplyr::bind_rows(Map(read_ddr_file, date, asset_class, curate))
}
