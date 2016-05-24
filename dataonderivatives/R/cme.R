if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("CME_ASSET_CLASSES_LONG", "CME_ASSET_CLASSES_SHORT"))
}

CME_ASSET_CLASSES_LONG <-
  c(CO = "commodities", CR = "credit", FX = "fx", IR = "rates")
CME_ASSET_CLASSES_SHORT <-
  c(CO = "commodity", CR = "CDS", FX = "fx", IR = "irs")

#' Get CME SDR data
#'
#' The CME Swap Data Repository (SDR) is a registered U.S. swap data repository
#' that allows market participants to fulfil their public disclosure
#' obligations under U.S. legislation. CME is required to make publicly
#' available price, trading volume and other trading data. It publishes this
#' data on an FTP site.
#'
#' @param date the date for which data is required as Date or DateTime object.
#'   It will only use the year, month and day elements to
#'   determine the set of trades to return. It will return the set of trades
#'   for the day starting on \code{date}.
#' @param asset_class the asset class for which you would like to download
#'   trade data. Valid inputs are \code{"CR"} (credit), \code{"IR"} (rates),
#'   \code{"FX"} (foreign exchange), \code{"CO"} (commodities). Can be a
#'   vector of these. Defaults to \code{NULL} which corresponds to all asset
#'   classes.
#' @param curate a logical flag indicating whether raw data should be returned
#'   or whether the raw data should be processed (default). The latter involves
#'   selecting particular fields and formatting these as seemed appropriate
#'   based on data reviews at the time the formatting was coded.
#' @return a data frame containing the requested data, or an empty data frame
#'   if data is unavailable
#' @importFrom dplyr %>%
#' @references \href{http://www.cmegroup.com/trading/global-repository-services/cme-swap-data-repository.html}{CME SDR}
#' @examples
#' \dontrun{
#' library (lubridate)
#' # All asset classes for day starting 6 May 2015
#' get_cme_data(ymd(20150506))
#' # Only IR and FX asset classes
#' get_cme_data(ymd(20150506), c("IR", "FX"))
#' }
#' @export


get_cme_data <- function (date, asset_class = NULL, curate = TRUE) {
  valid_asset_classes <- c('CR', 'FX', 'IR', 'CO')
  if (is.null(asset_class)) {
    asset_class <- valid_asset_classes
  }
  assertthat::assert_that(lubridate::is.instant(date), length(date) == 1,
    all(asset_class %in% valid_asset_classes))
  Map(download_cme_zip, date, asset_class)
  on.exit(clean_cme_files())
  dplyr::bind_rows(Map(read_cme_file, date, asset_class, curate))
}

download_cme_zip <- function (date, asset_class) {
  ftp_url <- cme_ftp_url(date, asset_class)
  tmpdir <- file.path(tempdir(), "cme")
  if (!dir.exists(tmpdir)) dir.create(tmpdir, recursive = TRUE)
  tmpfile_pattern <- cme_file_name(date, asset_class)
  tmpfile <- tempfile(tmpfile_pattern, tmpdir, fileext = ".zip")
  files <- character()
  try({
    suppressWarnings({
      # On Windows only (?) this downloads a file with non-zero byte size even if
      # no data (and file!) exists on the date / asset class requested. So work
      # around this.
      downloader::download(url = ftp_url, destfile = tmpfile, quiet = TRUE)
      # res is an em  pty character vector if no data exists for that date and asset
      # class, even though tmpfile is a non-zero sized zip file (hence unlinked).
      # Mimic download.file return value (with 0 = success). Also need to
      files <- utils::unzip(tmpfile, exdir = tmpdir)
    })
  }, silent = TRUE)
  if (identical(files, character())) res <- -1 else res <- 0
  unlink(tmpfile)
  invisible(res)
}

#' @importFrom dplyr %>%
read_cme_file <- function (date, asset_class, curate) {
  tmpdir <- file.path(tempdir(), 'cme/')
  csvfile <- list.files(tmpdir, cme_file_name(date, asset_class),
    full.names = TRUE)
  if (length(csvfile) < 1) {
    return(dplyr::data_frame())
  } else {
    # Should only have one file per day. Use first if multiple matches
    if (!curate) {
      return(readr::read_csv(csvfile[1]))
    } else {
      col_types <- specify_cme_col_types(asset_class)
      df <- readr::read_csv(csvfile[1], col_types = col_types)
      if (identical(asset_class, "IR")) {
        # IR file doesn't use title case for this field.
        df <- df %>% dplyr::rename_("Block/Off Facility" = "Block/Off facility")
      }
      return(df)
    }
  }
}

specify_cme_col_types <- function (asset_class) {
  if (identical(asset_class, "IR")) {
    return(list(
      `Event` = readr::col_character(),
      `Execution Timestamp` = readr::col_datetime("%m/%d/%Y %H:%M:%S"),
      `Dissemination Time` = readr::col_datetime("%m/%d/%Y %H:%M:%S"),
      `Cleared` = readr::col_character(),
      `Collateralization` = readr::col_character(),
      `End-User Exception` = readr::col_character(),
      `Bespoke` = readr::col_character(),
      `Block/Off facility` = readr::col_character(),
      `Execution Venue` = readr::col_character(),
      `UPI` = readr::col_character(),
      `Product` = readr::col_character(),
      `Contract Type` = readr::col_character(),
      `Effective Date` = readr::col_date("%m/%d/%Y"),
      `Maturity Date` = readr::col_date("%m/%d/%Y"),
      `Upfront Payment` = readr::col_double(),
      `Upfront Payment Currency` = readr::col_character(),
      `Upfront Payment Date` = readr::col_date("%m/%d/%Y"),
      `Settlement Currency` = readr::col_character(),
      `Leg 1 Type` = readr::col_character(),
      `Leg 1 Fixed Rate` = readr::col_double(),
      `Leg 1 Floating Index` = readr::col_character(),
      `Leg 1 Designated Maturity` = readr::col_character(),
      `Leg 1 Spread` = readr::col_double(),
      `Leg 1 Day Count Convention` = readr::col_character(),
      `Leg 1 Notional` = readr::col_double(),
      `Leg 1 Notional Currency` = readr::col_character(),
      `Leg 1 Payment Frequency` = readr::col_character(),
      `Leg1 Reset Frequency` = readr::col_character(),
      `Leg 2 Type` = readr::col_character(),
      `Leg 2 Fixed Rate` = readr::col_double(),
      `Leg 2 Floating Index` = readr::col_character(),
      `Leg 2 Designated Maturity` = readr::col_character(),
      `Leg 2 Spread` = readr::col_double(),
      `Leg 2 Day Count Convention` = readr::col_character(),
      `Leg 2 Notional` = readr::col_double(),
      `Leg 2 Notional Currency` = readr::col_character(),
      `Leg 2 Payment Frequency` = readr::col_character(),
      `Leg 2 Reset Frequency` = readr::col_character(),
      `Embedded Option` = readr::col_character(),
      `Option Strike Price` = readr::col_double(),
      `Option Type` = readr::col_character(),
      `Option Family` = readr::col_character(),
      `Option Currency` = readr::col_character(),
      `Option Premium` = readr::col_double(),
      `Option Lockout Period` = readr::col_character(),
      `Option Expiration Date` = readr::col_date("%m/%d/%Y"),
      `Asset Class` = readr::col_character(),
      `Rpt ID` = readr::col_character(),
      `Prev Rpt ID` = readr::col_character()))
  } else if (identical(asset_class, "FX")) {
    return (list(
      `Event` = readr::col_character(),
      `Execution Timestamp` = readr::col_datetime("%m/%d/%Y %H:%M:%S"),
      `Dissemination Time` = readr::col_datetime("%m/%d/%Y %H:%M:%S"),
      `Cleared` = readr::col_character(),
      `Collateralization` = readr::col_character(),
      `End-User Exception` = readr::col_character(),
      `Bespoke` = readr::col_character(),
      `Block/Off Facility` = readr::col_character(),
      `Execution Venue` = readr::col_character(),
      `UPI` = readr::col_character(),
      `Product` = readr::col_character(),
      `Asset Class` = readr::col_character(),
      `Contract Type` = readr::col_character(),
      `Effective Date` = readr::col_date("%m/%d/%Y"),
      `Maturity Date` = readr::col_date("%m/%d/%Y"),
      `Currency 1` = readr::col_character(),
      `Currency 2` = readr::col_character(),
      `Notional Amount 1` = readr::col_double(),
      `Notional Amount 2` = readr::col_double(),
      `Exchange Rate` = readr::col_double(),
      `Settlement Date` = readr::col_date("%m/%d/%Y"),
      `Settlement Currency` = readr::col_character(),
      `Embedded Option` = readr::col_character(),
      `Option Strike Price` = readr::col_double(),
      `Option Type` = readr::col_character(),
      `Option Family` = readr::col_character(),
      `Option Currency` = readr::col_character(),
      `Option Premium` = readr::col_double(),
      `Option Lockout Period` = readr::col_character(),
      `Option Expiration Date` = readr::col_date("%m/%d/%Y"),
      `Rpt ID` = readr::col_character(),
      `Prev Rpt ID` = readr::col_character()))
  } else if (identical(asset_class, "CO")) {
    return(list(
      `Event` = readr::col_character(),
      `Execution Timestamp` = readr::col_datetime("%m/%d/%Y %H:%M:%S"),
      `Dissemination Time` = readr::col_datetime("%m/%d/%Y %H:%M:%S"),
      `Cleared` = readr::col_character(),
      `Collateralization` = readr::col_character(),
      `End-User Exception` = readr::col_character(),
      `Bespoke` = readr::col_character(),
      `Block/Off Facility` = readr::col_character(),
      `Execution Venue` = readr::col_character(),
      `UPI` = readr::col_character(),
      `Product` = readr::col_character(),
      `Asset Class` = readr::col_character(),
      `Contract Type` = readr::col_character(),
      `Effective Date` = readr::col_date("%m/%d/%Y"),
      `Maturity Date` = readr::col_date("%m/%d/%Y"),
      `Settlement Currency` = readr::col_character(),
      `Leg 1 Type` = readr::col_character(),
      `Leg 1 Fixed Payment` = readr::col_double(),
      `Leg 1 Index` = readr::col_character(),
      `Leg 1 Index Location` = readr::col_character(),
      `Leg 1 Spread` = readr::col_double(),
      `Leg 1 Averaging Method` = readr::col_character(),
      `Leg 1 Delivery Point` = readr::col_character(),
      `Leg 2 Type` = readr::col_character(),
      `Leg 2 Fixed Payment` = readr::col_double(),
      `Leg 2 Index` = readr::col_character(),
      `Leg 2 Index Location` = readr::col_character(),
      `Leg 2 Spread` = readr::col_double(),
      `Leg 2 Averaging Method` = readr::col_character(),
      `Leg 2 Delivery Point` = readr::col_character(),
      `Embedded Option` = readr::col_character(),
      `Option Strike Price` = readr::col_double(),
      `Option Type` = readr::col_character(),
      `Option Family` = readr::col_character(),
      `Option Currency` = readr::col_character(),
      `Option Premium` = readr::col_double(),
      `Option Lockout Period` = readr::col_character(),
      `Option Expiration Date` = readr::col_date("%m/%d/%Y"),
      `Rpt ID` = readr::col_character(),
      `Prev Rpt ID` = readr::col_character(),
      `Price Unit` = readr::col_character()))
  } else {
    # Must be Credit class. And I can't find a sample file for this. So
    # no col type spec is provided
    return (NULL)
  }
}

clean_cme_files <- function () {
  unlink(file.path(tempdir(), 'cme'), recursive = TRUE)
}

cme_ftp_url <- function (date, asset_class) {
  # Eg URLS
  # ftp://ftp.cmegroup.com/sdr/fx/2015/03/RT.FX.20150301.csv.zip
  # ftp://ftp.cmegroup.com/sdr/rates/2013/07/RT.IRS.20130702.csv.zip
  # ftp://ftp.cmegroup.com/sdr/commodities/2015/02/RT.COMMODITY.20150201.csv.zip
  asset_class_long <- CME_ASSET_CLASSES_LONG[asset_class]
  paste0("ftp://ftp.cmegroup.com/sdr/",
    tolower(asset_class_long), "/",
    lubridate::year(date), "/",
    formatC(lubridate::month(date), width = 2, flag = "0"), "/",
    cme_file_name(date, asset_class), ".zip")
}

cme_file_name <- function (date, asset_class) {
  paste0("RT.", toupper(CME_ASSET_CLASSES_SHORT[asset_class]),
    format(date, ".%Y%m%d.csv"))
}
