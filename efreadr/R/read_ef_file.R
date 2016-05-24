#' Load a CSV file from the European Flux Database
#'
#' File name is parsed to extract year, site identification and aggregation type.
#' The file name must point to a valid European Fluxes file, in CSV format and must
#' resolve to a valid file format name.
#'
#' Year, file name and site identification are added as fields in the returned data frame
#' as 'efreadr_year', 'efreader_file_name' and 'efreader_site_id'.
#'
#' @note For semi-hourly L4 aggregation (i.e. "h" aggregation in file name) the last row is
#' reported as month 1, day 1, hour 00:00. A normal date conversion would convert this date to
#' be the very first half-hour in January 1st of the current year whereas it should be the first half-hour of the
#' January 1st of the following year.
#' Therefore a class date field ('efreader_date') is added to the returned data frame holding the correct
#' date (ie: 1st January of the following year).
#'
#' @param file_name Full path to 1 fluxes file
#' @importFrom readr read_csv
#' @importFrom ensurer ensure_that
#' @importFrom magrittr %>%
#' @export
#' @examples
#' file_name <- system.file("extdata", "CEIP_EC_L4_d_FABar_2015_v02.txt", package = "efreadr")
#' read_ef_file(file_name)
#' @return a data frame as loaded from the file, added with 'efreadr_year', 'efreadr_file_name' and 'efreadr_site_id' columns, and 'efreadr_date' column for half-hourly fluxes
read_ef_file <- function(file_name) `: dataframe_with_filename_and_siteid` ({
  file_name %>% length() %>% ensure_that(. == 1, err_desc = "Trying to load too many files at once or none at all, try with one at a time...")

  message(sprintf("Loading file '%s'...", file_name))

  file_name %>% ensure_that(file.exists(.), err_desc = "File does not exist.")

  flux_data     <- readr::read_csv(file_name)
  file_metadata <- c(
    sapply(
      strsplit(
        basename(file_name),
        "_",
        fixed = TRUE),
      `[`,
      c(4, 5, 6)))

  file_metadata %>% length() %>% ensure_that(. == 3, err_desc = "Fluxes file name malformed; is it really a fluxes file from European Fluxes Database?")
  names(file_metadata) <- c("aggregation", "site_code", "year")

  if (file_metadata["aggregation"] == "h") {
    c("Month", "Day") %in% colnames(flux_data) %>% sum() %>% ensure_that(. == 2, err_desc = "'Month' and/or 'Day' columns are missing from fluxes file; is it really a fluxes file from European Fluxes Database?")

    flux_data$efreadr_date <- as.Date(paste(file_metadata["year"], flux_data$Month, flux_data$Day, sep = "-"))
    flux_data$efreadr_date[nrow(flux_data)] <- as.Date(
      paste(
        as.numeric(file_metadata["year"]) + 1,
        01,
        01,
        sep = "-"))
  }

  flux_data$efreadr_year      <- file_metadata["year"]
  flux_data$efreadr_site_id   <- paste(substr(file_metadata["site_code"], 1, 2), substr(file_metadata["site_code"], 3, 5), sep = "-")
  flux_data$efreadr_file_name <- file_name

  message(sprintf("Imported flux data for site '%s', year %s", flux_data$efreadr_site_id[1], file_metadata["year"]))

  flux_data
})
