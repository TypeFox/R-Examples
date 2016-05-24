#' Load all European Fluxes CSV files in one or more directories, bind all observations and fill missing data with NA_real
#'
#' European fluxes CSV files are distributed as one or more zip-compressed files from http://gaia.agraria.unitus.it .
#' Once unzipped, all CSV files are to be found in uniquely identifed directories.
#'
#' All CSV files in that or those directories will be loaded and returned as a single row-wise bound data frame.
#' The function assumes the file name pattern is like "^CEIP_EC_Ln_a_[a-zA-Z0-9]{5}_20[0-9]{2}_v[0-9]{2}\\.txt$"
#' where n in level and a is aggregation period (optionally given as function arguments)
#' Year, file name and site identification are added as fields in the returned data frame
#' as 'efreadr_year', 'efreader_file_name' and 'efreader_site_id'.
#'
#' Optionally (but default behaviour), unavailable measures are converted to 'NA_real' values.
#'
#' @note All files in the same directory should belong to the same aggregation/level
#' combination in order to row-wise build a consistent dataframe. Note that when 'level_l'
#' and/or 'aggregation' arguments are not given, all files will be loaded for in the directory,
#' regardless of their level/aggregation.
#'
#' @note For semi-hourly L4 aggregation (i.e. "h" aggregation in file name) the last row is
#' reported as month 1, day 1, hour 00:00. A normal date conversion would convert this date to
#' be the very first half-hour of the current year whereas it should be the very first half-hour of the
#' following year. Therefore a class date field ('efreader_date') is added to the returned data frame holding the correct
#' date (ie: January 1st of the following year).
#'
#' @param dirs    a vector of directories where fluxes files are looked for. Defaults to current directory.
#' @param level_l level of fluxes files (defaults to NULL). Allowed levels are (currently) 3 and 4. When NULL, either L3 and L4 files are looked for.
#' @param aggregation aggregation of data (defaults to NULL) Allowed aggregations are (currently) "h" (half-hourly) and "d" (daily). When NULL, either "d" and "h" files are looked for.
#' @param fill_value a code for a not available (NA) observation in CSV file. When the argument is given (default behaviour) the observations with 'fill_value' values are converted to NA_real
#'
#' @importFrom ensurer  ensure_that
#' @importFrom ensurer  check_that
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom dplyr    bind_rows
#' @importFrom dplyr    mutate_each
#' @examples
#' dir_name <- system.file("extdata", package = "efreadr")
#' read_ef_files(dir_name)
#' @export
#' @return a data frame as loaded from the file, added with 'efreadr_year', 'efreadr_file_name' and 'efreadr_site_id' columns, and 'efreadr_date' column for half-hourly fluxes
read_ef_files <- function(dirs = getwd(), level_l = NULL, aggregation = NULL, fill_value = "-9999") `: dataframe_with_filename_and_siteid` ({

  allowed_levels <- c(3, 4)
  allowed_aggr   <- c("h", "d")

  dirs %>% dir.exists() %>% sum() %>% ensure_that(. > 0, err_desc = "At least one directory must exist!")

  if (is.null(level_l)) {
    level_l <- allowed_levels
    message(sprintf("Using %s as level", paste(level_l, collapse=",")))
  }

  if (is.null(aggregation)) {
    aggregation <- allowed_aggr
    message(sprintf("Using %s as aggregation", paste(aggregation, collapse = ",")))
  }

  level_l     %in% allowed_levels %>% sum() %>% ensure_that(. > 0, err_desc = "Level not allowed")
  aggregation %in% allowed_aggr   %>% sum() %>% ensure_that(. > 0, err_desc = "Aggregation not allowed")

  level_l = paste0("[", paste(allowed_levels, collapse = ""), "]")
  aggregation <- paste0("[", paste(allowed_aggr, collapse = ""), "]")

  file_pattern <- paste0(
    "^CEIP_EC_L",
    level_l,
    "_",
    aggregation,
    "_[a-zA-Z0-9]{5}_20[0-9]{2}_v[0-9]{2}\\.txt$")

  file_list <- unlist(lapply(
    dirs,
    list.files,
    pattern = file_pattern,
    full.names = TRUE))
print(dirs)
  file_list %>% length(.) %>% check_that(. > 0, err_desc = "Trying to load too many files at once or none at all, try with one at a time...")

  fluxes <- dplyr::bind_rows(
    lapply(
      file_list,
      read_ef_file))

  if (!is.null(fill_value)) {
    fluxes %<>%
      mutate_each(
        funs(. = ifelse(. == fill_value, NA_real_, .)),
        -starts_with("efreadr_"))
  }

  fluxes
})
