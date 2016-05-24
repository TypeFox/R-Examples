if (getRversion() >= "2.15.1")
  utils::globalVariables(c('REGION', 'SETTLEMENTDATE', 'PERIODTYPE'))

#' Download data from AEMO website
#'
#' Trading interval level price and demand data is available on the AEMO
#' website. This function downloads these files for each region and period
#' requested.
#'
#' \href{http://www.aemo.com.au/About-AEMO/Legal-Notices/Copyright-Permissions}{AEMO's copyright permissions}
#' specify that:
#'
#' "AEMO confirms its general permission for anyone to use AEMO Material
#' for any purpose, but only with accurate and appropriate attribution of the
#' relevant AEMO Material and AEMO as its author".
#'
#' @param regions a case invariant character vector of abbreviated AEMO regions.
#' Must be one of \code{"nsw"}, \code{"vic"}, \code{"sa"}, \code{"qld"},
#' \code{"snowy"} or \code{"tas"}. Snowy region data is only available up
#' to and including Jun 2008 while Tasmania region data is available from
#' and including May 2005. Otherwise, all regions have data available from and
#' including Dec 1998 up to the end of the month previous to the current month.
#' @param years a numeric vector of years in the form \code{YYYY}. See the
#' \code{regions} parameter for more info.
#' @param months a numeric vector of months in the form \code{MM}.See the
#' \code{regions} parameter for more info.
#' @param path defaults to \code{"."} (current working directory) and denotes
#' the path to which the AEMO files are downloaded.
#' @return data files whose number corresponds to the maximum lenght of one of
#' the parameters. The file names have the pattern: "[REGION][YYYY][MM].csv".
#' @examples
#' \dontrun{
#' get_aemo_data('nsw', 2014, 1)
#' }
#' @references
#' \href{http://www.aemo.com.au/Electricity/Data/Price-and-Demand/Aggregated-Price-and-Demand-Data-Files}{AEMO website}
#' @importFrom assertthat assert_that
#' @export

get_aemo_data <- function (regions, years, months, path = '.')
{
  regions <- toupper(regions)
  periods <- as.numeric(paste0(years, formatC(months, width = 2, flag = '0')))
  assert_that(all(regions %in% aemo_regions()),
    all(periods %in% aemo_periods()))
  url <- aemo_data_url(regions, years, months)
  destfile <- file.path(path, aemo_data_file_name(regions, years, months))
  invisible(Map(download.file, url = url, destfile = destfile))
}

#' Collate AEMO data
#'
#' This function collates AEMO data files that are on the \code{path}. These
#' are identified by having a filename of the form
#' \code{[REGION][YYYY][MM].csv} (e.g. \code{"NSW201405.csv"}). If no such
#' file exists on \code{path}, all available data files are downloaded from the
#' AEMO website using \code{\link{get_aemo_data}}.
#'
#' @param path a string containing the location of the AEMO (CSV) data files
#' that you wish to collate. These must have filenames of the form specified
#' above.
#' @param remove_files a boolean flag indicating whether the AEMO CSV files
#' should be removed once the data sets are collated.
#' @return a data frame containing five fields: \code{"REGION"} (factor),
#' \code{"SETTLEMENTDATE"} (POSIXct), \code{"TOTALDEMAND"} (numeric),
#' \code{"RRP"} (numeric) and \code{"PERIODTYPE"} (factor).
#' @examples
#' \dontrun{
#' collate_aemo_data()
#' }
#' @importFrom dplyr rbind_all %.% mutate
#' @importFrom lubridate ymd_hms
#' @export

collate_aemo_data <- function (path = '.', remove_files = TRUE)
{
  aemo_files <- list_aemo_data_files(path)
  if (identical(aemo_files, character()))
  {
    info <- all_regions_periods_combos()
    message('Downloading data from AEMO....')
    get_aemo_data(info$region, info$years, info$months, path)
    aemo_files <- list_aemo_data_files(path)
  }
  message('Reading AEMO data files...')
  aemo_dfs <- lapply(aemo_files, read.csv, stringsAsFactors = FALSE)
  if (remove_files) {
    message('Removing AEMO data files...')
    clean_up_aemo_data_files(path)
  }
  message('Collating AEMO data...')
  aemo <- rbind_all(aemo_dfs)
  message('Formatting data frame...')
  aemo %.% mutate(REGION = as.factor(REGION),
    SETTLEMENTDATE = ymd_hms(SETTLEMENTDATE, truncated = 1),
    PERIODTYPE = as.factor(PERIODTYPE))
}

#' Remove AEMO CSV files
#'
#' This function removes AEMO data files that are on the \code{path}. These
#' are identified by having a filename of the form
#' \code{[REGION][YYYY][MM].csv} (e.g. \code{"NSW201405.csv"}).
#'
#' @inheritParams collate_aemo_data
#' @return no return value. However, AEMO files that are identified are removed.
#' @export

clean_up_aemo_data_files <- function (path = '.')
  file.remove(list_aemo_data_files(path))

#' AEMO data set May 2009 - May 2014
#'
#' A sample data set containing trading interval price and demand data sourced
#' from the AEMO website on or about 21 June 2014. This was downloaded and
#' collated using \code{\link{collate_aemo_data}}.
#'
#' @format \code{"REGION"} (factor), \code{"SETTLEMENTDATE"} (POSIXct),
#' \code{"TOTALDEMAND"} (numeric), \code{"RRP"} (numeric) and
#' \code{"PERIODTYPE"} (factor).
#' @references
#' \href{http://www.aemo.com.au/Electricity/Data/Price-and-Demand/Aggregated-Price-and-Demand-Data-Files}{AEMO website}
#' @examples
#' data(aemo)
#' head(aemo)
"aemo"

aemo_regions <- function ()
  c('NSW', 'QLD', 'VIC', 'SA', 'TAS', 'SNOWY')

#' @importFrom lubridate year today month
aemo_periods <- function ()
{
  yy_seq <- 1998:year(today())
  period_min <- 199812
  period_max <- as.numeric(paste0(year(today()),
    formatC(month(today()) - 1, width = 2, flag = '0')))
  all_periods <- outer(yy_seq, formatC(1:12, width = 2, flag = '0'), paste0)
  all_periods <- sort(as.numeric(all_periods))
  in_period <- all_periods >= period_min & all_periods <= period_max
  all_periods[in_period]
}

aemo_data_url_stub <- function (regions, years, months)
{
  stubs <- vector('character', NROW(paste0(regions, years, months)))
  stubs[] <- 'http://www.nemweb.com.au/mms.GRAPHS/data/DATA'
  stubs[years == 1999 & months == 11] <-
    'http://www.nemmco.com.au/mms/data/DATA'
  return (stubs)
}

aemo_data_url <- function (regions, years, months)
{
  paste0(aemo_data_url_stub(regions, years, months), years,
    formatC(months, width = 2, flag = '0'), '_', toupper(regions), '1.csv')
}

aemo_data_file_name <- function (regions, years, months)
  paste0(toupper(regions), years, formatC(months, width = 2, flag = '0'), '.csv')

list_aemo_data_files <- function (path = '.')
  list.files(path, '(NSW|QLD|VIC|SA|TAS|SNOWY)[[:digit:]]{6}\\.csv')

#' @importFrom stringr str_sub
all_regions_periods_combos <- function ()
{
  all_periods <- aemo_periods()
  regions <- rep(aemo_regions(), each = NROW(all_periods))
  years <- as.numeric(rep(str_sub(all_periods, 1, 4), NROW(aemo_regions())))
  months <- as.numeric(rep(str_sub(all_periods, 5), NROW(aemo_regions())))
  is_invalid_tas <- (years * 100 + months < 200505) & regions == 'TAS'
  is_invalid_snowy <- (years * 100 + months > 200806) & regions == 'SNOWY'
  is_invalid_tas_snowy <- is_invalid_tas | is_invalid_snowy
  list(regions = regions[!is_invalid_tas_snowy],
    years = years[!is_invalid_tas_snowy],
    months = months[!is_invalid_tas_snowy])
}
