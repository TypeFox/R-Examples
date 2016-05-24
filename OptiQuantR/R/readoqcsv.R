#'Read OptiQuant's exported .csv file.
#'
#'\code{readoqcsv} returns a keyed data.table and data.frame ready for fast and
#'easy analysis and reporting.
#'
#'Optimized for speed and preparation of data for custom analysis through simple
#'user input. Reads in .csv eliminating unnecessary columns. Sets column classes
#'as appropriate for subsequent use in analysis and reporting. Based on the data
#'imported, it creates new columns to facilitate simple user input for
#'customizing the analysis. Also creates column to facilitate the simple
#'inclusion of a margin of error when filtering data.
#'
#'readoqcsv("/Users/jpinelo/JPLab/R Projects/opti_557d819b88209.csv")
#'Note that full path is necessary when the file is not in the current working
#'directory.
#'
#'readoqcsv("opti_557d819b88209.csv")
#'Note that when file is in current working directory, the name and extension of
#'the file are enough.
#'
#'
#'@param
#'x    Single string. Name / path of the file to read in. See details.
#'
#'
#'@examples
#'readoqcsv(system.file("extdata", "testdata.csv", package = "OptiQuantR"))
#'
#'@return
#' Classes ‘data.table’ and 'data.frame':	n obs. of  17 variables:
#'
#' Where n is the number of rows no the csv file - 1.
#'
#' Names of variables (columns): "timeStampO", "session_id", "session_started",
#' "session_finished", "session_total_counted", "gate_id", "person_kind", "id",
#' "year", "month", "day", "weekday", "hour", "week", "min_start_off",
#' "session_length", "timeStamp"
#'@export
readoqcsv <- function(x) {
# Load data --------------------------------------------------------------------
  # colclasses and comment.char improve speed.
  # date cols as character for coercion as posix date
      dataImport <- utils::read.table(file = x,
                                      header = TRUE,
                                      sep = "," ,
                                      dec = "." ,
                                      colClasses = c("character",
                                                     "integer",
                                                     "character",
                                                     "character",
                                                     "integer",
                                                     "NULL",
                                                     "NULL",
                                                     "NULL",
                                                     "integer",
                                                     "NULL",
                                                     "NULL",
                                                     "NULL",
                                                     "NULL",
                                                     "factor",
                                                     "NULL",
                                                     "NULL",
                                                     "NULL",
                                                     "NULL",
                                                     "NULL"),
                                      comment.char = "")
  data.table::setDT(dataImport)

# Prepare data -----------------------------------------------------------------
  # Rename col count_datetime to time_stampO to save name time)tamp for analysis
  data.table::setnames(dataImport, "count_datetime", "timeStampO")

  # POSIXlt is needed for analysis stage
  dataImport$timeStampO <- strptime(dataImport$timeStampO,
                                    "%m/%d/%y  %H:%M")

  dataImport$session_started <- strptime(dataImport$session_started,
                                         "%m/%d/%y  %H:%M")

  dataImport$session_finished <- strptime(dataImport$session_finished,
                                          "%m/%d/%y  %H:%M")

  # Create id for each event for analysis stage
  dataImport$id <- c(1:nrow(dataImport))
#  data.table::setkey(dataImport, id)

  # copy data to coerce dates to POSIXct
  dataImportCt <- dataImport

  # coerce all time cols to POSIXct
  dataImportCt$timeStampO <- as.POSIXct(dataImport$timeStampO,
                                        format = "%m/%d/%y  %H:%M")

  dataImportCt$session_started <- as.POSIXct(dataImport$session_started,
                                             format = "%m/%d/%y  %H:%M")

  dataImportCt$session_finished <- as.POSIXct(dataImport$session_finished,
                                              format = "%m/%d/%y  %H:%M")

  # extract date time elements for easier use
  dataImportCt$year <- as.integer(lubridate::year(dataImportCt$timeStampO))

  dataImportCt$month <- as.integer(lubridate::month(dataImportCt$timeStamp))

  dataImportCt$day <- as.integer(lubridate::day(dataImportCt$timeStamp))

  dataImportCt$weekday <- as.factor(weekdays(dataImportCt$timeStampO,
                                             abbreviate = TRUE))

  dataImportCt$hour <- as.integer(lubridate::hour(dataImportCt$timeStampO))

  dataImportCt$week <- as.integer(lubridate::week(dataImportCt$timeStampO))

  # extract mins from session_started to identify sessions
  # which started off late or early
  dataImport$min_start_off <- lubridate::minute(dataImport$session_started)

  dataImportCt$min_start_off <- lubridate::minute(dataImportCt$session_started)

  # create col session_length = session_finished - session_started
   dataImport$session_length <- with(dataImport, difftime(session_finished,
                                                          session_started,
                                                          unit = "mins"))

  dataImportCt$session_length <- with(dataImportCt, difftime(session_finished,
                                                             session_started,
                                                             unit = "mins"))

  # insert data filtering here
  dataImportCt$timeStamp <- dataImportCt$timeStampO

    # Keep dataImport and dataImportCt for further use in later releases
  events <- dataImportCt
# Dataset to be used by user ---------------------------------------------------
  return(events)
}
