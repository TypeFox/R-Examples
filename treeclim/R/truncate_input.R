##' Truncate climate and tree data to common interval
##' 
##' Truncate climate and tree data either to common shared interval
##' (default) or to specified range of years. Depending on minmonth,
##' the climate data will need to start one year earlier than the tree
##' data, if data from the previous year should be used.
##' @param chrono a tree-ring chronology
##' @param climate the climate data as returned by as_tcclimate
##' @param timespan the timespan for truncating as vector with min and max year
##' @param minmonth the earliest month used for the calibration, as
##' returned by check_months
##' @param moving moving or not? (logical)
##' @return a list of truncated data.frames for climate and tree data
##' @keywords manip internal
truncate_input <- function(chrono, climate, timespan = NULL, minmonth,
                           moving, silent = FALSE) {

  ## get time spans of both input data
  chrono_years <- as.numeric(row.names(chrono))
  climate_years <- sort(unique(climate[, 1]))

  ## calculate overlap
  if (chrono_years[1] <= climate_years[1]) {
    overlap <- na.omit(climate_years[match(chrono_years, climate_years)])
  } else {
    overlap <- na.omit(chrono_years[match(climate_years, chrono_years)])
  }

  ## take maximum overlap, when timespan is not set
  if (is.null(timespan)) {
    start_year <- overlap[1]
    end_year <- tail(overlap, 1)
  } else {
    if (minmonth > 0) {
      if (!is.element(timespan[1], overlap) ||
          !is.element(timespan[2], overlap)) {
        stop(paste("`timespan` has to be between ", overlap[1],
                   " and ", tail(overlap, 1),
                   " for start dates in current year.",
                   sep = ""))
      } else {
        start_year <- timespan[1]
        end_year <- timespan[2]
      }
    } else {
      if (!is.element(timespan[1], overlap) ||
          !is.element(timespan[2], overlap)) {
        stop(paste("`timespan` has to be between ", overlap[1] +
                   1, " and ", tail(overlap, 1),
                   " for start dates in previous year.", sep = ""))
      } else {
        start_year <- timespan[1]
        end_year <- timespan[2]
      }
    }
  }
  
  ## check if a previous year is available in climatic data; otherwise
  ## set start_year + 1
  if (minmonth < 0 && !((start_year - 1) %in% climate_years)) { 
    offset <- 1
    if (is.null(timespan) & !silent) {
      message(paste("treeclim tries to use the maximum overlap in timespan for chronology and climate data. The overlap starts in ",
                    start_year,
                    ", but to be able to use climate data from the previous year (as you chose by setting 'selection' accordingly), the analysis starts in ",
                    start_year + 1, ".", sep = ""))
    } else {
        if (!silent) {
            message(paste("`start_year` is set from", start_year, "to",
                          start_year + 1,
                          "to be able to use climate data from the previous year."))
        }
    }
    
  } else {
    offset <- 0
  }

  ## make sure that data get truncated properly	
  if (minmonth < 0) { 
    interval_climate <- (start_year - 1 + offset):end_year
    interval_chrono <- (start_year + offset):end_year
    pad <- FALSE
  } else {
    interval_climate <- (start_year + offset):end_year
    interval_chrono <- (start_year + offset):end_year
    pad <- TRUE
  }

  a <- chrono_years %in% interval_chrono
  b <- climate[, 1] %in% interval_climate
  
  ## report time span used
  ## calculate timespan for analysis for reporting
  if (!moving & !silent) {
    run_years <- chrono_years[a]
    cat("Running for timespan ", run_years[1], " - ", tail(run_years, 1), "...\n",
        sep = "")
  }

  ## finally truncate data
  chrono_trunc <- chrono[a, 1]
  climate_trunc <- climate[b, ]

  list(chrono = chrono_trunc,
       climate = climate_trunc,
       pad = pad)
}
