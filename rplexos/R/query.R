# Get one table from a SQLite database
get_table <- function(filename, table) {
  # Open connection
  thesql <- src_sqlite(filename, create = FALSE)
  
  if (table %in% src_tbls(thesql)) {
    out <- tbl(thesql, table) %>% collect
  } else {
    out <- data.frame()
  }
  
  # Close connection
  DBI::dbDisconnect(thesql$con)
  
  # Return result
  out
}

# Get one table from a SQLite database
get_list_tables <- function(filename) {
  # Open connection
  thesql <- src_sqlite(filename, create = FALSE)
  
  # Read names
  out <- data.frame(tables = src_tbls(thesql))
  
  # Close connection
  DBI::dbDisconnect(thesql$con)
  
  # Return result
  out
}

#' Get a query for a SQLite file
#'
#' Used internally by \code{\link{query_sql}}. Use that function instead to access data. The
#' use of this function is not recommended
#'
#' @param filename SQLite file location
#' @param sql String containing the SQL query to be performed
#'
#' @seealso \code{\link{query_sql}} to perform standard queries of data
#' 
#' @keywords internal
#' @export
get_query <- function(filename, sql) {
  out <- data.frame()
  thesql <- src_sqlite(filename, create = FALSE)
  try(out <- RSQLite::dbGetQuery(thesql$con, sql))
  DBI::dbDisconnect(thesql$con)
  out
}

# Get a table for all scenarios
get_table_scenario <- function(db, from) {
  # Check inputs
  check_rplexos(db)
  stopifnot(is.character(from), length(from) == 1L)
  
  db %>%
    group_by(scenario, position, filename) %>%
    do(get_table(.$filename, from)) %>%
    ungroup()
}

#' Get query for all scenarios
#'
#' Send a SQL query to all the files in a PLEXOS database object.
#'
#' @inheritParams query_master
#' @param sql String containing the SQL query to be performed
#'
#' @seealso \code{\link{plexos_open}} to create the PLEXOS database object
#' @seealso \code{\link{query_master}} to perform standard queries of data
#' @family special queries
#' 
#' @examples
#' # Process the folder with the solution file provided by rplexos
#' location <- location_solution_rplexos()
#' process_folder(location)
#' 
#' # Query data
#' db <- plexos_open(location)
#' query_sql(db, "SELECT * from day")
#' query_sql(db, "SELECT * from time")
#' 
#' @export
query_sql <- function(db, sql) {
  # Check inputs
  check_rplexos(db)
  stopifnot(is.character(sql), length(sql) == 1L)
  
  # Make sure that columns are reported
  db <- db %>%
    group_by(scenario, position, filename)
  
  # Get query data
  if (!is_parallel_rplexos()) {
    out <- db %>%
      do(get_query(.$filename, sql))
  } else {
    out <- foreach(i = db$position, .combine = rbind_list,
                   .packages = c("dplyr", "rplexos", "DBI", "RSQLite")) %dopar% {
      db %>%
        filter(position == i) %>%
        do(get_query(.$filename, sql))
    }
  }
  
  out %>% ungroup
}

#' Get list of available properties
#'
#' Produce a list of available properties, their units and the collections they belong to.
#' Additionally, a column is created for each scenario that indicates in how many databases
#' the property appears.
#'
#' @inheritParams query_master
#' @seealso \code{\link{plexos_open}} to create the PLEXOS database object
#' @family special queries
#' 
#' @examples
#' # Process the folder with the solution file provided by rplexos
#' location <- location_solution_rplexos()
#' process_folder(location)
#' 
#' # Query data
#' db <- plexos_open(location)
#' query_property(db)
#' 
#' @export
query_property <- function(db) {
  get_table_scenario(db, "property") %>%
    add_phase_names %>%
    group_by(phase_id, phase, is_summary, class_group, class, collection, property, unit, scenario) %>%
    summarize(n = n()) %>%
    tidyr::spread(scenario, n) %>%
    as.data.frame
}

#' Query configuration tables
#'
#' Get information from the \code{config} table, which includes: PLEXOS version, solution
#' date and time, machine and location of PLEXOS input database, model description and user
#' name. Additionally, it stores the version of rplexos used to process the PLEXOS database.
#'
#' @inheritParams query_master
#' @seealso \code{\link{plexos_open}} to create the PLEXOS database object
#' @family special queries
#'
#' @examples
#' # Process the folder with the solution file provided by rplexos
#' location <- location_solution_rplexos()
#' process_folder(location)
#' 
#' # Query data
#' db <- plexos_open(location)
#' query_config(db)
#'
#' @export
query_config <- function(db) {
  data <- get_table_scenario(db, "config")
    tidyr::spread(data, element, value) %>%
    as.data.frame
}

#' Query log file information
#'
#' During the processing of the PLEXOS databases, information from the log file is saved
#' into the database. This includes solution times and infeasibilities for the different phases.
#'
#' @inheritParams query_master
#' @seealso \code{\link{plexos_open}} to create the PLEXOS database object
#' @family special queries
#'
#' @examples
#' # Process the folder with the solution file provided by rplexos
#' location <- location_solution_rplexos()
#' process_folder(location)
#' 
#' # Query data
#' db <- plexos_open(location)
#' query_log(db)
#' query_log_steps(db)
#' 
#' @export
query_log <- function(db) {
  get_table_scenario(db, "log_info") %>%
    select(-position)
}

#' @rdname query_log
#' @export
query_log_steps <- function(db) {
  get_table_scenario(db, "log_steps") %>%
    select(-position)
}

# Query databases ***********************************************************************

#' Query data and aggregate data
#' 
#' This collection of functions retrieves data from the processed PLEXOS solutions and
#' returns it in a convenient format.
#' 
#' The family \code{query_*} returns the raw data in the databases, while \code{sum_*}
#' aggregates the data according to \code{columns}.
#'
#' The functions \code{*_day}, \code{*_week}, \code{*_month} and \code{*_year} are
#' shortcuts for the corresponding, \code{*_master} function.
#' 
#' The following is a list of valid items for \code{columns} and filtering. Additionally,
#' \code{time} can be specified for summary data (interval data always includes \code{time}).
#' \itemize{
#'   \item{\code{category}}
#'   \item{\code{property}}
#'   \item{\code{name} (default for columns)}
#'   \item{\code{parent} (automatically selected when \code{name} is selected)}
#'   \item{\code{category}}
#'   \item{\code{region} (only meaningful for generators)}
#'   \item{\code{zone} (only meaningful for generators)}
#'   \item{\code{period_type}}
#'   \item{\code{band}}
#'   \item{\code{sample}}
#'   \item{\code{timeslice}}
#' }
#' 
#' If defined, the \code{filter} parameter must be a \code{list}. The elements must be chracter
#' vectors and need to have a valid column name (see previous bullet points). For example, one
#' could define it as follows:
#' 
#' \code{filter = list(name = c("Generator1", "Generator2"), region = "Region1")}
#' 
#' To filter by time use the \code{time.range} parameter, instead of adding it as an entry in the
#' \code{filter} parameter. For example use \code{c("2015-03-14", "2015-03-15")} in your query.
#' Please note that the year/month/date starts at midnight (00:00:00).
#' 
#' If a scenario has multiple databases, the data will be aggregated automatically. If two or more
#' databases within the same scenario have overlapping time periods, the default is to select the
#' data from the last database (execute \code{summary(db)} so see the order). To change this behavior
#' set the global option \code{rplexos.tiebreak} to \code{first}, \code{last}, or \code{all} to
#' select data from the first database, the last one or keep all of them.
#' 
#' Multiple properties can be queried within a collection. If \code{prop} equals the widcard
#' \code{"*"}, all the properties within a collection are returned.
#' 
#' The parameter \code{multiply.time} allows to multiply values by interval duration (in hours) when
#' doing the sum of interval data. This can be used, for example, to obtain total energy (in MWh)
#' from power time series (in MW).
#' 
#' @param db PLEXOS database object
#' @param time character. Table to query from (interval, day, week, month, year)
#' @param col character. Collection to query
#' @param prop character vector. Property or properties to query
#' @param columns character. Data columns to query or aggregate by (defaults to \code{name})
#' @param time.range POSIXt or character. Range of dates of length 2 (given as date, datetime or character in 'ymdhms' or 'ymd' format)
#' @param filter list. Used to filter by data columns (see details)
#' @param phase integer. PLEXOS optimization phase (1-LT, 2-PASA, 3-MT, 4-ST)
#' @param multiply.time boolean. When summing interval data, provide the value multiplied by interval duration (See details).
#' @param ... parameters passed from shortcut functions to master (all except \code{time})
#' 
#' @return A data frame that contains data summarized/aggregated by scenario.
#' 
#' @seealso \code{\link{plexos_open}} to create the PLEXOS database object
#' @seealso \code{\link{query_sql}} to perform custom queries
#' 
#' @examples
#' # Process the folder with the solution file provided by rplexos
#' location <- location_solution_rplexos()
#' process_folder(location)
#' 
#' # Query data
#' db <- plexos_open(location)
#' query_day(db, "Generator", "Generation")
#' query_day(db, "Region", "*")
#' query_interval(db, "Generator", "Generation")
#' 
#' @export
#' @importFrom data.table data.table CJ
#' @importFrom foreach foreach %dopar%
query_master <- function(db, time, col, prop, columns = "name", time.range = NULL, filter = NULL, phase = 4) {
  # Check inputs
  check_rplexos(db)
  stopifnot(is.character(time), length(time) == 1L)
  stopifnot(is.character(col), length(col) == 1L)
  stopifnot(is.character(prop), is.character(columns))
  stopifnot(is.numeric(phase), length(phase) == 1L)
  if (!time %in% c("interval", "day", "week", "month", "year"))
    stop("'time' must be one of: interval, day, week, month or year", call. = FALSE)
  if(!phase %in% 1:4)
    stop("'phase' must be one of: 1 (LT), 2 (PASA), 3 (MT) or 4 (ST)", call. = FALSE)
  if(!all(columns %in% valid_columns()))
    stop("Incorrect column parameter. Use valid_columns() to get the full list.", call. = FALSE)
  
  # Key filter checks
  if (!is.null(filter)) {
    stopifnot(is.list(filter))
    if ("time" %in% names(filter))
      stop("time should not be an entry in filter. Use time.range instead.", call. = FALSE)
    if(!all(names(filter) %in% valid_columns()))
      stop("The names in 'filter' must correspond to correct columns. Use valid_columns() to get the full list.", call. = FALSE)
  }
  
  # Time range checks
  if (!is.null(time.range)) {
    stopifnot(length(time.range) == 2L)
    
    if (inherits(time.range, "POSIXt")) {
      time.range2 <- time.range
    } else {
      time.range2 <- c(NA, NA)
      
      if (inherits(time.range, "character")) {
        time.range2 <- lubridate::parse_date_time(time.range, c("ymdhms", "ymd"), quiet = TRUE)
      }
      
      if(any(is.na(time.range2)))
        stop("time.range must be POSIXt or character with 'ymdhms' or 'ymd' formats", call. = FALSE)
    }
    
    # Convert dates to ymdhms format, so that queries work correctly
    time.range <- format(time.range2, "%Y-%m-%d %H:%M:%S")
  }
  
  ### BEGIN: Master query checks
  
  # Get list of properties for the collection
  is.summ <- ifelse(identical(time, "interval"), 0, 1)
  is.summ.txt <- ifelse(identical(time, "interval"), "interval", "summary")
  res <- rbind_all(db$properties) %>%
    filter(collection == col, is_summary == is.summ, phase_id == phase)
  
  # Check that collection is valid
  if (nrow(res) == 0L) {
    stop("Collection '", col, "' is not valid for ",
         is.summ.txt, " data and phase '", phase, "'.\n",
         "   Use query_property() for list of collections and properties.",
         call. = FALSE)
  }
  
  # Checks if properties are valid
  if (!identical(prop, "*")) {
    invalid.prop <- setdiff(prop, res$property)
     if (length(invalid.prop) > 0L) {
      stop("Properties ", paste0("'", invalid.prop, "'", collapse = ", "), " in collection '", col,
           "' are not valid for ", is.summ.txt, " data and phase '", phase, "'.\n",
           "   Use query_property() for list of available collections and properties.",
           call. = FALSE)
     }
    
    # Filter properties
    res <- res %>%
      filter(property %in% prop)
  }
  
  # Find if the data is going to have multiple sample, timeslices or bands
  res2 <- res %>%
    ungroup() %>%
    summarize(is_multi_band      = max(count_band) > 1,
              is_multi_sample    = max(count_sample) > 1,
              is_multi_timeslice = max(count_timeslice) > 1)
  if (res2$is_multi_timeslice)
    columns <- c(setdiff(columns, "timeslice"), "timeslice")
  if (res2$is_multi_band)
    columns <- c(setdiff(columns, "band"), "band")
  if (res2$is_multi_sample)
    columns <- c(setdiff(columns, "sample"), "sample")
  
  # Columns should not include collection and property; they are always reported
  columns <- setdiff(columns, c("collection", "property"))
  
  # If columns include name, add parent automatically
  if ("name" %in% columns)
    columns <- c("name", "parent", setdiff(columns, c("name", "parent")))
  
  ### END: Master query checks
  
  # Query data for each property
  db2 <- db %>%
    group_by(scenario, position)
  
  if (!is_parallel_rplexos()) {
    out <- db2 %>%
      do(query_master_each(., time, col, prop, columns, time.range, filter, phase))
  } else {
    out <- foreach(i = db2$position, .combine = rbind_list,
                   .packages = c("dplyr", "rplexos", "DBI", "RSQLite")) %dopar% {
      db2 %>%
        filter(position == i) %>%
        do(query_master_each(., time, col, prop, columns, time.range, filter, phase))
    }
  }
  
  # Ungroup results
  out <- out %>% ungroup
  
  # Return empty dataframe if no results were returned
  if (nrow(out) == 0) {
    warning("Query returned no results", call. = FALSE)
    return(data.frame())
  }
  
  # Check if any scenario is missing from the results
  missing.scenario <- setdiff(unique(db$scenario), unique(out$scenario))
  if (length(missing.scenario) >= 1L) {
    warning("Query returned no results for scenario(s): ",
            paste(missing.scenario, collapse = ", "),
            call. = FALSE)
  }
  
  # Solve ties if they exist
  out <- out %>%
    solve_ties()
  
  out
}

#' Open one database and query data
#'
#' Used internally by \code{\link{query_master}}. Use that function instead to access data. The
#' use of this function is not recommended
#'
#' @inheritParams query_master
#'
#' @seealso \code{\link{query_master}} to query PLEXOS solutions
#' 
#' @keywords internal
#' @export
query_master_each <- function(db, time, col, prop, columns = "name", time.range = NULL, filter = NULL, phase = 4) {
  # Open connection
  thesql <- src_sqlite(db$filename, create = FALSE)
  
  if (!identical(time, "interval")) {
    # Query interval data
    if (identical(prop, "*")) {
      out <- tbl(thesql, time)
    } else if (length(prop) == 1L) {
      # Workaround for a problem with dplyr
      out <- tbl(thesql, time) %>% filter(property == prop)
    } else {
      out <- tbl(thesql, time) %>% filter(property %in% prop)
    }
    
    out <- out %>%
      filter(collection == col, phase_id == phase) %>%
      filter_rplexos(filter) %>%
      filter_rplexos_time(time.range) %>%
      select_rplexos(columns, add.key = FALSE) %>%
      collect %>%
      mutate(time = lubridate::ymd_hms(time, quiet = TRUE))
  } else {
    # Query interval data
    # Get the table names that store the data
    if (identical(prop, "*")) {
      t.name <- db$properties[[1]]
    } else if(length(prop) == 1L) {
      t.name <- db$properties[[1]] %>% filter(property == prop)
    } else {
      t.name <- db$properties[[1]] %>% filter(property %in% prop)
    }
    t.name <- t.name %>%
      filter(collection == col, phase_id == phase, is_summary == 0) %>%
      select(collection, property, table_name) %>%
      mutate(table_name = gsub("data_interval_", "", table_name))
    
    # If t.name is empty (data is not available), return an empty data frame
    if (nrow(t.name) == 0L) {
      DBI::dbDisconnect(thesql$con)
      return(data.frame())
    }
    
    # Get max/min time existing in the table to be queried
    #   In case time table has more time stamps than those in the dataset
    time.limit <- t.name %>%
      group_by(collection, property) %>%
      do(
        tbl(thesql, .$table_name) %>%
          filter(phase_id == phase) %>%
          summarize(time_from = min(time_from), time_to = max(time_to)) %>%
          collect
      )
    min.time.data <- min(time.limit$time_from)
    max.time.data <- max(time.limit$time_to)
    
    # Collect time data
    time.data <- tbl(thesql, "time") %>%
      filter(phase_id == phase) %>%
      filter(between(time, min.time.data, max.time.data)) %>%
      filter_rplexos_time(time.range) %>%
      select(time) %>%
      collect
    
    # If time data is empty, return an empty data frame
    if (nrow(time.data) == 0L) {
      DBI::dbDisconnect(thesql$con)
      return(data.frame())
    }
    
    # Convert into R time-data format
    time.data$time <- lubridate::ymd_hms(time.data$time, quiet = TRUE)
    
    # Get interval data
    out1 <- t.name %>%
      group_by(collection, property) %>%
      do(tbl(thesql, .$table_name) %>%
           filter(phase_id == phase) %>%
           filter_rplexos(filter) %>%
           filter_rplexos_time(time.range, modified = TRUE) %>%
           select(-time_to) %>%
           rename(time = time_from) %>%
           select_rplexos(columns, add.key = TRUE) %>%
           collect
      ) %>%
      ungroup %>%
      mutate(time = lubridate::ymd_hms(time, quiet = TRUE))
  
    # Expand data
    #   This will be easier when dplyr supports rolling joins
    out2 <- data.table(out1, key = "key,time")
    cj2 <- CJ(key = unique(out1$key), time = time.data$time)
    
    out3 <- out2[cj2, roll = TRUE]
    out <- out3 %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      select(-key)
    
    # Restore time zone
    attributes(out$time) <- attributes(time.data$time)
  }
  
  # Disconnect database
  DBI::dbDisconnect(thesql$con)
  
  # Return value
  return(out)
}


# Deal with repeats
solve_ties <- function(x, opt = getOption("rplexos.tiebreak")) {
  # Get option to see how to deal with ties (defaults to last)
  if (is.null(opt)) {
    opt <- "last"
  } else if (!opt %in% c("first", "last", "all")) {
    warning("Invalid 'rplexos.tiebreak' option (must be one of: first, last, all). Using last instead", call. = FALSE)
    opt <- "last"
  }
  
  if (opt %in% c("first", "last")) {
    # Group by time
    x2 <- x %>%
      ungroup() %>%
      group_by(scenario, time)
    
    if (identical(opt, "last")) {
      # If there are repeats, use the latter entry
      x2 <- x2 %>%
        filter(position == max(position))
    } else {
      # If there are repeats, use the latter entry
      x2 <- x2 %>%
        filter(position == min(position))
    }
    
    # Ungroup and delete path column
    x2 <- x2 %>%
      ungroup() %>%
      select(-position)
  }
  
  x2
}

#' @rdname query_master
#' @export
query_interval <- function(db, ...) query_master(db, "interval", ...)
#' @rdname query_master
#' @export
query_day      <- function(db, ...) query_master(db, "day", ...)
#' @rdname query_master
#' @export
query_week     <- function(db, ...) query_master(db, "week", ...)
#' @rdname query_master
#' @export
query_month    <- function(db, ...) query_master(db, "month", ...)
#' @rdname query_master
#' @export
query_year     <- function(db, ...) query_master(db, "year", ...)


# Aggregation ***************************************************************************

#' @rdname query_master
#' @export
sum_master <- function(db, time, col, prop, columns = "name", time.range = NULL, filter = NULL, phase = 4, multiply.time = FALSE) {
  # Check inputs to unique
  stopifnot(is.logical(multiply.time), length(multiply.time) == 1L)

  # Make sure to include time
  columns2 <- c(setdiff(columns, "time"), "time")
  
  # Run query_master to get the raw data
  out <- query_master(db, time, col, prop, columns2, time.range, filter, phase)
  
  # If empty query is returned, return empty data.frame
  if(nrow(out) == 0L)
    return(data.frame())
  
  if (identical(time, "interval") && (!"time" %in% columns) && multiply.time) {
    # Get length of intervals in hours
    times <- get_table_scenario(db, "time")
    delta <- times %>%
      group_by(scenario) %>%
      mutate(time = lubridate::ymd_hms(time, quiet = TRUE)) %>%
      summarize(interval = difftime(lead(time), time, units = "hours") %>%
                  min(na.rm = TRUE) %>%
                  as.numeric)
    
    # Add interval duration to the sum
    out <- out %>%
      inner_join(delta, by = "scenario") %>%
      group_by_char(c("scenario", "collection", "property", columns)) %>%
      summarise(value = sum(value * interval))  
    
    # If unit is a column, modify column
    if ("unit" %in% names(out)) {
      out <- out %>%
        mutate(unit = paste(unit, "* h"))
    }
  } else {
    # Sum values
    out <- out %>%
      group_by_char(c("scenario", "collection", "property", columns)) %>%
      summarise(value = sum(value))  
  }
  
  out
}

#' @rdname query_master
#' @export
sum_interval <- function(db, ...) sum_master(db, "interval", ...)
#' @rdname query_master
#' @export
sum_day      <- function(db, ...) sum_master(db, "day", ...)
#' @rdname query_master
#' @export
sum_week     <- function(db, ...) sum_master(db, "week", ...)
#' @rdname query_master
#' @export
sum_month    <- function(db, ...) sum_master(db, "month", ...)
#' @rdname query_master
#' @export
sum_year     <- function(db, ...) sum_master(db, "year", ...)


# Filtering *****************************************************************************

# Time filter
filter_rplexos_time <- function(out, time.range, modified = FALSE) {
  # Do nothing if time.range is empty
  if (!is.null(time.range)) {
    if (modified) {
      out <- filter(out, time_from <= time.range[2], time_to >= time.range[1])
    } else {
      out <- filter(out, between(time, time.range[1], time.range[2]))
    }
  }
  
  out
}

# Other filters
filter_rplexos <- function(out, filt) {
  # Do nothing if filter is empty
  if (is.null(filt))
    return(out)
  if (length(filt) == 0L)
    return(out)
  
  # Write the condition as text
  vals <- lapply(filt, function(x)
    paste0("\"", x, "\"", collapse = ", ")) %>%
    paste0("c(", ., ")")
  cons <- ifelse(lapply(filt, length) == 1L, "==", "%in%")
  cond <- paste(names(filt), cons, vals)
  
  # Apply condition
  out %>% filter_(.dots = cond)
}

# Dynamically select the columns
select_rplexos <- function(x, columns, add.key) {
  if (add.key) {
    columns.dots <- c("key", "unit", setdiff(columns, "time"), "time", "value")
  } else {
    columns.dots <- c("collection", "property", "unit", setdiff(columns, "time"), "time", "value")
  }
  
  columns.dots <- columns.dots %>%
    as.list %>%
    lapply(as.symbol)
  
  select_(x, .dots = columns.dots)
}
