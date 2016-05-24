#' Extract a Time Series from a Data.table 
#' 
#' This function extracts time series from data.table columns and returns object of class ts. 
#' @author Matthias Bannert, Gabriel Bucur
#' @param data_table a data.table 
#' @param time_column character name of the column which contains the time index
#' @param group_list list or NULL
#' @param freq integer value either 4 denoting quarterly frequency or 12 denoting quarterly frequency
#' @param item character name of the column which contains the item
#' that is extracted from the data.table
#' @param variable character name of the variable selected 
#' @param prefix character prefix attached to the dynamically generated key string to identify the time series. Recommend key format:
#' ISOcountry.provider.source.aggregationLevel.selectedGroup.variable.item
#' @export 
extractTimeSeries <- function(data_table, time_column, group_list, freq,
                              item, variable, prefix = "CH.KOF.IND") {
  
  old_key <- key(data_table)
  setkeyv(data_table, names(group_list))
  
  # Cross Join is used because J is not a function
  # and cannot be do.called
  if (is.null(group_list)) {
    dt_subset <- data_table
  } else {
    dt_subset <- data_table[do.call("CJ", group_list)]
  }
  
  start_date <- min(as.Date(dt_subset[[time_column]]))
  year <- as.numeric(format(start_date, "%Y"))
  
  if (freq == 4) {
    period <- (as.numeric(format(start_date, "%m")) - 1) / 3 + 1
  } else if(freq == 12) {
    period <- as.numeric(format(start_date, "%m"))
  } else {
    stop("Not a standard frequency.")
  }
  
  setkeyv(dt_subset, time_column)
  setkeyv(data_table, old_key)
  
  # create output time series in order to add attributes
  out_ts <- ts(dt_subset[[item]],
               start = c(year, period),
               frequency = freq)
  
  # generate key parts seperately to make key pasting easier to read
  if (is.null(group_list)) {
    group_level <- NULL
    group_value <- NULL
  } else {
    group_level <- paste(setdiff(names(group_list),
                                 time_column), collapse = "_")
    group_value <- paste(unlist(group_list), collapse = "_")
  }
  
  ts_key <- paste(c(prefix, group_level, group_value,
                    variable, item), collapse = ".")
  
  # substitute : with _
  ts_key <- gsub("\\:","_",ts_key)
  
  attr(out_ts,"ts_key") <- toupper(ts_key)
  out_ts  
  
}
