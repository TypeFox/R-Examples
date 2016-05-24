
# Common functions across the package

#' @importFrom data.table rbindlist
#' @importFrom stringr str_split
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract
#' @importFrom utils tail
#' @importFrom plyr ldply
#' @importFrom dplyr arrange_
#' @importFrom dplyr '%>%' 
#' @importFrom dplyr mutate 
#' @importFrom dplyr select_
#' @importFrom dplyr group_by
#' @importFrom dplyr do
#' @importFrom dplyr select
#' @import methods
NULL

#' Example data set showing MTC Log data
#'
#' A manually created dataset showing a log data file, parsed and read into R. The columns are 
#' \itemize{
#'   \item timestamp. Timestamp of the event
#'   \item data_item_name Name of the data Item from the delimited MTC data. Can be empty.
#'   \item value of the data item
#' }
#'
#' @format A data frame some rows and 3 variables
"example_dmtcd"

#' Example data set showing Xpaths from a device XML
#'
#' Dataset showing a parsed DeviceXML file showing all the XPaths and the properties
#' \itemize{
#'   \item id ID of the data item
#'   \item name Name of the data Item from the delimited MTC data. Can be empty.
#'   \item type MTC Type of the data item
#'   \item category MTC Category of the data item
#'   \item subType MTC subType of the data item. Can be emoty
#'   \item xpath xpath showing the truncated path to the particular data item in the device XML
#' }
#'
#' @format A data frame some rows and 6 variables
"example_xpath_info"

#' Example data set showing a MTConnect Device
#'
#' The data can be accessed using the @ function. The slots are:
#' \itemize{
#'   \item rawdata Original delimited MTC data (parsed from which the data was created)
#'   \item metadata Metadata (if any) for the device
#'   \item data_item_list processed data showing each data item as a separate device
#'   \item device_uuid UUID of the device
#' }
#'
#' @format An MTCDevice data item
"example_mtc_device"

#' A bigger example data set showing a MTConnect Device with path position and conditions
#'
#' The data can be accessed using the @ function. The slots are:
#' \itemize{
#'   \item rawdata Original delimited MTC data (parsed from which the data was created)
#'   \item metadata Metadata (if any) for the device
#'   \item data_item_list processed data showing each data item as a separate device
#'   \item device_uuid UUID of the device
#' }
#'
#' @format An MTCDevice data item
"example_mtc_device_2"


#' Example data set showing a MTConnect DataItem
#'
#' The data can be accessed using the @ function. The slots are:
#' \itemize{
#'   \item data Data for a single data item at a data.frame in timestamp, value format
#'   \item data_type Type of Data - can be even or sample
#'   \item path XML Xpath
#'   \item data_source source from which the data item was created
#'   \item xmlID id of the data item in the devices XML
#' }
#'
#' @format An MTCDevice data item
"example_mtc_data_item"



#' Convert Time Series to Intervals
#' 
#' Function to convert a continuous time series data to interval data.
#' The last row which goes to infinity can be deleted, else will be given dump value.
#' 
#' @param df a data frame with continuous time series data
#' @param endtime_lastrow POSIXct value for the last row. Defaults to NA
#' @param arrange_cols whether to rearrange the columns adding the new columns at the front
#' @param time_colname column name of the timestamp variable
#' @param round_duration number of decimals to rounds the duration to. Defaults
#'  to 2. If no rounding required, give NULL.
#' @seealso \code{\link{convert_interval_to_ts}}
#' @export
#' @examples 
#' 
#' ts_data = data.frame(ts = as.POSIXct(c(0.5, 1, 1.008, 1.011),  tz = 'UTC', origin = "1970-01-01"),
#'                      x = c("a", "b", "c", "d"), y = c("e", "e", "e", "f"))
#' convert_ts_to_interval(ts_data, time_colname = "ts", endtime_lastrow = ts_data$ts[1] + 10)
convert_ts_to_interval <- function(df, endtime_lastrow = as.POSIXct(NA), arrange_cols = T,
                                   time_colname = 'timestamp', round_duration = 2)
{
  start_col = which(colnames(df) == time_colname)
  if (!is.null(endtime_lastrow)) df$end = endtime_lastrow else
    df$end = df[,start_col]  # Dump values for the the End times
  
  df = df[order(df[,start_col]), ]
  
  if (nrow(df) > 1) {
    df$end[1:(nrow(df) - 1)] = df[,time_colname][2:nrow(df)]
    if (is.null(endtime_lastrow))
    {
      df = df[-nrow(df),] #Deleting the last row which goes to infinity
    }else df$end[nrow(df)] = endtime_lastrow
  }
  df$duration = as.numeric(df$end) - as.numeric(df[,time_colname]) #Duration of each process
  if (!is.null(round_duration))
    df$duration <- round(df$duration, round_duration)
  
  colnames(df)[start_col] = 'start'
  if (arrange_cols == T) df = df[,c(start_col, (ncol(df) - 1 ), ncol(df), setdiff(1:(ncol(df) - 2), start_col))]
  rownames(df) = NULL
  return(df)
}

#' Convert Interval to Time Series
#' 
#' Basically reverse the effect of \code{\link{convert_ts_to_interval}}
#' 
#' @param df data.frame in start, end, duration, value1, value2,...
#' @param time_colname name of the time column
#' 
#' @seealso \code{\link{convert_ts_to_interval}}
#' @export
#' @examples 
#' test_interval = 
#'   data.frame(start = as.POSIXct(c(0.5, 1, 1.008, 1.011),  tz = 'CST6CDT', origin = "1970-01-01"),
#'              end   = as.POSIXct(c(1, 1.008, 1.011, 2),  tz = 'CST6CDT', origin = "1970-01-01"),
#'              duration = c(0.50, 0.01, 0.00, 0.99),
# '             x     = c("a", "b", "c", "d"), 
#'              y     = c("e", "e", "e", "f"))
#' convert_interval_to_ts(test_interval)
convert_interval_to_ts <- function(df, time_colname = 'start')
{
  data_n = df
  data_n$end = data_n$duration = NULL
  data_n[nrow(data_n) + 1,] = NA
  rownames(data_n) = NULL
  data_n$start[nrow(data_n)] = tail(df$end, 1)
  
  data_n[nrow(data_n), is.na( data_n[nrow(data_n), ])] = NA
  colnames(data_n)[1] = "timestamp"
  
  if (all(is.na(data_n[nrow(data_n),])) & nrow(data_n) > 1) data_n = data_n[1:(nrow(data_n)-1),]
  
  return(data_n)
}


#' Removes Redundant Rows in a data frame assuming statefullness
#' 
#' @param df data.frame in timestamp, value1, value2,...
#' @param clean_colname name of the column to clean as basis
#' @param echo whether to return messages or not
#' 
#' @export
#' @examples 
#' test_interval = 
#'   data.frame(timestamp = as.POSIXct(c(0.5, 1, 1.008, 1.011), origin = "1970-01-01"),
#'             x     = c("a", "b", "b", "b"), 
#'              y     = c("e", "e", "e", "f"))
#' clean_reduntant_rows(test_interval, "x")
clean_reduntant_rows = function(df, clean_colname = "value", echo = F) {
  df = data.frame(df)
  clean_col = grep(clean_colname, names(df))
  
  if (echo) message(paste("Cleaning table with ", paste(names(df)[clean_col], collapse=","), " as basis..."))
  
  if (length(clean_col) == 0 ) message("No Columns match the required pattern for cleaning!")
  if (length(clean_col) > 1 )  pasted_vector = do.call(paste, df[clean_col]) else
    pasted_vector = df[[clean_col]]
  
  data_n = diff(as.numeric(as.factor(pasted_vector)))
  data_n[is.na(data_n)] = -100
  selected_rows = c(T, abs(data_n)!= 0)
  df = df[selected_rows,]
  rownames(df) = NULL
  return(df)
}