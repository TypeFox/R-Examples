# Author: David Barker, Justin Hemann <support@causata.com>

Connect <- function(sql.server.host, sql.server.port, sql.username, sql.password, group=NULL, verbose=FALSE) {
  # hostname: the host of the query broker
  # port: the port that the MySql interface is running on
  # username: causata username
  # password: causata password.
  # group: group of data to load from config file
  
  # create a list of arguments supplied in the function call
  arglist.fun <- as.list(match.call()[-1])
  # create another argument list extracted from causata config file
  arglist.yaml <- LoadCausataConfig(group)
  # combine the arguments from these two sources into a single list
  # if an argument exists in both lists then the function arguments take precendence since they are first
  arglist <- c(arglist.fun, arglist.yaml)
  
  this <- list();
  class(this) <- "Connect"
  if (verbose){
    cat(paste("Connecting to:", arglist$sql.server.host, arglist$sql.server.port, arglist$sql.username, arglist$sql.password, "\n", collapse=", "))
  }
  
  this$conn <- dbConnect(MySQL(), user=arglist$sql.username, password=arglist$sql.password, dbname="causata", host=arglist$sql.server.host, port=arglist$sql.server.port, client.flag=CLIENT_MULTI_STATEMENTS)
  this <- GetMetadata(this)
  return(this)
}
is.Connect <- function(this) inherits(this, "Connect")


LoadCausataConfig <- function(group){
  config.file <- path.expand("~/.causata-config.yaml")
  # if the group was supplied then load default values from the yaml file
  if (!is.null(group)){
    # group was provided, try to load config data from YAML file
    result <- try(yaml.load_file(config.file))
    if(class(result) == "try-error"){
      # loading the YAML file resulted in an error
      stop("Could not load configuration data from ", config.file)
    } else {
      # the YAML config file was loaded successfully. See if it has the right group name in the list returned by yaml.load_file
      if (!(group %in% names(result))){
        stop("Could not load the group ", group, " from ", config.file)
      }
      # the group exists, so copy it to argument list
      return(result[[group]])
    }
  } else {
    # no group was provided, return an empty list
    return(list())
  }
}


#
# define generic functions
#
GetRawData <- function(conn, ...) {
  # generic to get a raw dataframe from a query
  UseMethod("GetRawData")
}
GetData <- function(conn, ...) {
  # generic to get a dataframe from a query
  UseMethod("GetData")
}
GetCausataData <- function(conn, ...) {
  # generic to get a dataframe from a query
  UseMethod("GetCausataData")
}
GetMetadata <- function(conn, ...) {
  UseMethod("GetMetadata")
}
Close <- function(conn, ...) {
  # Closes this connection
  UseMethod("Close", conn)
}
GetNames <- function(this, ... ) {
  UseMethod("GetNames", this)
}

GetNames.Connect <- function(this, kind, ...){
  # return the list of causata variable display names
  # Variables with an archived flag are NOT returned
  if (kind=="system"){
    return(this$variables$NAME[this$variables$IS_ARCHIVED == 0])
  } else if (kind=="display"){
    return(this$variables$DISPLAY_NAME[this$variables$IS_ARCHIVED == 0])
  } else {
    stop("Invalid kind argument, must be 'system' or 'display': ", kind)
  }
}

# Gets the raw result of this query
#
GetRawData.Connect <- function(conn, query, ...) {
  sql <- as.character(query)
  #cat(paste("Executing query:", sql))
  dataframe <- dbGetQuery(conn$conn, sql)
  if (nrow(dataframe)==0){
    stop("Data frame with zero rows returned, query may not be valid.")
  } else {
    return(dataframe)
  }
}

# Gets the result of the given query with character and boolean data converted to factors.
#
GetData.Connect <- function(conn, query, ...) {
  Reclass.Connect(conn, GetRawData.Connect(conn, query))
}

# Returns a CausataData object from the given query
#
GetCausataData.Connect <- function(conn, query, dependent.variable, ...) {
  reclassed <- GetData.Connect(conn, query)
  data <- CausataData(
    dataframe=reclassed,
    dependent.variable,
    query=query
  )
  return(data)
}

#
# private functions
#


GetMetadata.Connect <- function(this, ...) {
  #
  # Returns a dataframe with variables and types.
  # both are strings, the possible types are:
  # date, string, numeric
  #
  
  # extract metadata from SQL
  this$variables    <- dbGetQuery(this$conn, "SELECT * FROM information_schema.causata_variables")
  this$labels       <- dbGetQuery(this$conn, "SELECT * FROM information_schema.causata_variable_label")
  this$timePoints   <- dbGetQuery(this$conn, "SELECT * FROM information_schema.causata_time_points")
  this$timeRanges   <- dbGetQuery(this$conn, "SELECT * FROM information_schema.causata_intervals")
  this$dependencies <- dbGetQuery(this$conn, "SELECT * FROM information_schema.causata_variable_dependencies")
  this$eventTypes   <- dbGetQuery(this$conn, "SELECT * FROM information_schema.causata_event_types")
  this$eventAttributes <- dbGetQuery(this$conn, "SELECT * FROM information_schema.causata_event_attributes")
  return(this)
}



Reclass.Connect <- function(this, df, unique.rate=0.95){
  # This function changes the classes of variables in a data frame.
  # SQL queries return Causata boolean data as numeric, and Causata text as characters.
  # We want Causata booleans and text to be returned as factors.  
  # - this: a Connect object
  # - df: a data frame extracted from Causata through the SQL interface
  
  # rename variables to R format
  names(df) <- CausataToRNames(names(df))
  
  # build a vector of variable names in the data frame with the time period removed.  The names will 
  # be joined to a time period separated by a dollar sign, name$period, so split at the $
  # and keep only the part of the string before the $
  variables.df <- sapply( strsplit(names(df),"__"), "[", 1)
  variables.schema <- CausataToRNames(this$variables$NAME)
  # determine which variables from the dataframe are found in the set of variables from the causata schema
  idx <- variables.df %in% variables.schema
  if (sum(idx) == 0){
    # the dataframe does not contain variables found in the schema, return unaltered dataframe
    return(df)
  }
  
  # at least one variable from the data frame was found in the schema, check types
  variables.match <- variables.df[idx]
  
  # set attribute for class.causata, loop for each column in data frame
  class.causata.vec <- rep('NULL', length(idx))
  for (i in 1:length(idx)){
    if (idx[i]){
      # this column is in the schema, look up the class
      var.schema <- variables.df[i] # variable as it appears in schema
      class.causata <- this$variables$DATA_TYPE[ variables.schema == var.schema ]
      
      # use the set function from data.table to copy by reference, this minimizes memory use
      set(df, j=i, value = ReclassVector(df[, i], class.causata, unique.rate))
    }
  }
  return(df)
}


ReclassVector <- function(xin, class.causata, unique.rate){
  # use the attribute to control how data is modified or passed unchanged
  if (class.causata == "NULL"){
    # no causata class assigned, return unaltered data
    return(xin)
    
  } else if (class.causata == "Boolean"){
    # SQL returns booleans as 1/0, so remap to a factor that is true/false
    xout <- rep("false", length(xin))
    # set missing values if present
    xout[is.na(xin)] <- NA
    xout[xin==1] <- "true"
    return(as.factor(xout))
    
  } else if (class.causata %in% c("Integer","Long","Float","Double")){
    # convert integer to numeric
    return(as.numeric(xin))
    
  } else if (class.causata == "String"){
    # determine if values are all (or almost all) unique
    if ((length(unique(xin)) / length(xin)) < unique.rate) {
      # the rate of unique values is below threshold so convert to factor
      return(as.factor(xin))
    } else {
      # rate of unique values is high, return string
      return(xin)
    }
  } else if (class.causata == "Date"){
    # return dates unaltered
    return(xin)
    
  } else {
    warning("Invalid data class ", class.causata)
    return(xin)
  }
}


Close.Connect <- function(conn, ...) {
  dbDisconnect(conn$conn)
}
