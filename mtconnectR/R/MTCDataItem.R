
## MTCDataItem Class methods
setClass("MTCDataItem",  representation(data = "data.frame", data_type = "character", path = "character", dataSource="character", xmlID="character"), prototype("data_type" = "Event", "path" = ""))

#' Get data from the object in a data frame form
#' @param .Object An MTC Object
#' @param pattern OPTIONAL can be used to query specific data items
#' @examples 
#' data("example_mtc_data_item")
#' getData(example_mtc_data_item)
setGeneric("getData", function(.Object, pattern){standardGeneric("getData")})

setValidity("MTCDataItem", function(object)
{
  if(!is(object@data$timestamp, "POSIXct")) {	
    stop("timestamp objects not of class POSIXct", structure = str(object))
  }
  if(object@data_type == "Sample") {
    if(!is(object@data$value, "numeric")) {
      stop("Sample objects in MTCDataItem not of class numeric", structure = str(object))
    }
  } else if(!is(object@data$value, "character")) {
    stop("Non - Sample objects in MTCDataItem not of class character", structure = str(object))
  }
  return(TRUE)
})

setMethod("initialize", "MTCDataItem", function(.Object, data, data_type="Event", path=NULL, dataSource="JSON", xmlID = "No XML ID found"){
  
  .Object@data_type = data_type
  .Object@path = path
  rownames(data) <- NULL
  .Object@data = as.data.frame(data)
  .Object@dataSource = dataSource
  
  if ((data_type == "Sample") && ("value" %in% names(.Object@data))) {
    .Object@data$value[.Object@data$value == "UNAVAILABLE"] = NA_real_
    .Object@data$value = as.numeric(.Object@data$value)
  }
  
  if(data_type == "Event")  .Object@data$value = as.character(.Object@data$value) 
  
  .Object@xmlID <- xmlID
  
  return(.Object)
})

#' Show a quick summary of the MTCDataItem
#' 
#' @param object The MTCDataItem object
#' @examples 
#' data("example_mtc_data_item")
#' summary(example_mtc_data_item)
#' @export
setMethod("summary", "MTCDataItem", function(object){
  output = data.frame((object@path), nrow(object@data), object@data$timestamp[1], tail(object@data$timestamp, 1), object@data_type)
  names(output) = c("path", "Records", "start", "end", "data_type")
  return(output)
})

#' Get Data from the Object as a data.frame
#' 
#' @param .Object object of MTCDataItem Class
#' @examples 
#' data("example_mtc_data_item")
#' getData(example_mtc_data_item)
#' @export
setMethod("getData", "MTCDataItem", function( .Object){
  return(.Object@data)
})
