setClass("MTCDevice", representation(rawdata = "list", metadata = "list"), contains = "MTCCycle")

setMethod("initialize", "MTCDevice", function(.Object, data_item_list = list(), rawdata = list(), 
                                              metadata = list(), device_uuid = NULL){
 
    .Object@rawdata = rawdata
    .Object@metadata = metadata
    .Object@data_item_list = data_item_list
    .Object@device_uuid = device_uuid
  return(.Object)
})


#' Create an MTC device object from a merged time series data frame
#' 
#' @param merged_device An existing object of MTCDevice Class
#' @param device_uuid UUID to be given to the device
#' @examples 
#' data("example_mtc_device")
#' merged_device = merge(example_mtc_device)
#' create_mtc_device_from_ts(merged_device)
#' 
#' @export
create_mtc_device_from_ts <- function(merged_device, device_uuid = "unmerged_device"){
  data_item_names = setdiff(names(merged_device), "timestamp")
  data_item_list = lapply(data_item_names, function(x){
    temp_df = data.frame(timestamp = merged_device$timestamp, value = merged_device[[x]]) %>% 
      clean_reduntant_rows()
    new("MTCDataItem", temp_df,
        data_type = ifelse(test = is.numeric(temp_df$value), yes = 'Sample', no = 'Event'),
        path = x, dataSource = "Unmerged", xmlID = "") 
  })
  names(data_item_list) = data_item_names
  new('MTCDevice', rawdata = list(), data_item_list = data_item_list, device_uuid = device_uuid)
  
}