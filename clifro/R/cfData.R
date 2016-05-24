#' @include dataFrame.R
NULL

# cfData Class ------------------------------------------------------------

#' @importFrom methods setClass
setClass("cfData", slots = c(dt_name = "character",
                             dt_type = "character"),
         contains = "dataFrame")

# cfData Objects ----------------------------------------------------------

#' @importFrom methods setClass
setClass("cfWind", slots = c(data_label = "character"),
         contains = "cfData")

#' @importFrom methods setClass
setClass("cfRain", slots = c(data_label = "character"),
         contains = "cfData")

#' @importFrom methods setClass
setClass("cfScreenObs", slots = c(data_label = "character",
                                  plot_label = "call"),
         contains = "cfData")

#' @importFrom methods setClass
setClass("cfTemp", slots = c(data_label = "character",
                             plot_label = "call"),
         contains = "cfData")

#' @importFrom methods setClass
setClass("cfEarthTemp", slots = c(data_label = "character",
                                  plot_label = "call"),
         contains = "cfData")

#' @importFrom methods setClass
setClass("cfSunshine", slots = c(data_label = "character"),
         contains = "cfData")

#' @importFrom methods setClass
setClass("cfPressure", slots = c(data_label = "character"),
         contains = "cfData")

#' @importFrom methods setClass
setClass("cfOther", slots = c(data_label = "character"),
         contains = "cfData")

# Create individual S4 objects based on the datatype
#
# This is the initialisation function for all the different types of data 
# available from CliFlo. 
#
# object a cfData object
#' @importFrom methods new

create_object = function(object){
  stopifnot(is(object, "cfData"))
  
  ## Wind
  if (tolower(object@dt_name) %in% c("surface wind", "max gust")){
    dl = if (grepl("wind run", tolower(object@dt_type))){
      "9am wind run (km)"
    } else {
      spd_col = pmatch("Speed", object@names)
      dt_units = strsplit(object@names[spd_col], "(", 
                          fixed = TRUE)[[1]][2]
      paste(object@dt_type, tolower(object@dt_name), 
            paste0("(", dt_units))
    }
    return(new("cfWind", data_label = dl, object))
  }
  
  ## Precipitation
  
  if (object@dt_name == "Rain")
    return(new("cfRain",
               data_label = paste(object@dt_type, "Rain (mm)"), object))
  
  if (object@dt_name == "Snow")
    return(new("cfRain", data_label = "Snow", object))
  
  ## ScreenObs
  if (tolower(object@dt_name) == "screenobs"){
    
    
    return(new("cfScreenObs",
               data_label = paste(object@dt_type, tolower(object@dt_name)),
               plot_label = bquote(.(object@dt_type) ~ screen ~ observations ~ (degree * C)),
               object))
  }
  
  ## Temperature and Humidity
  if (tolower(object@dt_name) == "max_min"){
    return(new("cfTemp", 
               data_label = paste(object@dt_type, "maximum/minimum temperature"),
               plot_label = bquote(.(object@dt_type) ~ temperature ~ (degree * C)),
               object))
  }
  
  if (tolower(object@dt_name) == "earth temp")
    return(new("cfEarthTemp",
               data_label = paste("Earth temperature at", object@dt_type, "depth"),
               plot_label = bquote(Earth ~ temperature ~ at ~ .(object@dt_type) ~ (degree * C)),
               object))
  
  if (tolower(object@dt_name) == "sunshine")
    return(new("cfSunshine", 
               data_label = paste(object@dt_type, "sunshine (hours)"),
               object))
  
  if (tolower(object@dt_name) == "radiation")
    return(new("cfSunshine", 
               data_label = paste(object@dt_type, "radiation (MJ/m2)"),
               object))
  
  if (tolower(object@dt_name) == "pressure")
    return(new("cfPressure",
               data_label = paste(object@dt_type, "atmospheric MSL pressure (hPa)"),
               object))
  
  return(new("cfOther",
             data_label = object@dt_name,
             object))
}

# Methods -----------------------------------------------------------------

#' @importFrom methods setClass
setMethod("show",
          signature(object = "cfData"),
          function (object)
          {
            cat(object@data_label, "\n")
            print(head(as(object, "data.frame"), 4))
            if (nrow(object) > 4){
              n_omitted = nrow(object) - 4
              cat(paste("[~~ omitted", n_omitted,"rows ~~]\n"))
            }
          }
)

#' @importFrom lubridate ymd_hm ymd_hms
#' @importFrom methods setAs
#' @importFrom stats setNames
setAs("cfData", "data.frame",
      function(from){
        df_names = gsub(")", "", from@names, fixed = TRUE)
        df_names = gsub("(", ".", df_names, fixed = TRUE)
        df_names = gsub("/", "", df_names, fixed = TRUE)
        cf_df = setNames(data.frame(from, stringsAsFactors = FALSE),
                         df_names)
        cf_df[, 1] = factor(cf_df[, 1])
        if (grepl("rain rate", tolower(from@dt_type)))
          cf_df[, 2] = ymd_hms(cf_df[, 2], tz = "NZST")
        else
          cf_df[, 2] = ymd_hm(cf_df[, 2], tz = "Pacific/Auckland")
        cf_df
      })

#' @importFrom methods setMethod
#' @importFrom utils head
setMethod("head",
          signature(x = "cfData"),
          function(x, n = 6L, ...) head(as(x, "data.frame"), n, ...))

#' @importFrom methods setMethod
#' @importFrom utils tail
setMethod("tail",
          signature(x = "cfData"),
          function(x, n = 6L, ...) tail(as(x, "data.frame"), n, ...))