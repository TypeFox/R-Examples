#' This function retrieves time series data.
#'
#' @author Claudia Vitolo
#'
#' @description Given the station identification number(s), this function retrieves data (time series in zoo format with accompanying metadata) from the WaterML2 service on the NRFA database. The time series can be of two types: \code{cmr} (catchment mean rainfall, monthly) or \code{gdf} (gauged daily flows, daily).
#'
#' @param id station identification number(s), each number should be in the range [3002,236051].
#' @param type This is character string that can have one of the two following values: "cmr" (to obtain catchment mean rainfall) or "gdf" (to obtain gauged daily flow).
#' @param metadata Logical, FALSE by default. If metadata = TRUE means that the result for a single station is a list with two elements: data (the time series) and meta (metadata).
#' @param parallel Logical, FALSE by default. If parallel = TRUE means that the function can be used in parallel computations.
#'
#' @return list composed of as many objects as in the list of station identification numbers. Each object can be accessed using their names or index (e.g. x[[1]], x[[2]], and so forth). Each object contains a zoo time series.
#'
#' @examples
#' # getTS(18019, type = "cmr")
#' # getTS(c(54022,54090,54091), type = "cmr")
#' # getTS(18019, type = "gdf")
#' # getTS(c(54022,54090,54091), type = "gdf")
#'

getTS <- function(id, type, metadata = FALSE, parallel = FALSE){

  # require(RCurl)
  # require(XML2R)
  # require(stringr)
  # require(zoo)
  # id <- c(54022,54090,54091)

  options(warn=-1)                                       # do not print warnings

  id <- as.character(id)       # in case it is a factor, convert to characters

  if (length(as.list(id)) == 0) {

    message("Please, enter valid id.")
    stop

  }else{

    if (length(as.list(id)) > 1 & parallel == FALSE){

      # multiple identification numbers
      tsList <- lapply(X = as.list(id), FUN = getTS_internal, type, metadata)
      names(tsList) <- id

    }else{

      # this is the case of a single identification number
      if (metadata) {
        tsList <- getTS_internal(id, type, metadata)
      }else{
        tsList <- unlist(getTS_internal(id, type, metadata))
      }

    }

  }

  return(tsList)

}


getTS_internal <- function(idx, type, metadata){

  website <- "http://nrfaapps.ceh.ac.uk/nrfa"

  myURL <- paste(website,"/xml/waterml2?db=nrfa_public&stn=",
                 idx, "&dt=", type, sep="")

  if ( url.exists(myURL) ){

    doc <- urlsToDocs(myURL)
    nodes <- docsToNodes(doc,xpath="/")
    myList <- nodesToList(nodes)

    data <- FindTS(myList)
    if (metadata) meta <- FindInfo(myList)

  }else{

    message(paste("For station", idx,
                  "there is no available online dataset in waterml format \n"))

    data <- NULL
    if (metadata) meta <- NULL

  }

  if (metadata) {
    return(list("data" = data, "meta" = meta))
  }else{
    return(data)
  }

}
