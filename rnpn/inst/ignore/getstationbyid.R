#   FOR SOME REASON THIS DOESN'T WORK, PROBLEM ON SERVER SIDE I THINK
#   FOR SOME REASON THIS DOESN'T WORK, PROBLEM ON SERVER SIDE I THINK
#   FOR SOME REASON THIS DOESN'T WORK, PROBLEM ON SERVER SIDE I THINK
#   FOR SOME REASON THIS DOESN'T WORK, PROBLEM ON SERVER SIDE I THINK
#   FOR SOME REASON THIS DOESN'T WORK, PROBLEM ON SERVER SIDE I THINK
#   FOR SOME REASON THIS DOESN'T WORK, PROBLEM ON SERVER SIDE I THINK


#' Get a list of all the networks, i.e. organizations, with which NPN is associated.
#' @import RJSONIO RCurl plyr
#' @param stationid Station ID (numeric).
#' @param downform Download format, one of 'json' or 'xml'.
#' @param url the PLoS API url for the function (should be left to default)
#' @param ... optional additional curl options (debugging tools mostly)
#' @param curl If using in a loop, call getCurlHandle() first and pass
#'  the returned value in here (avoids unnecessary footprint)
#' @return Data frame of each network’s name and id.
#' @export
#' @examples \dontrun{
#' getstationbyid(5122)
#' getstationbyid(c(5122,1915))
#' getstationbyid(5122, xml)
#' }
# getstationbyid <- 
# 
# function(stationid = NA, downform = 'json',
#   url = 'http://www.usanpn.org/npn_portal/stations/getStationById',
#   ..., 
#   curl = getCurlHandle() ) 
# {
#   url2 <- paste(url, '.', downform, sep='')
#   args <- list()
#   if(!is.na(stationid[1]))
#     for(i in 1:length(stationid)) {
#       args[paste('station_ids[',i,']',sep='')] <- stationid[i]
#     }
#   tt <- getForm(url2,
#           .params = args,
#           ...,
#           curl = curl)
#   fromJSON(tt)
# }