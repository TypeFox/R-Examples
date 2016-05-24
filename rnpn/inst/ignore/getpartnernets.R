#' Get a list of all the networks, i.e. organizations, with which NPN is associated.
#' @import RJSONIO RCurl plyr
#' @param url the PLoS API url for the function (should be left to default)
#' @param ... optional additional curl options (debugging tools mostly)
#' @param curl If using in a loop, call getCurlHandle() first and pass
#'  the returned value in here (avoids unnecessary footprint)
#' @return Data frame of each network’s name and id.
#' @export
#' @examples \dontrun{
#' getpartners()
#' }
getpartners <- 

function(url = 'http://www.usanpn.org/npn_portal/networks/getPartnerNetworks',
  ..., 
  curl = getCurlHandle() ) 
{
  url2 <- paste(url, '.json?', sep='')
  out <- fromJSON(getURLContent(url2))
  ldply(out, identity)
}