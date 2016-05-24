##' Get basic data on all counties
##' 
##' This function is essentially a  wrapper around Local.getCounties().
##' @usage getAllCounties()
##' @return A data frame with a row for each county and columns with the following variables describing the county:\cr counties.county*.localId,\cr counties.county*.name,\cr counties.county*.url.
##' @references http://api.votesmart.org/docs/Local.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of all counties 
##' \dontrun{counties <- getAllCounties()}
##' \dontrun{head(counties)}
##' @export



getAllCounties <-
  function() {
    
    states <- State.getStateIDs()
    nocounties <- c("NA", "PR", "GU", "AS", "VI", "DC")
    states <- states[!(states$stateId %in% nocounties),]
    Local.getCounties(states$stateId) # should be extended with Rwebapi-dl functions! 
  }