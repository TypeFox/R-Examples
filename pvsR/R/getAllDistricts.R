##' Get basic data on all districts
##' 
##' This function is essentially a wrapper around District.getByOfficeState().
##' @usage getAllDistricts()
##' @return A data frame with a row for each district and columns with the following variables describing the district:\cr districtList.district*.districtId,\cr districtList.district*.name,\cr districtList.district*.officeId,\cr districtList.district*.stateId.
##' @references http://api.votesmart.org/docs/District.html\cr
##' See http://api.votesmart.org/docs/semi-static.html for a list of office IDs or use Office.getOfficesByType(), Office.getOfficesByLevel(), Office.getOfficesByTypeLevel() or Office.getOfficesByBranchLevel().\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of all districts
##' \dontrun{districts <- getAllDistricts()}
##' \dontrun{head(districts)}
##' @export



getAllDistricts <-
  function() {
    
    of <- getOffices()
    states <- State.getStateIDs()
    
    District.getByOfficeState(of$officeId, states$stateId, all=TRUE) 
    
  }
