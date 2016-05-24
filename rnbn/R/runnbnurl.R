#' Run the constructed url and get data
#' 
#' The URL to call a web-service is constructed using makenbnurl()
#' This is run and, hopefully, returns a JSON object
#' 
#' @export
#' @import httr
#' @param service the service you want to call. One of \code{"obs"} for the 
#'   taxonObservations service, \code{"feature"} for the features service, 
#'   \code{"taxon"} for the taxa service, \code{"list"} for listing services,
#'   \code{"ancestry"} for taxonomy service and \code{"species"} for the 
#'   species service. The first letter is sufficient
#' @param tvks a list of TVKs which are strings of 16 alphanumeric characters
#' @param datasets a list of dataset keys which are strings of 8 alphanumeric 
#'   characters
#' @param feature a featureID which is a string of 8 alphanumeric characters
#' @param startYear a 4 digit integer year
#' @param endYear a 4 digit integer year
#' @param list url segment as a string as required to append to the base url to 
#'        give the list required as a part of the \code{"list"} service.
#' @param VC a string giving a vice-county name (see \code{\link{listVCs}})
#' @param group a string giving the name of a group (see \code{\link{listGroups}})
#' @param query a string used to perform a taxa search
#' @param gridRef a string giving a gridreference in which to search for occurrences
#' @param attributes if \code{TRUE} then attribute data is returned
#' @return a JSON object resulting from the call
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @examples \dontrun{ 
#'  json <- runnbnurl(service="obs", tvks="NBNSYS0000007073", datasets="SGB00001",
#'                    startYear="1990", endYear="2010")
#' }
#' 
runnbnurl <- function(service=NULL, tvks=NULL, datasets=NULL, feature=NULL,
                      startYear=NULL, endYear=NULL, list=NULL, VC=NULL, group=NULL,
                      query=NULL, gridRef=NULL, attributes=FALSE) {
    
    url <- makenbnurl(service=service, tvks=tvks, datasets=datasets, feature=feature,
                      startYear=startYear, endYear=endYear, list=list, VC=VC,
                      group=group, query=query, gridRef=gridRef, attributes=attributes)
    
    
    # Run login script, this checks whether the user has a cookie that the webservice 
    # knows and if not gets them to log in and stores the cookies
    nbnLogin()
    
    #print(url)
    
    resp_url <- GET(url)
    
    if(!grepl("^2[[:digit:]]{2}$", resp_url$status_code)) stop(paste('Error accessing', url, '-', http_status(resp_url)$message))
        
    return(content(resp_url))
}
