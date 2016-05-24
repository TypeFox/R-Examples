#' Build the URL to call NBN web-services
#' 
#' The base URL for the RESTful services is:
#' \url{https://data.nbn.org.uk/api/} 
#' Parameters can be added, depending on the service called. \cr\cr
#' The following services are available at the time of writing: \cr
#' \code{taxonObservations?ptvk=<tvk>&datasetKey=<dataset>&startYear=<year> 
#' &endYear=<year>&featureID=<feature>} - get observations for the given
#' (list of) Taxon Version Keys (TVKs) for the given (list of) dataset keys
#' optionally with dates between startYear and endYear and within the spatial
#' feature, featureID (e.g. a vice-county). If no dataset keys are given then all 
#' publically available data for the taxa will be returned \cr\cr
#' \code{features/<featureID>} - get details about a particular location given
#' its <featureID> \cr\cr
#' \code{taxa/<tvk>} - get details about a taxon given its <tvk>
#' 
#' @export
#' @param service the service you want to call. One of \code{"obs"} for the 
#'   taxonObservations service, \code{"feature"} for the features service, 
#'   \code{"taxon"} for the taxa service, \code{"list"} for listing services,
#'   \code{"ancestry"} for taxonomy service, \code{"species"} for the 
#'   species service and \code{"query"} for the search service.
#' @param tvks a list of TVKs which are strings of 16 alphanumeric characters
#' @param datasets a list of dataset keys which are strings of 8 alphanumeric 
#'   characters
#' @param feature a featureID an integer number
#' @param startYear a 4 digit year
#' @param endYear a 4 digit year
#' @param list url segment as a string as required to append to the base url to 
#'        give the list required as a part of the \code{"list"} service 
#' @param VC a string giving a vice-county name (see \code{\link{listVCs}})
#' @param group a string giving the name of a group (see \code{\link{listGroups}})
#' @param query a string used to search for taxa
#' @param gridRef a string giving a gridreference in which to search for occurrences
#' @param attributes if \code{TRUE} then attribute data is returned
#' @return the URL to call - a character string
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @examples
#' makenbnurl(service="obs", tvks="NBNSYS0000007073")
#' makenbnurl(service="obs", tvks="NBNSYS0000007073", datasets="SGB00001")
#' makenbnurl(service="obs", tvks="NBNSYS0000007073", datasets="SGB00001",
#'            startYear="1990", endYear="2010")
#' makenbnurl(service="feature", feature="284443")
#' makenbnurl(service="taxon", tvks="NBNSYS0000007073")
#' 
makenbnurl <- function(service=NULL, tvks=NULL, datasets=NULL, feature=NULL,
                       startYear=NULL, endYear=NULL, list=NULL, VC=NULL, group=NULL,
                       query=NULL, gridRef=NULL, attributes=FALSE) {

    ##----------------------------------------------------------------------
    ## function to check that parameters are correctly formatted
    ## id     - the parameter to check
    ## list   - is it allowed to be a list, or should it be a single value?
    ## len    - required length - number of characters
    ## num    - should it be a number (chars 0-9 only)?
    ## Returns: TRUE/FALSE
    checkID <- function(id, list=FALSE, len=0, num=FALSE) {
        if (!list & length(id) > 1) {
            ok <- FALSE
        } else {
            if (num) {
                if (len > 0) {
                    pattern <- paste("^[0-9]", "{", len, ",", len, "}$", sep="")    
                } else {
                    pattern <- "^[0-9]+$"
                }
            } else {
                pattern <- paste("^[0-9,A-Z,a-z]", "{", len, ",", len, "}$", sep="")
            }
            ok = (length(grep(pattern, id)) == length(id))
        }
        ok
    }
    ##----------------------------------------------------------------------
    
    ## This is the base url for the NBN web-services as of April 2013
    baseURL <- "https://data.nbn.org.uk/api/"
    url <- baseURL
    
    ## cope with year/feature parameters being given as numbers
    if (!is.null(startYear)) startYear <- as.character(startYear)
    if (!is.null(endYear)) endYear <- as.character(endYear)
    if (!is.null(feature)) feature <- as.character(feature)
    
    ## which service?
    if (is.character(service)) {
        ## get single lower case character from start of parameter. Should be
        ## o (taxon Occurrences), f (feature details) or t (taxon details)
        svc <- tolower(substr(as.character(service),1,1))
        
        switch(svc, 
            ## taxon occurrences ---------------------------------------------
            ## MUST have tvks (which can be a list)
            ## datasets (can be a list), startYear and endYear are optional
            o ={
                url <- paste(url, "taxonObservations?", sep="")     
                if (is.null(gridRef) & is.null(tvks)) stop("One of tvks or gridRef is required")
                if (is.character(tvks)) {
                    if (checkID(tvks, list=TRUE, len=16)) {
                        url <- paste(url, "ptvk=", paste(unlist(tvks), collapse="&ptvk="), sep="")
                    } else {
                        stop("tvks parameter is incorrect")
                    }
                } else {
                    #stop("tvks parameter is required")
                }
                if (is.character(datasets)) {
                    if (checkID(datasets, list=TRUE, len=8)) {
                        url <- paste(url, "&datasetKey=", paste(datasets, collapse="&datasetKey="), sep="")
                    } else {
                        stop("datasets parameter is incorrect")    
                    }
                }
                if (is.character(startYear)) {
                    if (checkID(startYear, num=TRUE, len=4)) {
                        url <- paste(url, "&startYear=", startYear, sep="")    
                    } else {
                        stop("startYear parameter is incorrect")
                    }
                }
                if (is.character(endYear)) {
                    if (checkID(endYear, num=TRUE, len=4)) {
                        if ((is.character(startYear)) && (as.numeric(startYear) > as.numeric(endYear))) {
                                stop("startYear cannot be later than endYear")
                        }
                        url <- paste(url, "&endYear=", endYear, sep="") 
                    } else {
                        stop("endYear parameter is incorrect")
                    }
                }
                if (!is.null(VC)) {
                    if(is.character(VC)) VCID <- getVCid(VC)
                    url <- paste(url, "&featureID=", VCID, sep="") 
                }
                if (!is.null(gridRef)) {
                    url <- paste(url, "&gridRef=", gridRef, sep="") 
                }
                if (attributes) {
                    url <- paste(url, "&includeAttributes=true", sep="") 
                }
            },
               
            ## details for a feature -------------------------------------
            ## MUST have a single value in feature
            f = {
                url <- paste(url, "features/", sep="")
                if (is.character(feature)) {
                    if (checkID(feature, list=FALSE, num=TRUE)) {
                        url <- paste(url, feature, sep="")
                    } else {
                        stop("feature parameter is incorrect")
                    }
                } else {
                    stop("feature parameter is required")
                }
            },
               
            ## details for a taxon ---------------------------------------
            ## MUST have single value in tvks
            t = {
                url <- paste(url, "taxa/", sep="")
                if (is.character(tvks)) {
                    if (checkID(tvks, list=FALSE, len=16)) {
                        url <- paste(url, tvks, sep="")    
                    } else {
                        stop("tvks parameter is incorrect")    
                    }
                } else {
                    stop("tvks parameter is required")
                }
            },
               
           ## details for ancestry (taxonomy) ---------------------------------------
           ## MUST have single value in tvks
           a = {
               url <- paste(url, "taxa/", sep="")
               if (is.character(tvks)) {
                   if (checkID(tvks, list=FALSE, len=16)) {
                       url <- paste(url, paste(tvks,'/taxonomy',sep=''), sep='')    
                   } else {
                       stop("tvks parameter is incorrect")    
                   }
               } else {
                   stop("tvks parameter is required")
               }
           },
               
           ## details for reference lists ---------------------------------------
           l = {
               url <- paste(url, list, sep="")
           },
           
           ## details of species in a given group -------------------------------
           s = {
               url <- paste(url, 'taxa?&taxonOutputGroupKey=', group, '&rows=5000', sep="")
           },
               
           ## details of taxa matching a search -------------------------------
           q = {
               url <- paste(url, 'search/taxa?q=', query, sep="")
           },
           
           ## details of taxa in a dataset -------------------------------
           d = {
               url <- paste(url, 'taxonDatasets/', datasets, '/taxa', sep="")
           },
           ## details of providers of a dataset -------------------------------
           p = {
               url <- paste(url, 'datasets/', datasets, sep="")
           },
        stop("service not recognised")) ## end of switch
        
    ## no value given for service
    } else {
        stop("no service specified")    
    }
    
    return(url)  ## return value
}
