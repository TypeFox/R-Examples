#' Get occurrences for a given species
#' 
#' Gets occurrence data from the NBN to which you have access. To get access to data
#' you must first register at \url{https://data.nbn.org.uk/User/Register}. You will
#' need your username and password when running this function for the first time.\cr
#' You can specify the data to retrieve by dataset, species, time, location and/or
#' group.
#'
#' @export
#' @param tvks a list of TVKs which are strings of 16 alphanumeric characters.
#' You can look these up using \code{getTVKQuery}
#' @param datasets a list of dataset keys which are strings of 8 alphanumeric 
#'   characters. Look up datasets here: \url{https://data.nbn.org.uk/Datasets}
#' @param startYear a 4 digit integer year
#' @param endYear a 4 digit integer year
#' @param VC a string giving a vice-county name (see \code{\link{listVCs}})
#' @param group a string giving the name of a group (see \code{\link{listGroups}}).
#' Using group will retireve data for all TVKs in this group. for example using group 'reptile'
#' will search using over 150 TVKs including TVKs for higher taxonomic groups such families
#' within reptiles. Therefore it may be preferrable to search using a list TVKs aquired
#' using getTVKQuery
#' @param gridRef a string giving a gridreference in which to search for occurrences
#' @param latLong logical, if TRUE latitude and longitude are returned as additional columns.
#' The conversion to latitude and longitude is currently accurate to about about ~20 meters,
#' greater than the vast majoring of records' precision.
#' @param acceptTandC if set to \code{TRUE} you accept the NBN gateway terms and 
#' conditions and privacy policy. These can be found at \url{https://data.nbn.org.uk/Terms}.
#' Accepting the terms and conditions supresses the corresponding warning message.
#' @param silent If \code{TRUE} batch request information is supressed
#' @param attributes If \code{FALSE} then attribute data is not returned, this may
#' improve the speed of large requests.
#' @return a data.frame of occurence records. Details of the data providers that 
#' contributed to the data returned is given as a 'providers' attribute
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk} and Tom August, CEH \email{tomaug@@ceh.ac.uk}
#' @seealso \code{\link{getFeature}}, \code{\link{getTVKQuery}}, \code{\link{listVCs}},
#' \code{\link{listDatasets}}, \code{\link{listGroups}}
#' @examples \dontrun{ 
#'  dt1 <- getOccurrences(tvks="NBNSYS0000002987", datasets="GA000373", 
#'                        startYear="1990", endYear="1991")
#'                       
#'  dt2 <- getOccurrences(tvks=c("NBNSYS0000002987","NHMSYS0001688296","NHMSYS0000080210"),
#'                        startYear="1990", endYear="1991")
#'                        
#'  dt3 <- getOccurrences(group="quillwort", startYear="1990", endYear="2010",
#'                        VC="Shetland (Zetland)")
#'  
#'  # Get the data providers information
#'  dp <- attr(dt1,'providers')                      
#'                        
#' }
#' 
getOccurrences <- function(tvks=NULL, datasets=NULL, startYear=NULL, 
                           endYear=NULL, VC=NULL, group=NULL, gridRef=NULL,
                           latLong = TRUE, acceptTandC=FALSE, silent=FALSE,
                           attributes = FALSE) {
    
    if(!is.null(tvks) & !is.null(group)) stop('group and tvks cannot be used at the same time')
    if(is.null(tvks) & is.null(group) & is.null(gridRef)) stop('One of group, tvks or gridRef must be given')
    
    # If we are searching by group get the group tvks
    if(!is.null(group)) tvks <- getGroupSpeciesTVKs(group)
    
    ## If you have more than 5 TVKs break it up into batches of 2
    # Set up parameters
    if(!is.null(tvks)){
        tvks <- unique(tvks)
        nTVK <- length(tvks)
    } else {
        tvks <- 1
        nTVK <- 1
    }
    
    start <- 1
    d_master <- NULL
    
    while(start <= nTVK){
        
        if(!silent) cat('Requesting batch', ceiling(start/2), 'of', ceiling(nTVK/2),'\n', sep=' ')
        end <- start + 1
        if(!is.null(tvks)){temp_tvks <-  na.omit(tvks[start:end])}else{temp_tvks=NULL}
        
        ## return a JSON object (list of lists)
        json <- runnbnurl(service="obs", tvks=temp_tvks, datasets=datasets, 
                          startYear=startYear, endYear=endYear, VC=VC,
                          gridRef=gridRef, attributes=attributes) 
        
        if (length(json) > 0) {
            ## find the unique names that are used in occ
            n <- unique(unlist(c(sapply(json, function(x) names(unlist(x))))))
            ## dimension a matrix for the required number of rows and cols
            d <- matrix(nrow = length(json), ncol = length(n), 
                        dimnames = list(seq(1:length(json)), n))
            ## now we can go through the list and insert
            ## the values into the correct cells of the matrix
            ## This should be quick because the matrix is pre-allocated
            ## The unlisting allows us to bring out the attributes fields
            
            for(i in 1:length(json)) {
                for (j in 1:length(unlist(json[[i]]))) {
                    k <- grep(names(unlist(json[[i]])[j]),n)
                    d[i,k] <- unlist(json[[i]])[[j]]
                }
            }
            
            ## cooerce the matrix to a data.frame
            d <- as.data.frame(d, stringsAsFactors = F)
            
            if(is.null(d_master)){d_master <- d} else{d_master <- merge(d_master, d, all = TRUE)}
            
        }

        start <- start + 2
        
    }#end of while
    
    if ("absence" %in% colnames(d_master)) {
        if(TRUE %in% d_master$absence) warning('NOTE: There are absence records in this dataset')
    }
    
    ## Format date columns as dates
    if("startDate" %in% colnames(d_master)) d_master$startDate <- as.Date(d_master$startDate)
    if("endDate" %in% colnames(d_master)) d_master$endDate <- as.Date(d_master$endDate)
    
    ## Add lat long if requested
    if(latLong & !is.null(d_master)){ 
        latlong <- gr2gps_latlon(d_master$location, centre=TRUE)
        d_master$latitude <- latlong$LATITUDE
        d_master$longitude <- latlong$LONGITUDE        
    }
    
    ## Add an attribute giving details of the data providers
    if(!is.null(d_master)){
        if(!silent) cat("Requesting data providers' information\n")
        datasets <- unique(d_master$datasetKey)
        providers <- dataProviders(datasets)
        attr(x=d_master,which='providers') <- providers
    }
        
    ## Write out a statement about the T's & C's
    if(!acceptTandC) message('IMPORTANT: By using this package you are agreeing to the Gateway Terms & Conditions and Privacy Policy (see https://data.nbn.org.uk/Terms). This message can be suppressed using the acceptTandC argument') 
    
    return(d_master)
}