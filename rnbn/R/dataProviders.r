#' Get data providers
#' 
#' Gets details of data providers for a given dataset (or list of datasets)
#'
#' @export
#' @param datasets A string (or list of strings) giving the dataset keys to search within.
#' Use \code{\link{listDatasets}} to find dataset keys, or visit the NBN gateway
#' \url{https://data.nbn.org.uk/Datasets}
#' @return a data.frame giving the details of data providers 
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk} and Tom August, CEH \email{tomaug@@ceh.ac.uk}
#' @seealso \code{\link{listDatasets}}
#' @examples \dontrun{ 
#' 
#' providers <- dataProviders(c('GA000426','GA000832'))
#' 
#' }
#' 
dataProviders <- function(datasets=NULL) {
    
    # A function for removing html from strings
    cleanFun <- function(htmlString) {
        htmlString <- gsub("<.*?>", "", htmlString)
        htmlString <- gsub("\r","",htmlString)
        htmlString <- gsub(" \n",", ",htmlString)
    }
    
    if(is.null(datasets)) stop('datasets parameter cannot be NULL')
   
    org_master <- NULL
    
    #Set the columns that we want from each provider
    columnNames <- c("id","name","address","postcode","contactName","contactEmail","website")
        
    for(dataset in datasets){
        
        json <- runnbnurl(service="p", datasets=dataset)
        
        if (length(json) > 0 & class(json) == 'list'){
            
            # Get lead org data
            organisation <- as.data.frame(json$organisation)[,columnNames[columnNames %in% names(as.data.frame(json$organisation))]]
              
            # Get all other orgs data
            if(length(json$contributingOrganisations)!=0){
                COjson <- json$contributingOrganisations
                if (length(COjson) > 0) {
                    ## find the unique names that are used
                    n <- unique(unlist(c(sapply(COjson, names))))
                    ## dimension a matrix for the required number of rows and cols
                    d <- matrix(nrow=length(COjson), ncol=length(n), 
                                dimnames=list(seq(1:length(COjson)),n))
                    ## now we can go through the list and insert
                    ## the values into the correct cells of the matrix
                    ## This should be quick because the matrix is pre-allocated
                    for (i in 1:length(COjson)) {
                        for (j in 1:length(COjson[[i]])) {
                            k <- grep(names(COjson[[i]][j]),n)
                            d[i,k] <- COjson[[i]][[j]]
                        }
                    }
                    
                    ## cooerce the matrix to a data.frame
                    contributingOrganisations <- as.data.frame(d, stringsAsFactors = F)[,columnNames[columnNames %in% names(as.data.frame(json$organisation))]] 
                    # Combine the organisations into one dataframe
                    organisation <- rbind(organisation,contributingOrganisations)
                }
            }
        
        if(is.null(org_master)){org_master <- organisation} else{org_master <- merge(org_master, organisation, all = TRUE)}
    
        } else {
            warning(paste('Dataset',dataset,'returned no taxa, check this is spelt correctly'))            
        }
    }    
    
    # If we have some data do some tidying up
    if(!is.null(org_master)){
        # Get rid of bits of HTML
        org_master <- unique(as.data.frame(lapply(org_master,FUN=cleanFun)))
            
        # Sort the master df by name
        if(!is.null(org_master)){
            org_master <- org_master[with(org_master, order(name)), ]
            row.names(org_master) <- 1:nrow(org_master)    
        }
    }
    
    return(org_master)
}