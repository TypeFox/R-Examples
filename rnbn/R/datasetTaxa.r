#' Get a list of taxa present in a given dataset
#' 
#' Gets a list of taxa that are included in a given dataset or list of datasets.
#'
#' @export
#' @param datasets A string (or list of strings) giving the dataset keys to search within.
#' Use \code{\link{listDatasets}} to find dataset keys, or visit the NBN gateway
#' \url{https://data.nbn.org.uk/Datasets}
#' @return a data.frame of taxa details 
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk} and Tom August, CEH \email{tomaug@@ceh.ac.uk}
#' @seealso \code{\link{listDatasets}}
#' @examples \dontrun{ 
#' 
#' taxa <- datasetTaxa(datasets=c('GA001044','GA001011'))
#' 
#' }
#' 
datasetTaxa <- function(datasets=NULL) {
    
    if(is.null(datasets)) stop('datasets parameter cannot be NULL')
   
    d_master <- NULL
    
    for(dataset in datasets){
        
        ## return a JSON object (list of lists)
        json <- runnbnurl(service="d", datasets=dataset) 

        if (length(json) > 0) {
            ## find the unique names that are used
            n <- unique(unlist(c(sapply(json, names))))
            ## dimension a matrix for the required number of rows and cols
            d <- matrix(nrow=length(json), ncol=length(n), 
                        dimnames=list(seq(1:length(json)),n))
            ## now we can go through the list and insert
            ## the values into the correct cells of the matrix
            ## This should be quick because the matrix is pre-allocated
            for (i in 1:length(json)) {
                for (j in 1:length(json[[i]])) {
                    k <- grep(names(json[[i]][j]),n)
                    d[i,k] <- json[[i]][[j]]
                }
            }
            
            ## cooerce the matrix to a data.frame
            d <- as.data.frame(d, stringsAsFactors = F)    
            if(is.null(d_master)){d_master <- d} else{d_master <- merge(d_master, d, all = TRUE)}
            
        } else { # if returned json contains nothing
            warning(paste('Dataset',dataset,'returned no taxa, check this is spelt correctly'))            
        }
    }
    
    # Sort the master df by name
    if(!is.null(d_master)){
        d_master <- d_master[with(d_master, order(name)), ]
        row.names(d_master) <- 1:nrow(d_master)    
    }
        
    return(d_master)
}
