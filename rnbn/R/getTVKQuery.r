#' Gets TVKs for a query
#' 
#' Given a search term this function returns taxon information, including pTVKs,
#' for the first 25 taxa that match that search on the NBN. 
#' 
#' @export
#' @param query A query string. This can range from latin binomials to partial english names.
#' @param species_only Logical, if \code{TRUE} pTVKs of species are returned (i.e.
#' sub-species and aggregates are removed). Defaults to \code{TRUE}
#' @param rec_only Logical, if \code{TRUE} pTVKs of recommended names are returned (i.e.
#' synonyms are removed). Defaults to \code{FALSE}. Remember, the pTVK of a synonym is a 
#' taxa with 'recommended' name status.
#' @param top Logical, if \code{TRUE} only the top answer is returned. This is what the
#' gateway thinks you are most likely to be after but may not always be right, use with
#' care!
#' @return A dataframe containing information on each taxa entry that matches the query 
#' string in rows. ptaxonVersionKey (preferred taxon version key) should be used when
#' searching for records.
#' @author Tom August, CEH \email{tomaug@@ceh.ac.uk}
#' @seealso \code{\link{getGroupSpeciesTVKs}}, \code{\link{getOccurrences}}
#' @examples \dontrun{ 
#'  t <- getTVKQuery('blue tit')
#' }
#' 
getTVKQuery<-function(query=NULL, species_only = TRUE, rec_only = FALSE, top = FALSE){
    
    if(is.null(q)) stop('query string must not be null')
    
    d_master <- NULL
    
    for(q in query){
        
        q <- tolower(gsub(' ','%20', q))
        
        json <- runnbnurl(service="query", query=q)
        
        json<-json$results
        
        if (length(json) > 0) {
            ## find the unique names that are used in occ
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
            d <- as.data.frame(d, stringsAsFactors = FALSE)
            
            # On reflection it is probably best to return everything # keep only the columns I need
            #d <- d[colnames(d) %in% c('searchMatchTitle','rank','nameStatus','ptaxonVersionKey')]            
           
            ## Keep only species if desired
            if(species_only) d <- d[d$rank == 'Species',]
            
            ## Keep only recommended if desired
            if(rec_only) d <- d[d$nameStatus == 'Recommended',]
            
            ## Keep top only
            if(top) d <- d[1,]
            
            if(is.null(d_master)){d_master <- d}else{d_master<-merge(d_master,d,all=T)}
        }
    }
    
    return(d_master)
    
}