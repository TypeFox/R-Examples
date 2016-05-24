#' Gets TVKs for a group
#' 
#' Given the name of a group (see \code{\link{listGroups}}) this function returns the pTVKs 
#' (preferred taxon version keys) for all members of that group. This is currently restricted
#' to returning up to 5000 results.
#' 
#' @export
#' @param name A string for a group on the NBN gateway (e.g. 'reptile')
#' @return A vector of TVKs as characters.
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @seealso \code{\link{getFeature}}, \code{\link{getOccurrences}}, \code{\link{getTVKQuery}}
#' @examples \dontrun{ 
#'  t <- getGroupSpeciesTVKs('reptile')
#' }
#' 
getGroupSpeciesTVKs<-function(name){
    
    groupID<-getGroupID(name)
    
    json <- runnbnurl(service="species", group=groupID)
    
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
        d <- as.data.frame(d)
        ## we are only interested in presences (not absences)
        if ("absence" %in% colnames(d)) {
            d <- d[which(d$absence == FALSE),]
        }
        return(as.character(d$ptaxonVersionKey))
    }
}