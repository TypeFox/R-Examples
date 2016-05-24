#' Get the taxonomical heirarchy for a given taxon
#' 
#' Given the Taxon Version Key (a 16-character string), this function get details
#' of the taxonomical heirarchy above the given TVK. This gives the rank, name and 
#' preferred Taxon Version Key for each level of the taxonomy
#' 
#' @export
#' @param tvks A Taxon Version Key (TVK) which is a 16-character string ID
#' @return A dataframe containing the JSON object returned by the NBN Gateway.
#' @author Tom August, CEH \email{tom.august@@ceh.ac.uk}
#' @seealso \code{\link{getFeature}}, \code{\link{getOccurrences}}
#' @examples \dontrun{ 
#'  t <- getTaxonomy("NHMSYS0000528028") # Myotis daubentonii (Daubenton's Bat)
#' }

getTaxonomy <- function(tvks=NULL) {
    
    ## return a JSON object (list of lists)
    json <- runnbnurl(service="ancestry", tvks=tvks) 
    
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
        return(d[c('rank','name','ptaxonVersionKey')])
    } else {
        return(NULL)
    }
}

