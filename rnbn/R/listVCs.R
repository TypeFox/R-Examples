#' Get list of vice-counties
#' 
#' Gives a dataframe of the Watsonian vice-counties and their keys for reference.
#' 
#' @export
#' @return A dataframe containing the name and featureID of each vice-county.
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @seealso \code{\link{getFeature}}, \code{\link{getOccurrences}}
#' @examples \dontrun{ 
#'  t <- listVCs()
#' }

listVCs <- function() {
    
    ## return a JSON object (list of lists)
    json <- runnbnurl(service="list",list='siteBoundaryDatasets/GA000344/siteBoundaries')
            
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
        return(d[c('name','identifier','featureID')])
    } else {
        return(NULL)
    }
}