#' Get group definitions
#' 
#' Gives a dataframe of the group definitions from the NBN Gateway for reference.
#' 
#' @export
#' @return A dataframe containing the definitions of groups on the NBN Gateway
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @seealso \code{\link{getFeature}}, \code{\link{getOccurrences}}
#' @examples \dontrun{ 
#'  t <- listGroups()
#' }

listGroups <- function() {
    
    ## return a JSON object (list of lists)
    json <- runnbnurl(service="list",list='taxonOutputGroups')
    
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
        
        # Order returned dataframe by name
        d<-d[with(d, order(d$name)),]
        
        return(d[c('name','key')])
        
    } else {
        return(NULL)
    }
}
