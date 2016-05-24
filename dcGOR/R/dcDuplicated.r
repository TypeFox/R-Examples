#' Function to determine the duplicated patterns from input data matrix
#'
#' \code{dcDuplicated} is supposed to determine the duplicated vectorised patterns from a matrix or data frame. The patterns can come from column-wise vectors or row-wise vectors. It returns an integer vector, in which the value indicates from which it duplicats.
#'
#' @param data an input data matrix/frame
#' @param pattern.wise a character specifying in which wise to define patterns from input data. It can be 'column' for column-wise vectors, and 'row' for row-wise vectors
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' an interger vector, in which an entry indicates from which it duplicats. When viewing column-wise patterns (or row-wise patterns), the returned integer vector has the same length as the column number (or the row number) of input data. 
#' @note
#' none
#' @export
#' @seealso \code{\link{dcAncestralMP}}, \code{\link{dcAncestralMP}}, \code{\link{dcAlgo}}
#' @include dcDuplicated.r
#' @examples
#' # an input data matrix storing discrete states for tips (in rows) X four characters (in columns)
#' data1 <- matrix(c(0,rep(1,3),rep(0,2)), ncol=1)
#' data2 <- matrix(c(rep(0,4),rep(1,2)), ncol=1)
#' data3 <- matrix(c(1,rep(0,3),rep(1,2)), ncol=1)
#' data <- cbind(data1, data2, data1, data3)
#' colnames(data) <- c("C1", "C2", "C3", "C4")
#' data
#' 
#' # determine the duplicated patterns from inut data matrix
#' res <- dcDuplicated(data, pattern.wise="column")
#' ## return an integer vector
#' res
#' ## get index for unique patterns
#' ind <- sort(unique(res))
#  data[,ind]
#' ## As seen above, the returned integer vector tells there are 3 unique patterns:
#' ## they are in columns (1, 2, 4). The column 3 is duplicated from column 1.

dcDuplicated <- function(data, pattern.wise=c("column","row"), verbose=T)
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    pattern.wise <- match.arg(pattern.wise)
    
    if(is.vector(data)){
        tmp_data <- matrix(data, ncol=1)
        if(!is.null(names(data))){
            rownames(tmp_data) <- names(data)
        }
        data <- tmp_data
    }
    if(is.data.frame(data)){
        data <- as.matrix(data)
    }
            
    if(pattern.wise=="row"){
        tdata <- t(data)
    }else if(pattern.wise=="column"){
        tdata <- data
    }
    
    colnames(tdata) <- 1:ncol(tdata)
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Merge %d patterns (%s) ...", ncol(tdata), as.character(now)), appendLF=T)
    }
    ## merge a column-wise vector into a string
    xx <- apply(tdata, 2, function(x){
        paste(x, collapse='')
    })
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Sort %d patterns (%s) ...", ncol(tdata), as.character(now)), appendLF=T)
    }
    ## sort merged strings
    xx_sort <- base::sort(xx)
    ## extract names (corresponding to column names of input data matrix)
    ind_sort <- as.numeric(names(xx_sort))

    if(verbose){
        now <- Sys.time()
        message(sprintf("Find index from %d patterns (%s) ...", ncol(tdata), as.character(now)), appendLF=T)
    }
    ## find index of those duplicated from merged strings
    x <- base::duplicated(xx_sort)
    ## create a vector to store duplicated index
    flag <- which(!x)
    tmp_all <- c(flag, length(x)+1)
    index_list <- lapply(1:(length(tmp_all)-1), function(i){
        rep(tmp_all[i], tmp_all[i+1]-tmp_all[i])
    })
    index_vec <- unlist(index_list)
    ## create a matrix storing the from-to mapping
    df <- data.frame(from=ind_sort, to=ind_sort[index_vec])
        
    ## return an index vector of the same length equaling the column number of 'tdata', each index indicating which column has the same pattern as it has
    res <- df[order(df$from),]$to
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("For the input data (with %d patterns), there are %d unique '%s'-wise patterns (%s)", ncol(tdata), length(unique(res)), pattern.wise, as.character(now)), appendLF=T)
    }
    
    return(invisible(res))
}
