#' Function to create a sparse matrix for an input file with three columns
#'
#' \code{dcSparseMatrix} is supposed to create a sparse matrix for an input file with three columns.
#'
#' @param input.file an input file containing three columns: 1st column for rows, 2nd for columns, and 3rd for numeric values. Alternatively, the input.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a list containing arguments and their default values
#' @note
#' None
#' @export
#' @seealso \code{\link{dcAlgoPredictMain}}
#' @include dcSparseMatrix.r
#' @examples
#' # create a sparse matrix of 4 X 2
#' input.file <- rbind(c('R1','C1',1), c('R2','C1',1), c('R2','C2',1), c('R3','C2',1), c('R4','C1',1))
#' res <- dcSparseMatrix(input.file)
#' res
#' # get a full matrix
#' as.matrix(res)

dcSparseMatrix <- function(input.file, verbose=T)
{
    
    if(is.matrix(input.file) | is.data.frame(input.file)){
        if(verbose){
            message(sprintf("Load the input data matrix of %d X %d ...", dim(input.file)[1], dim(input.file)[2]), appendLF=T)
        }
        
        if(ncol(input.file)==2){
            input.file <- cbind(input.file, rep(1,nrow(input.file)))
        }
        
        if(is.data.frame(input.file)){
            input <- cbind(as.character(input.file[,1]), as.character(input.file[,2]), as.character(input.file[,3]))
        }else{
            input <- input.file
        }
    }else if(is.character(input.file) & input.file!=''){
        if(verbose){
            message(sprintf("Read the input file '%s' ...", input.file), appendLF=T)
        }       
        input <- utils::read.delim(input.file, header=T, sep="\t", colClasses="character")
        if(ncol(input)==2){
            input <- cbind(input, rep(1,nrow(input)))
        }
    }else{
        return(NULL)
    }
    
    x <- input
    if(!is.null(x)){
        x_row <- sort(unique(x[,1]))
        x_col <- sort(unique(x[,2]))
        ind_row <- match(x[,1], x_row)
        ind_col <- match(x[,2], x_col)
        x.sparse <- Matrix::sparseMatrix(i=ind_row, j=ind_col, x=as.numeric(x[,3]), dims=c(length(x_row),length(x_col)))
        rownames(x.sparse) <- x_row
        colnames(x.sparse) <- x_col
        
        if(verbose){
            message(sprintf("There are %d entries, converted into a sparse matrix of %d X %d.", nrow(x), dim(x.sparse)[1], dim(x.sparse)[2]), appendLF=T)
        }
    }else{
        return(NULL)
    }
    
    invisible(x.sparse)
}
