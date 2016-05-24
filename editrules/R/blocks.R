#' Decompose a matrix or edits into independent blocks
#'
#' \code{blocks} returns a \code{list} of independent blocks \eqn{M_i} such that
#' \eqn{M=M_1\oplus M_2\oplus\cdots\oplus M_n}. 
#'
#' @param M \code{matrix}, \code{\link{editmatrix}}, editarray or editset to be decomposed into independent blocks
#' @return list of independent subobjects of \code{M}.
#' @example ../examples/blocks.R
#' @export
blocks <- function(M){
    blocks <- blockIndex(contains(M))
    lapply( blocks, 
        function(b){
            reduce(M[b,,drop=FALSE])
        }
    )
}


#' \code{blockIndex} returns a list of row indices in a \code{logical} matrix \code{D} designating independent blocks.
#' 
#'  
#' @param D matrix of type \code{logical}
#' @return \code{list} of row indices in \code{D} indicating independent blocks.
#'      Empty rows (i.e. every column \code{FALSE}) are ignored.
#' 
#' @rdname blocks
#' @export
blockIndex <- function(D){

    block <- function(B){
        x1 <- FALSE
        x <- B[1,]
        while (sum(x1 != x)){
            x1 <- x
            b <- sapply( 1:nrow(B)
                  , function(i){
                    any(B[i,] & x)
                   }
                  )
            x <- colSums(B[b,,drop=FALSE]) > 0 #this is another way of "or"ring all found rows
        }
        b
    }
    orignames <- row.names(D)
    row.names(D) <- 1:nrow(D)
   
    #remove empty rows
    b <- rowSums(D) == 0
    D <- D[!b,,drop=FALSE]
   
    # create a list which will contain the blocks
    blocks <- list()
    L <- 1
   
    # detect and remove blocks until no blocks are left
    while (nrow(D) > 0){
      
        # find block
        b <- block(D)
      
        # store the original row numbers of the detected block
        blocks[[L]] <- as.integer(row.names(D)[b])
        L <- L + 1
      
        # remove the detected block
        D <- D[!b,,drop=FALSE]
    }
    lapply(blocks,function(b) {names(b)<-orignames[b]; b})
}
