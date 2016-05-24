#
# Export a bds matrix in "list mode".
#  This has one row for each non-zero element
# Input: a bdsmatrix
# Output: a data frame containing "row", "col", "value" as variables
#
# Options:
#    id:  True: row/col contain the subject id (dimnames of the matrix)
#         False:row/col contain integers
#    diag: True -- the output contains the diagonal of the matrix
#          False-- the output does not contain the diagonal
#
listbdsmatrix <- function(x, id=TRUE, diag=FALSE) {
    if (!inherits(x, 'bdsmatrix')) stop("Invalid argument")
    
    nblock <- length(x@blocksize)
    bsize  <- length(x@blocks)
    indx <- .C("bdsmatrix_index2",
               as.integer(nblock),
               as.integer(x@blocksize),
               rows= integer(bsize),
               cols= integer(bsize))
    
    # toss any zeros, and optionally the diagonal
    if (diag) toss <- (x@blocks==0)
    else      toss <- (x@blocks==0 | indx$rows== indx$cols)
    
    dd <- dimnames(x)[[1]]
    if (id && !is.null(dd)) {
        xr <- dd[indx$rows]
        xc <- dd[indx$cols]
        }
    else {
        xr <- indx$rows
        xc <- indx$cols
        }

    if (any(toss)) 
         data.frame(row=xr[!toss], col=xc[!toss], value=x@blocks[!toss])
    else data.frame(row=xr, col=xc, value=x@blocks)
    }
