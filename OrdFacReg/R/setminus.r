setminus <- function(A, B){

    # compute A \ B
    res <- A[((A %in% B) == FALSE)]
    return(res)
}
