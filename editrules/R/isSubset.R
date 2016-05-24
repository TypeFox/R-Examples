
#' Check which edits are dominated by other ones.
#'
#' An edit defines a subregion of the space of all possible value combinations 
#' of a record. Records in this region are interpreted as invalid. An edit rule
#' which defines a region equal to or contained in the region defined by another
#' edit is redundant. (In data editing literature, this is often referred to as 
#' a \emph{domination} relation.)
#'
#' @param E \code{\link{editarray}} 
#' @return \code{logical} vector indicating if an edit is a subset of at least one other edit.
#' @export
isSubset <- function(E){
    if ( !is.editarray(E) ) stop('argument is not an editarray')
    isSubset.boolmat(getArr(E))
}

isSubset.boolmat <- function(A){
    if (nrow(A)==0) return(logical(0))
    if (nrow(A)==1) return(FALSE)

    d <- duplicated(A)
    wd <- which(d)
    m <- nrow(A)
    if ( m == 0 ) return(logical(0))
    M <- (1:m)[!d]
    if ( length(M) == 1 ) return(d)
    m1 <- length(M)-1
    s <- logical(m)
    s[M] <- vapply(M, 
        function(i){
            I <- c(i,wd)
            any(rowSums(A[-I,,drop=FALSE] - (A[-I,,drop=FALSE] | A[rep(i,m1),,drop=FALSE]) ) == 0)
        },
        FUN.VALUE=FALSE 
    )
    s | d
}

# check if edits in A are subset of edits in B: (returns boolean vector)
isSubsetWrt.boolmat <- function(A,B){
    m <- nrow(A)
    n <- nrow(B)
    if ( m == 0 ) return(logical(0))
    if ( n == 0 ) return(rep(FALSE,n))

    vapply(1:m,function(i){
        any(rowSums(abs({A[rep(i,n),,drop=FALSE] | B} - B)) == 0)
    },FUN.VALUE=FALSE)
}





