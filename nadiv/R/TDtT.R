TDtT <- function(A, ...){
    ch <- chol(A)
    dd <- diag(ch)
 return(list(T = t(drop0(zapsmall(ch / dd, ...))), D = Diagonal(nrow(A), dd^2)))
}

