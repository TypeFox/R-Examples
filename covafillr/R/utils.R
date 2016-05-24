
is.numsca <- function(x) is.numeric(x) & (length(x) == 1)
is.numint <- function(x) is.integer(x) & (length(x) == 1)
is.numvec <- function(x) is.vector(x) & is.numeric(x)
is.nummat <- function(x) is.matrix(x) & is.numeric(x)
is.samelength <- function(x,y){
    if(is.nummat(x) & is.nummat(y))
        return(all(dim(x)==dim(y)))
    if(is.numvec(x) & is.nummat(y))
        return(length(x) == dim(y)[1])
    if(is.nummat(x) & is.numvec(y))
        return(dim(x)[1] == length(y))
    if(is.numvec(x) & is.numvec(y))
        return(length(x) == length(y))
    return(FALSE)
}
