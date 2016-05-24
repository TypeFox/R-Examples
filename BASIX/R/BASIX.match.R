BASIX.match <- function(elements, vec){

ids <- .Call("my_match_C",elements,vec, PACKAGE="BASIX")

#ids[ids==-1] <- NA

return(ids)

}
