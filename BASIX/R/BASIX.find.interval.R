BASIX.find.interval <- function(vec,from,to,start=1){

ids <- .Call("find_windowC",vec,from,to,start,PACKAGE="BASIX")

return(ids)

}
