rad <- function(x){
    if(is(x, "fitsad"))
        y <- new("rad", data.frame(rank=1:length(x@data$x), abund=sort(x@data$x, decreasing=T)))
    else if(is(x,"fitrad"))
        y <- x@rad.tab
    else if(is(x,"numeric")){
        z <- x[x>0]
        y <- new("rad", data.frame(rank=1:length(z), abund=sort(z, decreasing=T)))
    }
    return(y)
}
