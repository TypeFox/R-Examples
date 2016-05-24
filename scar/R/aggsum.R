aggsum <- function(x,reverse=FALSE){
    n = length(x)
    if(n>1){
        if(reverse==FALSE){y=.C("AggSum",as.integer(n), x = as.double(x))$x} else
        {y =rev(x); y=.C("AggSum",as.integer(n), x = as.double(y))$x; y = rev(y)}
    }
    return(y)
}
