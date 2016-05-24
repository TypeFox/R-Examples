insert <-
function(x,a){
    p <- length(a)
    i <- which(a>x)
    len <- length(i)
    if(len==p) return(c(x,a))
    if(len==0) return (c(a,x))
    i1 <- i[1]
    return(c(a[1:(i1-1)],x,a[i1:p]))
}
