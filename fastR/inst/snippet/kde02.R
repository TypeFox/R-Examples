kde <- function(data,kernel=K1,...) {
    n <- length(data)
    scalingConstant=integrate(function(x){kernel(x,...)},-Inf,Inf)$value
    f <- function(x) {
        mat <- outer(x,data, FUN=function(x,data) {kernel(x-data,...)} )
        val <- apply(mat,1,sum)
        val <- val/(n*scalingConstant)
        return(val)
    }
    return(f)
}
