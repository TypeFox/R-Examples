"wawotest.ltraj" <- function(x,...)
{
    if (!inherits(x, "ltraj"))
        stop("x should be of class \"ltraj\"")
    if ((!is.regular(x))&attr(x,"typeII"))
        stop("x should be regular or of type I")
    foo <- function(df){
        res=apply(df[,c("dx","dy","dist")],2,wawotest)
        return(res)
    }
    return(lapply(x,foo))
}

