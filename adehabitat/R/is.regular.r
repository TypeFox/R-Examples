is.regular <- function(ltraj)
{
    if (!inherits(ltraj,"ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    if (!attr(ltraj, "typeII")) {
        return(FALSE)
    } else {
        return(sum(unlist(lapply(ltraj, function(x) {
            ddt <- x$dt[-nrow(x)]
            return(sum(abs(ddt-x$dt[1]) > 1e-7))
        }))) ==0 )
    }
}
