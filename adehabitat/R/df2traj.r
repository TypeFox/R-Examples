"df2traj" <- function(df)
{
    ## Verifications
    x<-df
    if (!inherits(x, "data.frame"))
        stop("x should be of class data.frame")

    ## Verification of the format:
    ok<-1
    if (is.null(x$id))
        ok<-0
    if (is.null(x$x))
        ok<-0
    if (is.null(x$y))
        ok<-0
    if (is.null(x$date))
        ok<-0
    if (is.null(x$burst))
        ok<-0
    if (!inherits(x$date, "POSIXct"))
        ok<-0
    if (!is.factor(x$id))
        ok<-0
    if (!is.factor(x$burst))
        ok<-0

    ## Output
    if (ok == 0)
        stop("non convenient format.\n please create the object with the function as.traj")
    class(x)<-c("traj", "data.frame")
    return(x)
}

