"getburst" <- function(x, burst=levels(x$burst),
                   id=levels(x$id), date=NULL)
{
    ## Verifications
    if (!inherits(x, "traj"))
        stop("should be an object of class traj")

    ## selection of dates
    if (!is.null(date))
        x<-x[(x$date>=date[1])&(x$date<date[2]),]

    ## selection of animals
    i<-split(x, x$id)
    x<-do.call("rbind", i[id])

    ## selection of the circuits
    i<-split(x, x$burst)
    x<-do.call("rbind", i[burst])

    ## Output
    x$burst<-factor(x$burst)
    x$id<-factor(x$id)
    return(x)
}

