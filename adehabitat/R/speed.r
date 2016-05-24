"speed" <- function(x, id=levels(x$id), burst=levels(x$burst),
                    date=NULL, units=c("seconds", "hours","days"))
{
    ## Verifications
    .Deprecated("as.ltraj")
    if (!inherits(x, "traj"))
        stop("should be an object of class traj")
    units<-match.arg(units)

    ## Selection of dates
    x<-getburst(x, burst=burst, id=id, date=date)

    ## distances between successives relocations
    li<-split(x, x$burst)
    foo<-function(x) {
        x1<-x[-1,]
        x2<-x[-nrow(x),]
        dist<-sqrt( (x1$x-x2$x)^2 + (x1$y-x2$y)^2)
        hour<-(unclass(x1$date)-unclass(x2$date))
        if (units=="hours")
            hour<-(unclass(x1$date)-unclass(x2$date))/3600
        if (units=="days")
            hour<-(unclass(x1$date)-unclass(x2$date))/(3600*24)
        disx<-(x1$x-x2$x)
        disy<-(x1$y-x2$y)
        so<-cbind.data.frame(id=x2$id,x=x2$x, y=x2$y, date=x2$date,
                             burst=x2$burst,
                             sp.x=disx/hour, sp.y=disy/hour,
                             speed=dist/hour, dt=hour)
        return(so)
    }

    ## Output
    lo<-do.call("rbind", lapply(li, foo))
    row.names(lo)<-1:nrow(lo)
    return(lo)
}

