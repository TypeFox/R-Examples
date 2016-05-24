"plot.traj" <-
function(x, id=levels(x$id), burst=levels(x$burst), date=NULL,
         asc=NULL, area=NULL,
         xlim=range(x$x), ylim=range(x$y),
         colasc=gray((256:1)/256), colpol="green",
         addpoints=TRUE, addlines=TRUE,
         perani=TRUE, final=TRUE,...)
{
    ## Verifications
    polygon<-area
    if (!is.null(area)) {
        if (!inherits(area, "area"))
            stop("x should be an object of class area")
    }
    if (!inherits(x, "traj"))
        stop("x should be an object of class traj")
    x <- x[!is.na(x$x),]
    if (any(is.na(x$x)))
        xlim <- range(x$x)
    if (any(is.na(x$y)))
        ylim <- range(x$y)


    ## select dates
    if (!is.null(date))
        x<-x[(x$date>=date[1])&(x$date<date[2]),]

    ## select des animals
    i<-split(x, x$id)
    x<-do.call("rbind", i[id])

    ## select bursts
    i<-split(x, x$burst)
    x<-do.call("rbind", i[burst])
    x$burst<-factor(x$burst)
    x$id<-factor(x$id)

    if (!perani)
        idc<-"burst"
    else
        idc<-"id"
    li<-split(x, x[[idc]])
    id<-levels(x[[idc]])
    opar<-par(mar=c(0.1,0.1,2,0.1),
              mfrow=n2mfrow(length(li)))
    m<-unlist(lapply(li, function(x) mean(x$date)))
    nli<-names(li)
    nli<-nli[order(m)]

    ## loop for each graph
    for (i in nli) {
        if (!is.null(asc))
            image(asc, col=colasc,
                  xlim=xlim, ylim=ylim, main=i, axes=FALSE,...)
        else
            plot(x$x,x$y, type="n", asp=1,
                 xlim=xlim, ylim=ylim, axes=FALSE,
                 main=i, ...)
        box()
        if (!is.null(polygon)) {
            pol<-split(polygon[,2:3], factor(polygon[,1]))
            for (j in 1:length(pol))
                polygon(pol[[j]], col=colpol)
        }
        if (addlines) {
            for (j in levels(factor(li[[i]]$burst))) {
                lines(x$x[x$burst==j], x$y[x$burst==j])
            }
        }
        if (addpoints) {
            for (j in levels(factor(li[[i]]$burst))) {
                points(x$x[x$burst==j],x$y[x$burst==j],pch=21,
                       col="black", bg="white")
            }
        }
        if (final) {
            for (j in levels(factor(li[[i]]$burst))) {
                points(x$x[x$burst==j][c(1,length(x$x[x$burst==j]))],
                       x$y[x$burst==j][c(1,length(x$y[x$burst==j]))],
                       pch=14, col=c("blue", "red"))
            }
        }
    }
    par(opar)
}

