`TimeSeriesPlot` <-
function(z, SubLength=Inf, aspect=0.25, type="l", xlab="Observation Number", ylab=NULL, main=NULL, ...){
if (SubLength==Inf) {
    pin <- par("pin")
    default.aspect <- pin[2]/pin[1]
    if(aspect < default.aspect)
        pin.new <- c(pin[1], pin[1] * aspect)
    else 
        pin.new <- c(pin[2]/aspect, pin[2])
    par(pin = pin.new)
    yl <- ylab
    ti <- main
    if (is.null(ti))
        ti <- attr(z, "title")
    plot(z, type=type,  main=ti, ylab=yl, xlab = xlab, ...)
    par(pin=pin)
    }
else { #this part requires lattice library
    n<-length(z)
    nblocks<-ceiling(n/SubLength)
    ti <- main
    if (is.null(ti))
        ti <- attr(z, "title")
    if (nblocks ==1)
        xyplot(z~(1:n) , aspect=aspect, type=type, ylab=ylab, xlab = xlab, main=ti, ...)
    else {
        y<-x<-numeric(0)
        u<-1:n
        if (n>SubLength){
            LastBlockz<-z[(n-SubLength+1):n]
            LastBlockx<-u[(n-SubLength+1):n]
            for (i in 1:(nblocks-1)) {
                ii<-seq((i-1)*SubLength+1,SubLength*i)
                y<-c(y,z[ii])
                x<-c(x,u[ii])
                }
            y<-c(y,LastBlockz)
            x<-c(x,LastBlockx)
            }
    epoch<-ordered(rep(1:nblocks,rep(SubLength,nblocks)))
    xyplot(y~x | epoch, aspect=aspect, type=type, ylab=ylab, xlab=xlab, scales=list(x="sliced",y="same"),strip=FALSE,main=ti, ...)
        }
    }
}

