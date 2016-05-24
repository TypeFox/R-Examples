## time=(tic_interval, show_delay, from, to)
## fades: 0--no fading, n--n steps of fading
animation<-function(xyt, poly=NULL, win=NULL, time=NULL,
                    fades=0)
{ 
    if(is.null(poly)) {
        xyrng<-apply(xyt[,1:2], 2, range, na.rm=T)
        poly<-cbind(xyrng[c(1,2,2,1),1], xyrng[c(1,1,2,2),2])
    }
    if(is.null(win)){
        polyrng<-apply(poly, 2, range, na.rm=T)
        win<-c(polyrng[1,1:2], polyrng[2,1]-polyrng[1,1], 
               polyrng[2,2]-polyrng[1,2])
    }
    if(is.null(time)){
        trng<-range(xyt[,3], na.rm=T)
        time<-c((trng[2]-trng[1])/10, 1, trng[1], trng[2])
    }
    .C("animation", as.double(xyt), as.integer(dim(xyt)[1]),
       as.double(poly), as.integer(dim(poly)[1]),
       as.double(win), as.double(time), as.integer(fades),
       PACKAGE="spatialkernel")
}
