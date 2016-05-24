SPLINE= function (rw, bandwidth=32, p=0.5) {
p=as.numeric(p)
series<-rw[!is.na(rw)]
curve<-suppressWarnings(ffcsaps(series, nyrs=bandwidth, f=p))
rw[!is.na(rw)]<-curve
return(rw)
}
