align2D <-
function(lat, long, x1, x2, wts=NULL){
    ## Get best fit in space of (latitude, longitude)
    if(is.null(wts))wts <- rep(1,length(x1))
    fitlat <- predict(lm(lat ~ x1+x2, weights=wts))
    fitlong <- predict(lm(long ~ x1+x2, weights=wts))
    list(fitlat = fitlat, fitlong=fitlong, lat=lat, long=long)
}
