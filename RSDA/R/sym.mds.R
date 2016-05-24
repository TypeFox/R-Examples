sym.mds <-
function(sym.data,distance=c("hausdorff","centers"),p=2,method=c("classic","INTERSCAL")) {
  distance<-match.arg(distance)
  method<-match.arg(method)
  if(method=="classic") 
    return(cmdscale(interval.dist(sym.data,distance)))
}
