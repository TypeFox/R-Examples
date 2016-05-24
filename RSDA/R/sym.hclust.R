sym.hclust <-
function(sym.data,distance=c("hausdorff","centers"),p=2,method=c("ward","single","complete",
                    "average","mcquitty","median","centroid"), members=NULL) {
  distance<-match.arg(distance)
  method<-match.arg(method)
  return(hclust(interval.dist(sym.data,distance),method))
}
