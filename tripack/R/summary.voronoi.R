summary.voronoi<-function(object,...)
{
  if(!inherits(object,"voronoi"))
    stop("object must be of class \"voronoi\"")
  ans<-list(nn=length(object$x),
            nd=length(object$dummy.x),
            call=object$call)
  class(ans)<-"summary.voronoi"
  ans
}
