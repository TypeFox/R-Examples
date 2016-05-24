distNode<-function(pts,decdeg=FALSE){
    long<-length(pts[,1])
    mydeg<-cbind(x1=pts[1:(long-1),1],y1=pts[1:(long-1),2],x2=pts[2:long,1],y2=pts[2:long,2])
    distSeg(mydeg,decdeg=decdeg)
}
