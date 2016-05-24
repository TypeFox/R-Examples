"polycirc" <-
function(radius, pts=c(0,0),nbr=50) {

# Giraudoux 02.2004 - return the polygon coordinates of
# a circle. Radius is the radius of the circle, 
# pts = vector of the coordinates of a center point
# nbr = the number of edges approximating the circle


circ=NULL
angles<-seq(0,2*pi,l=nbr)
angles<-angles[-length(angles)]
    
    for(i in angles) {
        circ<-rbind(circ,cbind(radius*sin(i),radius*cos(i)))
    }
    
    x<-cbind(circ[,1]+pts[1],circ[,2]+pts[2])
    x<-rbind(x,x[1,])
    return(x) 
}

