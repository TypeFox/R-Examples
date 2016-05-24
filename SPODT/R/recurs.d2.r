recurs.d2 <- function(points.d2,d2){
    new.d2<-base.d2(points.d2,d2)
    if (dim(new.d2)[1]>dim(d2)[1]){
        new.d2<-recurs.d2(points.d2,new.d2)
    }
return(new.d2)
}