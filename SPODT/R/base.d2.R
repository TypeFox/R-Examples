base.d2 <- function(points.d2,d2){
    nbLig <- dim(points.d2)[1]
    nbLig.d2 <- dim(d2)[1]
    new.d2<-d2
    for(i in 1:nbLig){
        for(k in 1:nbLig.d2){
            if(is.in(d2[k,1], d2[k,2], d2[k,3], d2[k,4], points.d2[i,1], points.d2[i,2])){
                new.d2[k,] <- c(d2[k,1], d2[k,2] , points.d2[i,1], points.d2[i,2])
                new.d2 <- rbind(new.d2,c(d2[k,3], d2[k,4] , points.d2[i,1], points.d2[i,2]))
            }
        }
    }
    return(new.d2)
}