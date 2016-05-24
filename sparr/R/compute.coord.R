compute.coord <- function(coord,data,h,n,WIN,type,q.opt=1,forAdapt=F,counts){
    XU <- (coord[1]-data[,1])/h
    YU <- (coord[2]-data[,2])/h

    result <- ((1/n)*sum((1/h^2)*counts*KERNEL(data.frame(cbind(XU,YU)),type)))/q.opt
    
    if(is.null(WIN)){
        return(result)
    } else if(inside.owin(coord[1],coord[2],WIN)){
        return(result)
    } else if(forAdapt){
        return(0)
    } else {    
        return(NA)
    }
}
