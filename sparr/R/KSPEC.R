KSPEC <- function(location,data,h,type,counts){
    
    if(any(is.na(location))) return(NA)
    
    d <- ncol(data)
    n <- nrow(data)
   
    datavec <- data[,1]
    if(d > 1) {
        for(i in 2:d) datavec <- append(datavec,data[,i])
    }
    
    umat <- t((location-matrix(datavec,d,n,byrow=T))/h)
    result <- (1/(n*h^d))*sum(counts*apply(umat,1,KERNEL,type=type,dim=d))
	
    return(result)
}
