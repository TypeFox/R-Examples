print.imageQT <-function (x, ...) {


    cat("\n2D quadtree decomposition\n") 
    cat(" object of class imageQT\n----------------------------------\n")
    cat("\nsummary\n=======\n")
    cat("data:", x$data.name,"\n")
    xindl<-unlist(x$indl)
    xindS<-unlist(x$indS)
    if(length(xindl)>0){
    	cat("Indices of non-stationary sub-images:\n",as.character(xindl),"\n")
    }
    else{
    	cat("Indices of non-stationary sub-images:\n","none","\n")
    }
    if(length(xindS)>0){
    	cat("Indices of stationary sub-images:",as.character(xindS),"\n")
    }
    else{
    	cat("Indices of stationary sub-images:","none","\n")
    }
    cat("\n")
    cat("minimum testing region:",x$minsize,"\n\n") 

}
