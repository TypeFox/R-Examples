
AREA_BASIC <- function(ndim = 2, lower, upper, minpts = 100, functn, ...){
    
    eval.len <- ceiling(sqrt(minpts))
    eval.x <- seq(lower[1],upper[1],length=eval.len)
    eval.y <- seq(lower[2],upper[2],length=eval.len)
    eval.x.dist <- eval.x[2]-eval.x[1]
    eval.y.dist <- eval.y[2]-eval.y[1]
    
    eval.frame <- data.frame(sort(rep(eval.x,eval.len)),rep(eval.y,eval.len))
    eval.fun <- apply(eval.frame,1,functn,...)
    
    eval.fun.s <- as.vector(na.omit(eval.fun))
 	
    return(list(val=sum(eval.x.dist*eval.y.dist*eval.fun.s),eval.fun=eval.fun,eval.frame=eval.frame))
}
