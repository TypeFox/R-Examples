cutPoints <-
function(x,y){
    od <- order(x)    
    xo <- x[od]
    yo <- y[od]    
    depth <- 1     
    
    gr <- function(low,upp,depth=depth){ 
       x <- xo[low:upp]  
       y <- yo[low:upp]  
       n <- length(y) 
       ct <- cutIndex(x,y)
       if(is.null(ct)) return (NULL) ## when cut index=NULL
       ci <- ct[1]; entropy <- ct[2]
       ret <- mdlStop(ci,y,entropy) # MDL Stop
       if(is.null(ret)) return(NULL)
       return(c(ci,depth+1)) 
} 
   
## xo: original x in ascending order of x; 
## yo: original y reordered in ascending order of x 
part <- function(low=1, upp=length(xo), cutPoints=NULL,depth=depth){
       x <- xo[low:upp]
       y <- yo[low:upp]
       n <- length(x)
       if(n<2) return (cutPoints)
       cc <- gr(low, upp, depth=depth)
       ci <- cc[1]
       depth <- cc[2]
       if(is.null(ci)) return(cutPoints)
       cutPoints <- c(cutPoints,low+ci-1)
       cutPoints <- as.integer(sort(cutPoints))
       return(c(part(low, low+ci-1,cutPoints,depth=depth), 
                part(low+ci,upp,cutPoints,depth=depth)))
    }
    
    res <- part(depth=depth)
    ci <- NULL ;cv <- numeric()
    if(!is.null(res)) {
        ci <- as.integer(res)
        cv <- (xo[ci]+xo[ci+1])/2
    }
    res <- unique(cv)## returns cutIndex and cutValues
    return(res)
}
