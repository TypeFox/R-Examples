topdown <-
function(data, method=1){
    maxP <- 28 ## maximum number of partitions
    p1 <- length(data[1,])
    p <- p1-1
    y <- as.integer(data[,p1])
    S <- dim(table(y))[1]
    n <- length(y)
    dCuts <- array(list(),p) ## cut points obtained for p X-variables
    
    for (i in 1:p){
        x <- data[,i]; od <- order(x)
        xo <- x[od]; yo <- y[od] 
        
        GlobalC <- 0; cacc <- 0; k <- 1; addCut <- NULL ;ci <- NULL        
        ci <- which(diff(xo)!=0)
        if(!is.null(ci)) cuts <- (xo[ci]+xo[ci+1])/2
        bd <- c(xo[1], xo[n])
        di <- cuts    
        
        while(length(di)>0) {
            ret <- findBest(xo,yo,bd,di,method)
            cacc <- ret$cacc 
            Dp <- insert(ret$addCut, bd)  ## sort(c(bd,ret$addCut))
            di <- ret$newDi 
                      
            if((cacc > GlobalC || k<S) & (k<n)){
                bd <- Dp
                GlobalC <- cacc 
                k <- k+1
            } else {
                Dp <- bd
                break
            }
        }      
        dCuts[[i]] <- bd        
    }
    return(dCuts)
}
