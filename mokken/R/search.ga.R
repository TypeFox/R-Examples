"search.ga" <- function(X, popsize, maxgens, alpha, critval, pxover, pmutation){
             
     nitem <- ncol(X)
     item.label <- dimnames(X)[[2]]
        
     if(nitem <= 10)      iter <- maxgens
     else if(nitem <= 20) iter <- 5000
     else if(nitem > 20)  iter <- round(4000/(nitem/20),0)
     
     if(maxgens < iter) iter <- maxgens
     
     npers <- nrow(X)
     variance <- c(var(X))
     max.variance <- c(var(apply(X, 2, sort)))
     SijMatrix <- c( outer(apply(X, 2, var), apply(X, 2, var), "*"))
     Population <- rep(0,nitem*(popsize+2))
     itercount <- 0
     fitness <- rep(0,(popsize+2)*3)
    
     while(1){
        Output <- .C("GeneticAlgorithm",as.integer(Population),as.integer(itercount),as.integer(popsize),
                  as.integer(nitem),as.integer(npers),as.integer(iter),as.double(pxover),as.double(pmutation),
                  as.double(critval),as.double(alpha),as.double(variance),as.double(max.variance),
                  as.double(SijMatrix), as.double(fitness))
        itercount <- Output[[2]]
        Population <- Output[[1]]
        fitness <- c(Output[[14]][1:(popsize+2)],rep(0,(popsize+2)*2))
        if(itercount == ceiling(maxgens/iter)) break   
     }
     InSet <- as.matrix(c(matrix(Output[[1]],popsize+2,nitem,byrow=T)[popsize+1,]))
     dimnames(InSet) <- list(item.label,"Scale")
     return(InSet)
}
