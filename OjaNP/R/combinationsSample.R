`combinationsSample` <- function(n,k,vec){
    quantity <- function(n,k,figure,digit){
        return(choose(n-figure,k-digit))
    }
    combi <- function(n,k,i){
        x <- numeric(k+1)
        for(y in 1:k){
            sumofquant <- 0
            for (z in (x[y]+1):n){
                quant <- quantity(n,k,z,y)
                sumofquant <- sumofquant + quant
                if(i <= sumofquant){x[y+1] <- z; i <- i-(sumofquant-quant); break}
            }
        }
        return(x[-1])    
    }
    return(t(sapply(X=vec,FUN=combi,n=n,k=k)))
}
