`createHyperplaneSample` <- function(N,pN,n,k){
    FailureMessage <- "Can't allocate enough memory: matrix 'X' is too large."
    op <- options(warn=-1)  
    on.exit(options(op))
    if (pN < N){  
        Sample <- try(sample(N,size=pN,replace=FALSE),silent=TRUE) # draws pN distinct numbers out of 1,...,N
        if (inherits(Sample,"try-error")){
            Sample <- try(floor((N*runif(n=pN))+0.5),silent=TRUE)
        } 
        if (inherits(Sample,"try-error")){
            stop(FailureMessage)         
        }
        else{
            CombinationsSample <- try(t(combn(x=n,m=k)[,Sample,drop=FALSE]),silent=TRUE)
            if (inherits(CombinationsSample,"try-error")){
                CombinationsSample <- try(combinationsSample(n=n,k=k,vec=Sample),silent=TRUE)
            }
        }
    }else{
        CombinationsSample <- try(t(combn(x=n,m=k)),silent=TRUE)          
    }
    if (inherits(CombinationsSample,"try-error")){
        stop(FailureMessage)
    }
    return(CombinationsSample)
}
