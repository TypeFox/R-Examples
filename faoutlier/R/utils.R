#common utility functions

faoutlierClusterEnv <- new.env()
faoutlierClusterEnv$ncores <- 1L

myApply <- function(X, MARGIN, FUN, ...){
    usepar <- if(MARGIN == 1L){
            ifelse(length(X[,1]) > 2*faoutlierClusterEnv$ncores, TRUE, FALSE)
        } else {
            ifelse(length(X[1,]) > 2*faoutlierClusterEnv$ncores, TRUE, FALSE)
        }
    if(!is.null(faoutlierClusterEnv$CLUSTER) && usepar){
        return(t(parallel::parApply(cl=faoutlierClusterEnv$CLUSTER, X=X, 
                                    MARGIN=MARGIN, FUN=FUN, ...)))
    } else {
        return(t(apply(X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    }
}

myLapply <- function(X, FUN, ...){
    usepar <- ifelse(length(X) > 2*faoutlierClusterEnv$ncores, TRUE, FALSE)
    if(!is.null(faoutlierClusterEnv$CLUSTER) && usepar){
        return(parallel::parLapply(cl=faoutlierClusterEnv$CLUSTER, X=X, 
                                   fun=FUN, ...))
    } else {
        return(lapply(X=X, FUN=FUN, ...))
    }
}