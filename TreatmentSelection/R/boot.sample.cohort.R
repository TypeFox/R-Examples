boot.sample.cohort <-
function( event, trt, rho){

    #boot1 <- sample(1:sum(trt)  ,replace=TRUE)

    #sample conditional on trt
    #trt.b <- c(trt[trt==0][boot0], trt[trt==1][boot1])
    #Y.b <- c(Y[trt==0][boot0], Y[trt==1][boot1])
    #event.b <- c(event[trt==0][boot0], event[trt==1][boot1])

    # trt.b <- c(trt[boot0])
    # Y.b <- c(Y[boot0])
    # event.b <- c(event[boot0])
 
    # Y2.b <- c(Y2[boot0])
    Resample <- TRUE

    tryN = 0

    while(Resample){
    
    ind <- sample(length(trt),replace=TRUE)
        
    Resample <- any(table(event[ind], trt[ind]) == 0) | length(unique(trt[ind]))==1 | length(unique(event[ind]))==1
    

    if(Resample){ warning("Bootstrap sample failed to sample individuals in each trt/event strata. Had to resample") ; tryN = tryN +1}
    if(tryN >100) stop("Had to resample too many bootstrap replicates due to too few observations in trt/event strata")
    }
   
return(c(rep(-9999, 7), ind))

}
