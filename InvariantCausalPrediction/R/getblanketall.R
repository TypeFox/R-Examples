getblanketall <-
function(X,Y,maxNoVariables=10,maxNoVariablesSimult=5){
    p <- ncol(X)
    usevar <- 1:p
    
    testsets <- list()
    if(length(usevar)>0){
        for (ic in ((1:2^length(usevar))-1)){
            testsets[[ic+1]] <- usevar[which( ((ic %/% 2^(0:(length(usevar)-1))) %% 2 )==1)]
        }
    }
    testsets <- unique(testsets)
    le <- sapply(testsets,length)
    testsets <- testsets[order(le)]
    return(testsets)
}
