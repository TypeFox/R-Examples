getblanketlasso <-
function(X,Y,maxNoVariables=10,maxNoVariablesSimult=5){
    p <- ncol(X)
    
    if(p <=maxNoVariables){
        usevar <- 1:p
    }else{
        lobs <- coef(glmnet(X,Y,family= if(!is.factor(Y)) "gaussian" else "binomial"))[-1,]
        nnz <- apply(lobs!=0,2,sum)
        nnzsel <- 0
        usevarcandidate <- numeric(0)
        usevar <- numeric(0)
        while(length(usevarcandidate)<maxNoVariables & nnzsel<max(nnz)){
            nnzsel <- nnzsel+1
            sel <- which( nnz==nnzsel)
            if(length(sel)>0) usevarcandidate <- sort(unique(c(usevarcandidate, which(apply(lobs[,sel,drop=FALSE]!=0,1,any)))))
            if(length(usevarcandidate)<=maxNoVariables) usevar <- usevarcandidate
        }
        # the following applies if several variables enter at once
    }
    testsets <- list()
    if(length(usevar)>0){
        for (ic in ((1:2^length(usevar))-1)){
            testsets[[ic+1]] <- usevar[which( ((ic %/% 2^(0:(length(usevar)-1))) %% 2 )==1)]
        }
    }
    testsets <- unique(testsets)
    le <- sapply(testsets,length)
    testsets <- testsets[ keep <- which(le>0 & le <= maxNoVariablesSimult) ]
    testsets <- testsets[order(le[keep])]
    return(testsets)
}
