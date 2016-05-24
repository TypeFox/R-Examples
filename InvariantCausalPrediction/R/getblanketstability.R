getblanketstability <-
function(X,Y,maxNoVariables=10,maxNoVariablesSimult=5){
    p <- ncol(X)
    
    if(p <=maxNoVariables){
        usevar <- 1:p
    }else{
        nsim <- 20
        usevarall <- numeric(0)
        for (sim in 1:nsim){
            ind <- sample(1:nrow(X),round(0.5*nrow(X)))
            lobs <- coef(glmnet(X[ind,],Y[ind],family= if(!is.factor(Y)) "gaussian" else "binomial"))[-1,]
            nnz <- apply(lobs!=0,2,sum)
            nnzsel <- 0
            usevar <- numeric(0)
            while(length(usevar)<maxNoVariables & nnzsel<max(nnz)){
                nnzsel <- nnzsel+1
                sel <- which( nnz==nnzsel)
                if(length(sel)>0) usevar <- sort(unique(c(usevar, which(apply(lobs[,sel,drop=FALSE]!=0,1,any)))))
            }
            usevarall <- c(usevarall,usevar)
        }
        tab <- table(usevarall)
        tab <- tab[order(tab,decreasing=TRUE)]
        usevar <- as.numeric(names(tab))[1:min(maxNoVariables,length(tab))]
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
