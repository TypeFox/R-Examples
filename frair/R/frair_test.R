# frair_test: Impliments the phenomenological test based on logistic regressions
frair_test <- function(formula, data){
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf_list <- as.list(mf)
    expandmod <- terms(formula(mf_list$formula), data=data)
    expandform <- formula(expandmod)
    leftside <- all.vars(expandform[[2]])
    rightside <- all.vars(expandform[[3]])
    if(length(leftside)!=1 || length(rightside)!=1) {
        stop('Currently only formulae with one dependent and one independent variable (e.g. y ~ x) are supported.')
    }
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    # moddata is the a nice clean frame
    moddata <- eval.parent(mf)
    names(moddata) <- c('eaten', 'density')
    moddata$noteaten <- moddata$density-moddata$eaten
    
    ## Go time!
    # Setup output
    out <- list('call' = call, 'x' = moddata$density, 'y'=moddata$eaten/moddata$density, 'xvar' = rightside, 'yvar' = 'Proportion')
    class(out) <- c('frtest', class(out))
    
    ## DO IT
    modT2 <- glm(cbind(eaten,noteaten)~density, data=moddata, family=binomial)
    modT3 <- glm(cbind(eaten,noteaten)~density+I(density^2), data=moddata, family=binomial)
    ## END
    out[['modT2']] <- modT2
    out[['modT3']] <- modT3
    
    # Finally return our object!
    return(out)
}

# A simple print method for the above
print.frtest <- function(x, ...){
    # Get T2 coefs
    T2Coef <- coef(summary(x$modT2))
    T2CoefOut <- T2Coef[2,]
    dim(T2CoefOut) <- c(1,ncol(T2Coef))
    dimnames(T2CoefOut) <- list(rownames(T2Coef)[2], colnames(T2Coef))
    # Get T3 coefs
    T3Coef <- coef(summary(x$modT3))
    T3CoefOut <- T3Coef[2:3,]
    dim(T3CoefOut) <- c(2,ncol(T3Coef))
    dimnames(T3CoefOut) <- list(rownames(T3Coef)[2:3], colnames(T3Coef))
    
    #t1sig <- t2sig <- t3sig <- FALSE
    t2sig <- t3sig <- FALSE
    
    cat('FUNCTIONAL RESPONSE TEST\n\n')
    
    # If the first coefficient isn't significant, it is flat and therefore a type-I
    #if(T2CoefOut[1,4]>=0.05){
    #    t1sig <- TRUE
    #    cat('Evidence for type-I reponse:\tYes\n')
    #    cat('Evidence for type-II reponse:\t-\n')
    #    cat('Evidence for type-III reponse:\t-\n\n')
    #
    # If the first coefficient is negative and significant, it is a type-II
    #} else if(T2CoefOut[1,1]<0 & T2CoefOut[1,4]<0.05){
    
    # If the first coefficient is negative and significant, it is a type-II
    if(T2CoefOut[1,1]<0 & T2CoefOut[1,4]<0.05){
        t2sig <- TRUE
        #cat('Evidence for type-I reponse:\tNo\n')
        cat('Evidence for type-II reponse:\tYes\n')
        cat('Evidence for type-III reponse:\t-\n\n')
    
    # If nothing has been triggered yet then it might be a type-III
    # The first must be positve, and significant
    # The second must be negative, and significant
    } else if(T3CoefOut[1,1]>0 & T3CoefOut[1,4]<0.05 & 
              T3CoefOut[2,1]<0 & T3CoefOut[2,4]<0.05){
        t3sig <- TRUE
        #cat('Evidence for type-I reponse:\tNo\n')
        cat('Evidence for type-II reponse:\tNo\n')
        cat('Evidence for type-III reponse:\tYes\n\n')
    
    # If nothing else has been triggered, we might have a problem
    } else {
        cat('No evidence for any response!\n')
    }
    
    #if(t1sig){
    #    cat('Type-I logistic regression output:\n')
    #    printCoefmat(T2CoefOut)
    #}
    
    if(t2sig){
        cat('Type-II logistic regression output:\n')
        printCoefmat(T2CoefOut)
    }
    
    if(t3sig){
        cat('Type-III logistic regression output:\n')
        printCoefmat(T3CoefOut)
    }
}