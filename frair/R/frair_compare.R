# frair_compare()
# Impliments the difference modelling using dummary variables promoted by Julliano in Juniano 2001, pg 193
# Seems sensible. Much better than t-tests on bootstrapped coefficients!

frair_compare <- function(frfit1, frfit2, start=NULL){
    if(!inherits(frfit1, 'frfit') | !inherits(frfit2, 'frfit')){
        stop('Both inputs must be of class frfit')
    }
    if(frfit1$response!=frfit2$response){
        stop('Both inputs must be fitted using the same response.')
    }
    # Get the name of the 'XX_nll_diff' function
    #fr_nll_difffunc <- get(paste0(unlist(frair_responses(show=FALSE)[[frfit1$response]])[1],'_nll_diff'), pos = "package:frair")
    fr_nll_difffunc <- get(paste0(frfit1$response,'_nll_diff'), pos = "package:frair")
    
    if(any(frfit1$optimvars!=frfit2$optimvars)){
        stop('Both inputs must have the same optimised variables.')
    }
    # Get resonable starting values for fitted coefs (the mean of the two)
    if(is.null(start)){
        start <- list()
        for(a in 1:length(frfit1$optimvars)){
            cname <- frfit1$optimvars[a]
            start[cname] <- mean(frfit1$coefficients[cname], frfit2$coefficients[cname])
        }
    } else {
        fr_checkstart(start, deparse(substitute(start)))
    }
    
    # Get the manes of the 'delta' variables (e.g. Da, Dh) and setup starting values for them
    varnames <- unlist(frair_responses(show=FALSE)[[frfit1$response]][4])
    deltavarnames <- NULL
    for(a in 1:length(varnames)){
        deltavarname <- paste0('D',varnames[a])
        start[deltavarname] <- 0
        deltavarnames <- c(deltavarnames,deltavarname)
    }
    
    if(any(frfit1$fixedvars!=frfit2$fixedvars)){
        stop('Both inputs must have the same fixed variables.')
    }
    # Get fixed values
    fixed=list()
    for(a in 1:length(frfit1$fixedvars)){
        fname <- frfit1$fixedvars[a]
        if(frfit1$coefficients[fname]!=frfit2$coefficients[fname]){
            stop('Fixed variables must have the same numerical value')
        }
        fixed[fname] <- frfit1$coefficients[fname]
    }
    
    # Get X and Y and setup dummy coding for the model
    Xin <- c(frfit1$x,frfit2$x)
    Yin <- c(frfit1$y,frfit2$y)
    grp <- c(rep(0,times=length(frfit1$x)), rep(1,times=length(frfit2$x)))
    
    # https://github.com/dpritchard/frair/issues/23
    if(length(unlist(start))>1){
        try_test <- try(bbmle::mle2(minuslogl=fr_nll_difffunc, start=start, fixed=fixed, 
                             data=list('X'=Xin, 'Y'=Yin, grp=grp), optimizer='optim', 
                             method='Nelder-Mead', control=list(maxit=5000)), 
                        silent=TRUE)
    } else {
        try_test <- try(bbmle::mle2(minuslogl=fr_nll_difffunc, start=start, fixed=fixed, 
                             data=list('X'=Xin, 'Y'=Yin, grp=grp), optimizer='optim', 
                             control=list(maxit=5000)), 
                        silent=TRUE)
    }
    ## End https://github.com/dpritchard/frair/issues/23
    
    if(inherits(try_test, 'try-error')){
        stop(paste0('Refitting the model for the test failed with the error: \n', try_test[1], '\nNo fallback exists, please contact the package author.'))
    }
    
    # Get output from bbmle::mle2 and calculate statistics
    cmatall <- cbind(Estimate = try_test@coef, 'Std. Error' = sqrt(diag(try_test@vcov)))
    zval <- cmatall[,'Estimate']/cmatall[,'Std. Error']
    pval <- 2*pnorm(-abs(zval))
    cmatall <- cbind(cmatall,'z value'=zval,'Pr(z)'=pval)
    
    # Extract the Delta parameters
    coefmatDeltas <- cmatall[deltavarnames,]
    if(is.null(dim(coefmatDeltas))){
        dim(coefmatDeltas) <- c(1,length(coefmatDeltas))
        dimnames(coefmatDeltas) <- list(deltavarnames, c('Estimate', 'Std. Error', 'z value', 'Pr(z)'))
    }
    
    # Get orginal coefs
    origcoef <- rbind(coef(frfit1)[frfit1$optimvars],coef(frfit2)[frfit1$optimvars])
    row.names(origcoef) <- c(deparse(substitute(frfit1)), deparse(substitute(frfit2)))
    
    # Print the output
    cat('FUNCTIONAL RESPONSE COEFFICIENT TEST\n\n')
    
    cat('Response:            ', frfit1$response, '\n', sep='')
    cat('Optimised variables: ', paste(frfit1$optimvars, collapse=','), '\n', sep='')
    cat('Fixed variables:     ', paste(frfit1$fixedvars, collapse=','), '\n\n', sep='')
    
    cat('Original coefficients: \n')
    print(round(origcoef,5))
    
    cat('\n')
    cat('Test: ', deparse(substitute(frfit1)), ' - ', deparse(substitute(frfit2)), '\n\n', sep='')
    
    printCoefmat(round(coefmatDeltas,5))
    
    output <- list(frfit1=frfit1, frfit2=frfit2, test_fit=try_test, result=coefmatDeltas)
    class(output) <- c('frcompare', class(output))
    invisible(output)
}
