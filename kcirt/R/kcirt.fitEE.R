kcirt.fitEE <-
function(model, mxHatLambda, maxIter=40, lambda.ridge=0.3, Seta.ridge=0.01) {
    
    
    preSS <- solve(crossprod(model$mxSlot)) %*% t(model$mxSlot)
    
    
    preDD <- pseudoinverse(crossprod(model$mxDelta)) %*% t(model$mxDelta)
    
    
    #mx.hat.largeLambda <- mx.largeLambda + matrix(rnorm( nBlocks^2*d^2, 0, 0.1 ), nBlocks*d, nBlocks*d)
    #mx.hat.largeLambda[ mx.largeLambda == 0 ] <- 0
    Zzero <- model$Z
    Zzero[ model$Yisna ] <- 0
    
    ddZ <- preDD %*% Zzero
    
    
    dim(ddZ)
    
    ## apply(ddZZ, 1, mean)
    
    
    ###ddZZc <- ddZZ - apply(ddZZ, 1, mean)
    
    ddZc <- ( ddZ - apply(ddZ, 1, mean) ) / apply(ddZ, 1, sd)
    
    ## apply(ddZZc, 1, mean)
    
    
    ## DDYYstar <- preDD %*% t(YYstar)
    
    
    
    hatSysCov <- mxHatLambda %*% diag(1, nrow(ddZc)) %*% t(mxHatLambda) + model$covShocks
    
    priorResErr <- Inf
    mxHatEta <- NA
    
    xboolKeepGoing <- TRUE
    kk <- 0
    while(xboolKeepGoing & kk < maxIter) {
        kk <- kk + 1
        
        PRIORmxHatLambda <- mxHatLambda
        PRIORhatSysCov <- hatSysCov
        PRIORmxHatEta <- mxHatEta
        
        LLcov <- crossprod(mxHatLambda)
        invLLcov <- solve(LLcov + diag(lambda.ridge, ncol(mxHatLambda)))
        ##invLLcov <- pseudoinverse(LLcov)
        
        hatSeta <- t( invLLcov %*% t(mxHatLambda) %*% (sqrt(diag(hatSysCov)) * ddZc))
        dim(hatSeta)
        
        hatXXSeta <- crossprod(hatSeta)
        hatCovSeta <- var(hatSeta)
        
        dim(hatCovSeta)
        
        invHatXXSeta <- pseudoinverse( hatXXSeta + diag(Seta.ridge*nrow(hatSeta), ncol(hatCovSeta)) )
        #invHatXXSeta <- pseudoinverse( hatXXSeta  )
        
        mxHatLambda <- invHatXXSeta %*% t(hatSeta) %*% t(sqrt(diag(hatSysCov)) * ddZc)
        
        ################################################### constraints
        #mxHatLambda[ model$mxLambda == 0 ] <- 0
        
        mxHatLambda <- mxHatLambda * diag(1, ncol(mxHatLambda))
        
        hatddZc <- mxHatLambda %*% t(hatSeta)
        
        hatSysCov <- mxHatLambda %*% hatCovSeta %*% t(mxHatLambda) + model$covShocks ###############
        
        hatEta <- t(preSS %*% t(hatSeta))
        
        hatddZcExpanded <- hatddZc / sqrt(diag(hatSysCov))
        
        #Resids <- ddZc - hatddZcExpanded
        
        Resids <- model$mxDelta %*% ( ddZc - hatddZcExpanded )
        
        resErr <- sqrt( mean( ( Resids[ !model$Yisna ] )^2 ) )
        
        
        
        cat( "iter: ", kk, " --- ", "res err: ", resErr,  "\n")
        
        if( abs(resErr - priorResErr) < 10^(-6) ) {
            xboolKeepGoing <- FALSE
        }
        if( resErr > priorResErr ) {
            xboolKeepGoing <- FALSE
            hatSysCov <- PRIORhatSysCov
            mxHatLambda <- PRIORmxHatLambda
            mxHatEta <- PRIORmxHatEta
            resErr <- priorResErr
        }
        
        priorResErr <- resErr
        
        
    }
    
    hatMu <- apply(ddZ, 1, mean) * sqrt(diag(hatSysCov))
    
    if(!is.null(model$mxEta)) {
        performance <- rep(NA, model$nuc)
        for(i in 1:model$nuc) {
            performance[i] <- cor(model$mxEta[ ,i], hatEta[ ,i])^2
        }
        model[["performance"]] <- performance
    }
    
    model[["hatMu"]] <- hatMu
    ###### model[["hatSysCov"]] <- hatSysCov #### this is NOT the system variance !!!!!!!!!!!!!!!!!!!!!
    model[["mxHatLambda"]] <- mxHatLambda
    model[["mxHatEta"]] <- hatEta
    model[["L2resErr"]] <- resErr
    
    
    return(model)
    
    
}
