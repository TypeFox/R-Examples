ikcirt.df1 <-
function(model, lambdaConstraint="self") {
    
    objvec <- c(ncol(model$Y), ncol(model$mxData), ncol(model$Z), nrow(model$mxEta), nrow(model$mxHatEta))
    if(length(objvec) < 1) { return("Cannot calculate df.  Model requires one of Y, Z, mxData, mxEta, or mxHatEta.") }
    totrq <- objvec[1]

    xdf <- model$nBlocks * totrq - ( length(model$mu) + ncol(model$Y)*model$nuc )
    
    lcti <- model$mxLambdaCTinfo
    
    iimss <- rep( I(1:nrow(lcti)), ncol(lcti) )
    jjmss <- rep( I(1:ncol(lcti)), each=nrow(lcti) )
    
    lcti <- as.vector(lcti)
    
    if(lambdaConstraint=="self") {
        lnp <- length(diag(model$mxLambda))
    }

    if(lambdaConstraint == "withinx") {
        lnp <- sum(lcti == "S" | lcti == "WF")
    }
    
    if(lambdaConstraint == "withini") {
        lnp <- sum(lcti == "S" | lcti == "WF" | lcti == "WT")
    }

    if(lambdaConstraint == "betweenx") {
        lnp <- sum(lcti == "S" | lcti == "WF" | lcti == "BF")
    }
    if(lambdaConstraint == "betweeni") {
        lnp <- sum(lcti == "S" | lcti == "WF" | lcti == "BF" | lcti == "WT" | lcti == "BT") 
    }
    
    if(lambdaConstraint == "priorx") {
        lnp <- sum(lcti == "S" | lcti == "WF" | (jjmss < iimss & lcti == "BF") )
    }
    if(lambdaConstraint == "priori") {
        lnp <- sum(lcti == "S" | lcti == "WF" | lcti == "WT" | (jjmss < iimss & ( lcti == "BF" | lcti == "BT" ) ) )
    }
    appxDf <- xdf - lnp
    
    names(appxDf) <- "appx-df"
    return(appxDf)

}
