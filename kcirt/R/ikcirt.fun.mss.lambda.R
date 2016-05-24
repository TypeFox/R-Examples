ikcirt.fun.mss.lambda <-
function(jj, iimss, jjmss, rndTrys, mxHatLambda, penalty, usetruesigma ) {
    
    valgp <- rndTrys[jj]
    
    mxHatLambda[iimss, jjmss] <- valgp
    
    mxStHE <- get("mxStHE")
    mxDelta <- get("mxDelta")
    hatMu <- get("hatMu")
    
    mxHatLSE <- mxHatLambda %*% mxStHE

    
    
    
    hatZstar <- mxDelta %*% ( hatMu  +  mxHatLSE )
    
    
    if(penalty == "L2" | penalty == "L2c") {
        if(!usetruesigma) {
            mxSlot <- get("mxSlot")
            mxHatEta <- get("mxHatEta")
            covStochastic <- get("covStochastic")
            mxHatDLS <- mxDelta %*% mxHatLambda %*% mxSlot
            useSysCov <- mxHatDLS %*% var(mxHatEta) %*% t(mxHatDLS)   +   covStochastic
        } else {
            mxSigma <- get("mxSigma")
            useSysCov <- mxSigma
        }
    }
    
    
    if(penalty == "logit") {
        
        
        
        covStochastic <- get("covStochastic")
        Y <- get("Y")
        xlambda.shrink <- get("xlambda.shrink")
        xbool.lambdaShrinkLL <- get("xbool.lambdaShrinkLL")
        
        if(xbool.lambdaShrinkLL) {
            xshrinkTerm <- xlambda.shrink * mean( diag(mxHatLambda %*% t(mxHatLambda)) )
        } else {
            xshrinkTerm <- xlambda.shrink * mean( diag(mxHatLambda)^2 )
        }
        
        hatWstar <- hatZstar / sqrt(diag(covStochastic))
        hatY <- matrix( pnorm( hatWstar ), nrow(hatWstar), ncol(hatWstar) )
        our.cost <- - mean( Y * log(hatY) + (1-Y) * log(1-hatY), na.rm=TRUE ) + xshrinkTerm ; our.cost
    }
    
    if(penalty == "L2") {
        Z <- get("Z")
        hatWstar <- hatZstar / sqrt(diag(useSysCov))
        our.cost <- sqrt( mean( (Z - hatWstar)^2, na.rm=TRUE ) ) ; our.cost
    }
    if(penalty == "L2c") {
        varZ <- get("varZ")
        Zconv <- get("Zconv")
        hatWstar <- varZ %*% hatZstar / sqrt(diag(useSysCov))
        our.cost <- sqrt( mean( (Zconv - hatWstar)^2, na.rm=TRUE ) ) ; our.cost
    }
    if(penalty == "miscat") {
        Y <- get("Y")
        our.cost <- 1 - sum( (hatZstar > 0 & Y == 1) | (hatZstar <= 0 & Y == 0) , na.rm=TRUE ) / sum(!is.na(Y)) ; our.cost

    }

    
    return(our.cost)
    
    
}
