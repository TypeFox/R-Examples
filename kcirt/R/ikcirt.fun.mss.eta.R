ikcirt.fun.mss.eta <-
function(jj, iimss, jjmss, rndTrys, mxHatEta, penalty, usetruesigma ) {
    
    valgp <- rndTrys[jj]
    mxHatEta[iimss, jjmss] <- valgp
    
    
    mxDelta <- get("mxDelta")
    mxHatLS <- get("mxHatLS")
    hatMu <- get("hatMu")
    xeta.shrink <- get("xeta.shrink")
    
    
    
    
    
    mxHatLSE <- mxHatLS %*% mxHatEta[ iimss, ]
    
    hatzzstar <- mxDelta %*% ( hatMu  +  mxHatLSE )
    
    if(penalty == "L2" | penalty == "L2c") {
        if(!usetruesigma) {
            mxHatLambda <- get("mxHatLambda")
            mxSlot <- get("mxSlot")
            covStochastic <- get("covStochastic")
            mxHatDLS <- mxDelta %*% mxHatLambda %*% mxSlot
            useSysCov <- mxHatDLS %*% var(mxHatEta) %*% t(mxHatDLS)   +   covStochastic
        } else {
            useSysCov <- get("mxSigma")
            #useSysCov <- mxSigma
        }
    }
    
    if(penalty == "logit") {
        Y <- get("Y")
        covStochastic <- get("covStochastic")
        yy=Y[ , iimss]
        hatwwstar <- hatzzstar / sqrt(diag(covStochastic))
        hatyy <- pnorm( hatwwstar )
        our.cost <- - mean( yy * log(hatyy) + (1-yy) * log(1-hatyy), na.rm=TRUE ) + xeta.shrink * mean( mxHatEta[iimss, ]^2 ) ; our.cost
    }
    
    if(penalty == "L2") {
        Z <- get("Z")
        hatwwstar <- hatzzstar / sqrt(diag(useSysCov))
        our.cost <- sqrt( mean( (Z[ , iimss] - hatwwstar)^2, na.rm=TRUE ) ) ; our.cost
    }
    
    if(penalty == "L2c") {
        #varZ <- get("varZ")
        #Zconv <- get("Zconv")
        #hatwwstar <- varZ %*% hatzzstar / sqrt(diag(useSysCov))
        #our.cost <- sqrt( mean( (Zconv[ , iimss] - hatwwstar)^2, na.rm=TRUE ) ) ; our.cost
        
        
        hatwwstar <- get("varZ") %*% hatzzstar / sqrt(diag(useSysCov))
        our.cost <- sqrt( mean( ( get("Zconv")[ , iimss] - hatwwstar)^2, na.rm=TRUE ) ) ; our.cost
        
    }
    
    if(penalty == "miscat") {
        Y <- get("Y")
        yy=Y[ , iimss]
        our.cost <- 1 - sum( (hatzzstar > 0 & yy == 1) | (hatzzstar <= 0 & yy == 0) , na.rm=TRUE ) / sum(!is.na(yy)) ; our.cost
        
    }
    
    return(our.cost)
    
    
}
