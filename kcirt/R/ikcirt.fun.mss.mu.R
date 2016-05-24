ikcirt.fun.mss.mu <-
function(jj, iimss, rndTrys, hatMu, useSysCov, penalty ) {
    
    valgp <- rndTrys[jj]
    
    hatMu[iimss] <- valgp
    
    mxDelta <- get("mxDelta")
    mxHatLSE <- get("mxHatLSE")
    covStochastic <- get("covStochastic")

    
    hatZstar <- mxDelta %*% ( hatMu  +  mxHatLSE )
    
    
    
    if(penalty == "logit") {
        Y <- get("Y")
        xmu.shrink <- get("xmu.shrink")
        
        hatWstar <- hatZstar / sqrt(diag(covStochastic))
        hatY <- matrix( pnorm( hatWstar ), nrow(hatWstar), ncol(hatWstar) )
        our.cost <- - mean( Y * log(hatY) + (1-Y) * log(1-hatY), na.rm=TRUE )  + xmu.shrink * mean( (hatMu)^2 ) ; our.cost
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
