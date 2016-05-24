kcirt.fitMSS <-
function(model, lambdaConstraint="self", kcpus=2, penalty="logit", usetruesigma=TRUE, mss.sd=0.2, nsearch=19, l2zvarpow=0, xmu.shrink=0, xlambda.shrink=0, xeta.shrink=0.4) {
    

    
    if(length(mss.sd)==1) { mss.sd <- rep(mss.sd, 3) }
    
    hatMu <- model$hatMu
    mxSlot <- model$mxSlot
    mxDelta <- model$mxDelta
    mxHatLambda <- model$mxHatLambda
    mxHatEta <- model$mxHatEta
    covShocks <- model$covShocks
    covStochastic <- model$covStochastic
    
    mxLambdaCTinfo <- model$mxLambdaCTinfo
    
    mxSigma <- model$mxSigma
    
    Y <- model$Y
    Z <- model$Z
    
    
    if(penalty == "L2c") {
        Z[is.na(Z)] <- 0
        if(usetruesigma) {
            varZ <- mxSigma
        } else {
            varZ <- mpower( (ncol(Z)^(-1) * tcrossprod(Z)), l2zvarpow )
        }
        Zconv <- varZ %*% Z
    }
    
    
    
    if(usetruesigma) {
        useSysCov <- mxSigma
    } else {
        mxDhLS <- mxDelta %*% mxHatLambda %*% mxSlot
        hatSysCov <- mxDhLS %*% var(mxHatEta) %*% t(mxDhLS)   +   covStochastic
        useSysCov <- hatSysCov
    }
    
    
    
    mxHatLSE <- mxHatLambda %*% mxSlot %*% t(mxHatEta)
    
    hatZstar <- mxDelta %*% ( hatMu  +  mxHatLSE )
    
    
    if(penalty == "logit") {
        hatWstar <- hatZstar / sqrt(diag(covStochastic))
        hatY <- matrix( pnorm( hatWstar ), nrow(hatWstar), ncol(hatWstar) )
        our.cost <- - mean( model$Y * log(hatY) + (1-model$Y) * log(1-hatY), na.rm=TRUE ) ; our.cost
    }
    if(penalty == "L2") {
        hatWstar <- hatZstar / sqrt(diag(useSysCov))
        our.cost <- sqrt( mean( (Z - hatWstar)^2, na.rm=TRUE ) ) ; our.cost
    }
    if(penalty == "L2c") {
        hatWstar <- varZ %*% hatZstar / sqrt(diag(useSysCov))
        our.cost <- sqrt( mean( (Zconv - hatWstar)^2, na.rm=TRUE ) ) ; our.cost
    }
    if(penalty == "miscat") {
        our.cost <- 1 - sum( (hatZstar > 0 & Y == 1) | (hatZstar <= 0 & Y == 0) , na.rm=TRUE ) / sum(!is.na(Y)) ; our.cost
        
    }
    cat("starting", penalty, "cost:", our.cost, "\n")
    
    
    iimss <- 1
    
    
    
    
    
    
    run.parallel <- TRUE
    sfInit(parallel=TRUE, cpus=kcpus)
    
    sfExport("Y", "Z", "mxDelta", "mxSlot", "covShocks", "covStochastic", "mxSigma")
    
    if(penalty == "L2c") {
        sfExport("varZ", "Zconv")
    }
    if(penalty == "logit") {
        sfExport("xmu.shrink", "xlambda.shrink", "xeta.shrink")
    }
    
    #sfExport("varZ", "Zconv", "halfTotRespDivNuc")
    
    
    sfExport("mxHatLSE")
    
    if(mss.sd[1] != 0) {
        
        cat("Locating mu ...", "\n")
        for(iimss in 1:length(hatMu)) {
            
            rndTrys <- hatMu[iimss] + c( 0, rnorm(nsearch, 0, mss.sd[1]) )
            
            if( run.parallel ) {
                sfOut <- sfClusterApplyLB( x=1:length(rndTrys), fun=ikcirt.fun.mss.mu, iimss=iimss, rndTrys=rndTrys, hatMu=hatMu,
                useSysCov=useSysCov,
                penalty=penalty )
                
            }
            cost.vec <- unlist(sfOut)
            cost.vec[is.na(cost.vec)] <- Inf
            cost.vec[is.nan(cost.vec)] <- Inf
            
            xndx.mincost <- which.min(cost.vec)
            best.cost <- cost.vec[ xndx.mincost ]
            
            hatMu[iimss] <- rndTrys[xndx.mincost]
            
            #cat(best.cost, "\n")
            
            txtProgressBar(min = 0, max = length(hatMu), initial = iimss, char = "=", width = 100, style=3)
        }
        cat("\n\n")
    }
    
    ######################################## MSS hatLambda search
    
    iimss <- 1 ; jjmss <- 1
    
    
    
    mxStHE <- mxSlot %*% t(mxHatEta)
    
    sfExport("mxStHE", "mxHatEta", "hatMu")
    
    if(mss.sd[2] != 0) {
        
        if(lambdaConstraint == "self") {
            xbool.lambdaShrinkLL <- FALSE
        } else {
            xbool.lambdaShrinkLL <- TRUE
        }
        sfExport("xbool.lambdaShrinkLL")


        cat("Locating Lambda ...", "\n")
        for(iimss in 1:nrow(mxHatLambda)) {
            for(jjmss in 1:ncol(mxHatLambda)) {
                
                this.Linfo <- mxLambdaCTinfo[iimss, jjmss]
                
                if( (lambdaConstraint == "self" & this.Linfo == "S") |
                
                ( lambdaConstraint == "withinx" & (this.Linfo == "S" | this.Linfo == "WF") ) |
                ( lambdaConstraint == "withini" & (this.Linfo == "S" | this.Linfo == "WF" | this.Linfo == "WT") ) |
                ( lambdaConstraint == "betweenx" & (this.Linfo == "S" | this.Linfo == "WF" | this.Linfo == "BF") ) |
                ( lambdaConstraint == "betweeni" & (this.Linfo == "S" | this.Linfo == "WF" | this.Linfo == "BF" | this.Linfo == "WT" | this.Linfo == "BT") ) |
                
                ( lambdaConstraint == "priorx" & (this.Linfo == "S" | this.Linfo == "WF" | (jjmss < iimss & this.Linfo == "BF") ) ) |
                ( lambdaConstraint == "priori" & (this.Linfo == "S" | this.Linfo == "WF" | this.Linfo == "WT" | (jjmss < iimss & ( this.Linfo == "BF" | this.Linfo == "BT" ) ) ) )
                
                ) {
                    
                    
                    rndTrys <- mxHatLambda[iimss, jjmss] + c( 0, rnorm(nsearch, 0, mss.sd[2]) )
                    
                    if( run.parallel ) {
                        sfOut <- sfClusterApplyLB( x=1:length(rndTrys), fun=ikcirt.fun.mss.lambda, iimss=iimss, jjmss=jjmss,
                        rndTrys=rndTrys,
                        mxHatLambda=mxHatLambda,
                        penalty=penalty, usetruesigma=usetruesigma )
                        
                    }
                    cost.vec <- unlist(sfOut) ; cost.vec
                    cost.vec[is.na(cost.vec)] <- Inf
                    cost.vec[is.nan(cost.vec)] <- Inf
                    
                    xndx.mincost <- which.min(cost.vec)
                    best.cost <- cost.vec[ xndx.mincost ]
                    
                    mxHatLambda[iimss, jjmss] <- rndTrys[xndx.mincost]
                    
                    #cat(best.cost, "\n")
                }
                
                txtProgressBar(min = 0, max = nrow(mxHatLambda), initial = iimss, char = "=", width = 100, style=3)
            }
        }
        cat("\n\n")
    }
    
    ######################################## MSS hatEta search
    
    mxHatDLS <- mxDelta %*% mxHatLambda %*% mxSlot
    mxHatLS <- mxHatLambda %*% mxSlot
    
    sfExport("mxHatDLS", "mxHatLambda", "mxHatLS")
    
    iimss <- 1 ; jjmss <- 1
    
    if(mss.sd[3] != 0) {
        
        cat("Locating Eta ...", "\n")
        for(iimss in 1:nrow(mxHatEta)) {
            for(jjmss in 1:ncol(mxHatEta)) {
                
                rndTrys <- mxHatEta[iimss, jjmss] + c( 0, rnorm(nsearch, 0, mss.sd[3]) )
                
                if( run.parallel ) {
                    sfOut <- sfClusterApplyLB( x=1:length(rndTrys), fun=ikcirt.fun.mss.eta, iimss=iimss, jjmss=jjmss,
                    rndTrys=rndTrys,
                    mxHatEta=mxHatEta,
                    penalty=penalty, usetruesigma=usetruesigma )
                    
                }
                
                cost.vec <- unlist(sfOut)
                cost.vec[is.na(cost.vec)] <- Inf
                cost.vec[is.nan(cost.vec)] <- Inf
                
                xndx.mincost <- which.min(cost.vec)
                best.cost <- cost.vec[ xndx.mincost ]
                
                mxHatEta[iimss, jjmss] <- rndTrys[xndx.mincost]
                
                #cat(iimss, jjmss, best.cost, "\n")
                
                
                txtProgressBar(min = 0, max = nrow(mxHatEta), initial = iimss, char = "=", width = 100, style=3)
                
            }
            #txtProgressBar(min = 0, max = length(hatMu), initial = iimss, char = "=", width = 100, style=3)
        }
        
        cat("\n\n")
        
    }
    
    if(sfIsRunning()) { sfStop() }
    
    for(i in 1:model$nuc) {
        #mxHatEta[ , i] <- qnorm( (rank(mxHatEta[ , i])-0.5) / nrow(mxHatEta) )
    }
    
    
    
    
    ###################### calculate final system cost
    mxHatLSE <- mxHatLambda %*% mxSlot %*% t(mxHatEta)
    
    
    mxDhLS <- mxDelta %*% mxHatLambda %*% mxSlot
    hatSysCov <- mxDhLS %*% var(mxHatEta) %*% t(mxDhLS)   +   covStochastic
    
    if(usetruesigma) {
        useSysCov <- mxSigma
    } else {
        
        useSysCov <- hatSysCov
    }
    
    
    
    
    
    
    
    hatZstar <- mxDelta %*% ( hatMu  +  mxHatLSE )
    
    if(penalty == "logit") {
        hatWstar <- hatZstar / sqrt(diag(covStochastic))
        hatY <- matrix( pnorm( hatWstar ), nrow(hatWstar), ncol(hatWstar) )
        our.cost <- - mean( model$Y * log(hatY) + (1-model$Y) * log(1-hatY), na.rm=TRUE ) ; our.cost
    }
    if(penalty == "L2") {
        hatWstar <- hatZstar / sqrt(diag(useSysCov))
        our.cost <- sqrt( mean( (Z - hatWstar)^2, na.rm=TRUE ) ) ; our.cost
    }
    if(penalty == "L2c") {
        hatWstar <- varZ %*% hatZstar / sqrt(diag(useSysCov))
        our.cost <- sqrt( mean( (Zconv - hatWstar)^2, na.rm=TRUE ) ) ; our.cost
    }
    if(penalty == "miscat") {
        our.cost <- 1 - sum( (hatZstar > 0 & Y == 1) | (hatZstar <= 0 & Y == 0) , na.rm=TRUE ) / sum(!is.na(Y)) ; our.cost
        
    }
    
    
    
    cat("ending", penalty, "cost:", our.cost, "\n")
    
    
    
    if(penalty=="logit") {
        model[["LogitResErr"]] <- our.cost
        model[["L2resErr"]] <- NA
    }
    if(penalty=="L2") {
        model[["LogitResErr"]] <- NA
        model[["L2resErr"]] <- our.cost
    }
    
    
    
    
    model$hatMu <- hatMu
    model$mxHatLambda <- mxHatLambda
    model$mxHatEta <- mxHatEta
    model$hatSysCov <- hatSysCov
    
    if(!is.null(model$mxEta)) {
        performance <- rep(NA, model$nuc)
        for(i in 1:model$nuc) {
            performance[i] <- cor(model$mxEta[ ,i], mxHatEta[ ,i])^2
        }
        model[["performance"]] <- performance
    }
    
    return(model)
    
}
