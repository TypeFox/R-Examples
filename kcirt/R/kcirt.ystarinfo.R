kcirt.ystarinfo <-
function(model) {
    
    sysvar <- model$covStochastic
    #   inv.sysvar <- solve(sysvar)
    inv.sysvar <- pseudoinverse(sysvar)
    
    mxDLS <- model$mxDelta %*% model$mxLambda %*% model$mxSlot
    #ii <- 1 ; jj <- 2
    FIM <- matrix(NA, model$nuc, model$nuc)
    for(ii in 1:model$nuc) {
        deriv.ii <- rep(0, model$nuc) ; deriv.ii[ii] <- 1
        for(jj in 1:model$nuc) {
            deriv.jj <- rep(0, model$nuc) ; deriv.jj[jj] <- 1
            FIM[ii,jj] <- t(mxDLS %*% deriv.ii) %*% inv.sysvar %*% (mxDLS %*% deriv.jj)
        }
    }
    
    return(FIM)
    
    
}
