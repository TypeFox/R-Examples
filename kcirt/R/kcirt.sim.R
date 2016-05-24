kcirt.sim <-
function(model, N, type="Eta") {
    
    mxEta <- rmvnorm(N, rep(0, model$nuc), model$covEta)
    
    mxShocks <- rmvnorm(N, rep(0, sum(model$ns)), model$covShocks)
    
    ##mxEta <- matrix(rnorm(N*model$nuc, 0, 1), N, model$nuc)
    
    ##mxShocks <- matrix(rnorm(N*sum(model$ns), 0, 1), N, sum(model$ns))
    
    
    Ustar <- matrix( model$mu, sum(model$ns), N )    +    model$mxLambda %*%  model$mxSlot %*% t(mxEta)  +  t(mxShocks)
    
    Ystar <- model$mxDelta %*% Ustar
    
    mxData <- ikcirt.Ustar2data(Ustar, model$qTypes, model$mxDelta, model$ns)
    
    Y <- ikcirt.data2Y(mxData, model$mxDelta)
    
    
    
    Z <- 2*Y - 1
    
    Yisna <- is.na(Y)
    
    model[["mxEta"]] <- mxEta
    model[["mxShocks"]] <- mxShocks
    model[["Ustar"]] <- Ustar
    model[["Ystar"]] <- Ystar
    model[["mxData"]] <- mxData
    model[["Yisna"]] <- Yisna
    model[["Y"]] <- Y
    model[["Z"]] <- Z
    
    return(model)
    
}
