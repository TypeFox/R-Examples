dffits.manylm <- 
function (model, infl = manylm.influence(model, do.coef = FALSE), 
    res = weighted.residuals(model)) {

    res <- res * sqrt(infl$hat)/(infl$sigma * (1 - infl$hat))
    res[is.infinite(res)] <- NaN

    res
}

# setGeneric("dffits")   
#setMethod("dffits", "manylm", dffits.manylm)
