################################################################################
# Regression Deletion Diagnostics for manylm objects                           #
################################################################################

covratio.manylm <-
function (model, infl = manylm.influence(model, do.coef = FALSE),
    res = weighted.residuals(model)) {
    
    n <- nrow(model$qr$qr)
    p <- model$rank
    omh <- 1 - infl$hat
    e.star <- res/(infl$sigma * sqrt(omh))
    e.star[is.infinite(e.star)] <- NaN
    1/(omh * (((n - p - 1) + e.star^2)/(n - p))^p)
    
}

# setGeneric("covratio")   
#setMethod("covratio", "manylm", covratio.manylm)


