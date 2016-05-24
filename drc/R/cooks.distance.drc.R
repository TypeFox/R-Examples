cooks.distance.drc <- function(model, ...)
{
    hatVal <- hatvalues(model)
    
    mse <- (rse(model))^2
    if (is.na(mse)) {mse <- 1}  
    # default for generalized linear models (assuming no overdispersion!)
    
    (residuals(model)^2) * (hatVal/((1-hatVal)^2)) / ((length(hatVal)- df.residual(model)) * mse)
}
