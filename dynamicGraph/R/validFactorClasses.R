"validFactorClasses" <-
function () 
{
    result <- cbind(factorTypes = c("Generator", "Discrete generator", 
        "Linear generator", "Quadratic generator"), factorClasses = c("dg.Generator", 
        "dg.DiscreteGenerator", "dg.LinearGenerator", "dg.QuadraticGenerator"))
    dimnames(result) <- list(result[, 1], c("Label", "Class"))
    return(result)
}
