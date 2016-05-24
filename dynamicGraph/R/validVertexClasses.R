"validVertexClasses" <-
function () 
{
    result <- cbind(vertexTypes = c("Discrete", "Ordinal", "Continuous"), 
        vertexClasses = c("dg.DiscreteVertex", "dg.OrdinalVertex", 
            "dg.ContinuousVertex"))
    dimnames(result) <- list(result[, 1], c("Label", "Class"))
    return(result)
}
