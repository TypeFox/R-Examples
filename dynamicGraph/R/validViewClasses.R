"validViewClasses" <-
function () 
{
    result <- cbind(vertexTypes = c("Simple", "Factor", "Moral", 
        "Essential"), vertexClasses = c("SimpleDynamicGraphView", 
        "FactorDynamicGraphView", "MoralDynamicGraphView", "EssentialDynamicGraphView"))
    dimnames(result) <- list(result[, 1], c("Label", "Class"))
    return(result)
}
