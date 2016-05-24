"validEdgeClasses" <-
function () 
{
    result <- cbind(edgeTypes = c("VertexEdge", "Dashed", "Dotted", 
        "DoubleArrow", "DoubleConnected", "TripleConnected"), 
        edgeClasses = c("dg.VertexEdge", "dg.DashedEdge", "dg.DottedEdge", 
            "dg.DoubleArrowEdge", "dg.DoubleConnectedEdge", "dg.TripleConnectedEdge"))
    dimnames(result) <- list(result[, 1], c("Label", "Class"))
    return(result)
}
