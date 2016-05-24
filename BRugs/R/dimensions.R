"dimensions" <-
function(node)
#   Get dimension information for quantity in OpenBUGS model
{
    nodeLabel <- as.character(node)
    if(!(nodeLabel %in% modelNames()))
        stop("node must be a variable name from the model")
    dimensions <- .OpenBUGS(c("BugsRobjects.SetVariable", "BugsRobjects.GetNumDimensions"),
                            c("CharArray", "Integer"),
                            list(nodeLabel, NA))[[2]]
    dimensions
}
