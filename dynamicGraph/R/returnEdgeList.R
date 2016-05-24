"returnEdgeList" <-
function (edge.list, vertices, width = 2, color = "DarkSlateGrey", 
    N = 3, oriented = NA, types = NULL, edgeClasses = validEdgeClasses()) 
{
    "newVertexEdgeList" <- function(list) return(new("dg.VertexEdgeList", 
        nodeList = list, N = N))
    vertex.names <- Names(vertices)
    if (is.null(edge.list)) 
        edge.list <- vector("list", length = 0)
    n <- length(edge.list)
    if (!is.null(oriented) && (length(oriented) == 1)) 
        oriented <- rep(oriented, n)
    if (!is.null(oriented) && !(length(oriented) == n)) 
        warning("Invalid length of argument 'oriented'")
    if (!is.null(width) && (length(width) == 1)) 
        width <- rep(width, n)
    if (!is.null(width) && !(length(width) == n)) 
        warning("Invalid length of argument 'width'")
    if (!is.null(color) && (length(color) == 1)) 
        color <- rep(color, n)
    if (!is.null(color) && !(length(color) == n)) 
        warning("Invalid length of argument 'color'")
    if (!.IsEmpty(types) && !(length(types) == n)) 
        warning("Invalid length of argument 'types'")
    if (n == 0) 
        result <- new("dg.VertexEdgeList")
    else {
        result <- vector("list", n)
        for (i in seq(along = edge.list)) {
            edge <- edge.list[[i]]
            if (!is.numeric(edge)) 
                edge <- match(edge, vertex.names)
            m <- length(edge)
            edge.vertices <- vector("list", m)
            for (j in seq(along = edge)) edge.vertices[[j]] <- vertices[[edge[j]]]
            if (is.null(types)) 
                type <- edgeClasses[, 1][1]
            else type <- types[i]
            class(edge.vertices) <- "dg.VertexList"
            if (!is.na(type) && type == "VertexEdge") 
                prototype <- "dg.VertexEdge"
            else prototype <- typeToPrototype(type = type, prototype = "dg.VertexEdge", 
                classes = edgeClasses)
            result[[i]] <- new(prototype, vertex.indices = edge, 
                vertices = edge.vertices, width = width[i], color = color[i], 
                oriented = oriented[i], N = N)
        }
        class(result) <- "dg.VertexEdgeList"
        names(result) <- Labels(result)
    }
    return(result)
}
