"returnFactorEdgeList" <-
function (edge.list, vertices, factorvertices = NULL, width = 2, 
    color = "DarkSlateGrey", N = 3, type = NULL) 
{
    "newFactorEdgeList" <- function(list) return(new("dg.FactorEdgeList", 
        nodeList = list))
    vertex.names <- Names(vertices)
    if (is.null(edge.list)) 
        edge.list <- vector("list", length = 0)
    n <- length(edge.list)
    if (!is.null(width) && (length(width) == 1)) 
        width <- rep(width, n)
    if (!is.null(width) && !(length(width) == n)) 
        warning("Invalid length of argument 'width'")
    if (!is.null(color) && (length(color) == 1)) 
        color <- rep(color, n)
    if (!is.null(color) && !(length(color) == n)) 
        warning("Invalid length of argument 'color'")
    if (n == 0) 
        result <- new("dg.FactorEdgeList")
    else {
        result <- vector("list", n)
        for (i in seq(along = edge.list)) {
            edge <- edge.list[[i]]
            if (!is.numeric(edge)) 
                edge <- match(edge, vertex.names)
            m <- length(edge)
            edge.vertices <- vector("list", m)
            for (j in seq(along = edge)) if (edge[j] > 0) 
                edge.vertices[[j]] <- vertices[[edge[j]]]
            else if (edge[j] < 0) 
                edge.vertices[[j]] <- factorvertices[[-edge[j]]]
            class(edge.vertices) <- "dg.VertexList"
            result[[i]] <- new("dg.FactorEdge", vertex.indices = edge, 
                vertices = edge.vertices, width = width[i], color = color[i], 
                N = N)
        }
        class(result) <- "dg.FactorEdgeList"
        names(result) <- Labels(result)
    }
    return(result)
}
