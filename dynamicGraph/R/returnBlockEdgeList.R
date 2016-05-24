"returnBlockEdgeList" <-
function (edge.list, vertices, blocks, visibleBlocks = 1:length(blocks), 
    width = 2, color = "default", N = 3, oriented = NA, type = NULL) 
{
    "newBlockEdgeList" <- function(list) return(new("dg.BlockEdgeList", 
        nodeList = list))
    "which.edge" <- function(e) 
        unlist(lapply(result, function(i) all(i@vertex.indices == e)))
    vertex.names <- Names(vertices)
    if (is.null(edge.list)) 
        edge.list <- vector("list", length = 0)
    result <- new("dg.BlockEdgeList")
    if (!is.null(oriented) && !(length(oriented) == 1)) 
        warning("Invalid length of argument 'oriented'")
    if (!is.null(width) && !(length(width) == 1)) 
        warning("Invalid length of argument 'width'")
    if (!is.null(color) && !(length(color) <= 2)) 
        warning("Invalid length of argument 'color'")
    if (!is.null(blocks) && !is.null(blocks[[1]])) {
        for (i in seq(along = edge.list)) {
            edge <- edge.list[[i]]
            if (!is.numeric(edge)) 
                edge <- match(edge, vertex.names)
            strata <- unlist(lapply(edge, function(i) stratum(vertices[[i]])))
            blockindex <- unlist(lapply(edge, function(i) 
                                                   blockindex(vertices[[i]])))
            if ((blockindex[1] != blockindex[2]) && (length(strata) > 
                1) && !any(is.na(strata))) {
                if (strata[1] < strata[2]) {
                  b1 <- blockindex[1]
                  b2 <- blockindex[2]
                  e1 <- edge[1]
                  e2 <- edge[2]
                }
                else {
                  b1 <- blockindex[2]
                  b2 <- blockindex[1]
                  e1 <- edge[2]
                  e2 <- edge[1]
                }
                "f" <- function(x) if (!(any(which.edge(x)))) {
                  block.vertices <- lapply(x, function(i) if (i < 0) 
                    blocks[[-i]]
                  else vertices[[i]])
                  if (color == "default") 
                    color <- ifelse(all(x < 0), "DarkGreen", 
                      "DarkBlue")
                  if (length(color) > 1) 
                    color <- ifelse(all(x < 0), color[1], color[2])
                  class(block.vertices) <- "dg.NodeList"
                  result <<- append(result, list(new("dg.BlockEdge", 
                    vertex.indices = x, vertices = block.vertices, 
                    width = width, color = color, oriented = oriented, 
                    N = N)))
                }
                b1.plus.ancestors <- numeric(0)
                if (b1 != 0) 
                  b1.plus.ancestors <- c(b1, blocks[[b1]]@ancestors)
                b2.plus.ancestors <- numeric(0)
                if (b2 != 0) 
                  b2.plus.ancestors <- c(b2, blocks[[b2]]@ancestors)
                if (b2 != 0) {
                  f(c(e1, -b2))
                  if (any(blocks[[b2]]@ancestors != 0)) 
                    for (i in blocks[[b2]]@ancestors) 
                      if ((i != 0) && !is.element(i, b1.plus.ancestors)) 
                      f(c(e1, -i))
                }
                if (b1 != 0) {
                  f(c(-b1, e2))
                  if (any(blocks[[b1]]@ancestors != 0)) 
                    for (i in blocks[[b1]]@ancestors) 
                      if ((i != 0) && !is.element(i, b2.plus.ancestors)) 
                      f(c(-i, e2))
                }
                if ((b1 != 0) && (b2 != 0)) {
                  f(c(-b1, -b2))
                  if ((any(blocks[[b1]]@ancestors != 0)) && 
                      (any(blocks[[b2]]@ancestors != 0))) 
                    for (i in b1.plus.ancestors) 
                      if ((i != 0) && !is.element(i, b2.plus.ancestors)) 
                      for (j in b2.plus.ancestors) 
                        if ((j != 0) && !is.element(j, b1.plus.ancestors)) 
                        f(c(-i, -j))
                }
            }
        }
        if (!is.null(result)) {
            class(result) <- "dg.BlockEdgeList"
            names(result) <- Labels(result)
        }
    }
    return(result)
}
