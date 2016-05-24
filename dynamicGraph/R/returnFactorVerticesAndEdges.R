"returnFactorVerticesAndEdges" <-
function (Vertices, factors = NULL, types = "Generator", 
          factorVertexColor = "default", 
          factorEdgeColor = "DarkOliveGreen", 
          fixedFactorPositions = FALSE, 
          factorClasses = validFactorClasses()) 
{
    "newFactorVertexList" <- function(list) 
        return(new("dg.FactorVertexList", nodeList = list))
    vertex.names <- Names(Vertices)
    "subReturnFactorList" <- function(factors, vertices, offset = 0, 
                                      types = "Generator", width = 2, 
                                      color = factorVertexColor) {
        if ((length(types) == 1) && !(is.list(types))) {
            Type <- grep(types, paste(factorClasses[, 1]))
            type <- paste(factorClasses[, 1][Type])
        }
        n <- length(factors)
        FactorVertices <- vector("list", n)
        FactorEdges <- NULL
        PairEdges <- NULL
        for (i in seq(along = factors)) {
            edge <- factors[[i]]
            if (!is.numeric(edge)) 
                edge <- match(edge, Names(vertices))
            m <- length(edge)
            edge.vertices <- vector("list", m)
            for (j in seq(along = edge)) if (edge[j] > 0) {
                edge.vertices[[j]] <- vertices[[edge[j]]]
                FactorEdges <- rbind(FactorEdges, c(-i - offset, edge[j]))
                if (j < length(edge)) 
                  for (k in (j + 1):length(edge)) if (edge[k] > 0) 
                    if (i == 1) 
                      PairEdges <- rbind(PairEdges, c(edge[j], edge[k]))
                    else {
                      x <- c(edge[j], edge[k])
                      if (!is.null(PairEdges)) 
                        if (!any(apply(PairEdges, 1, 
                          function(i) all(i == x)))) 
                          PairEdges <- rbind(PairEdges, x)
                    }
            }
            if (is.list(types)) {
                Type <- types[[i]] == paste(factorClasses[, 2])
                if (!any(Type)) 
                  Type <- types[[i]] == paste(factorClasses[, 1])
                type <- paste(factorClasses[, 1][Type])
            }
            else if (length(types) > 1) {
                Type <- types[i] == paste(factorClasses[, 2])
                if (!any(Type)) 
                  Type <- types[i] == paste(factorClasses[, 1])
                type <- paste(factorClasses[, 1][Type])
            }
            else {
                type <- types
            }
            class(edge.vertices) <- "dg.VertexList"
            prototype <- "dg.Generator"
            x <- match(type, factorClasses[, 1])
            if (!is.null(x)) 
                prototype <- paste(factorClasses[, 2][x])
            if (color == "default") 
                color <- c("yellow", "cyan", "magenta", "blue")[x]
            FactorVertices[[i]] <- new(prototype, vertex.indices = edge, 
                                       vertices = edge.vertices, 
                                       index = -i - offset, color = color, 
                                       fixed.positions = fixedFactorPositions)
        }
        class(FactorVertices) <- "dg.FactorVertexList"
        if (is.null(names(factors))) 
            names(FactorVertices) <- Names(FactorVertices)
        else names(FactorVertices) <- names(factors)
        return(list(FactorVertices = FactorVertices, 
                    FactorEdges = FactorEdges, 
                    PairEdges = PairEdges))
    }
    "two.to.pairs" <- function(from, to) {
        edge.list <- vector("list", length(to))
        for (j in seq(along = to)) edge.list[[j]] <- c(from[j], to[j])
        return(edge.list)
    }
    FactorVertices <- NULL
    PairEdges <- NULL
    FactorEdges <- NULL
    if (!(is.null(factors)) && (length(factors) > 0)) {
        if (!is.list(factors[[1]])) {
            result <- subReturnFactorList(factors, Vertices, 
                                          types = types, width = 2, 
                                          color = factorVertexColor)
            FactorVertices <- result$FactorVertices
            PairEdges <- result$PairEdges
            FactorEdges <- returnFactorEdgeList(
                                      two.to.pairs(result$FactorEdges[, 1], 
                                                   result$FactorEdges[, 2]), 
                                      Vertices, FactorVertices, 
                                      color = factorEdgeColor)
        }
        else {
            for (j in seq(along = factors)) {
                result <- subReturnFactorList(factors[[j]], Vertices, 
                                              offset = length(FactorVertices), 
                                              types = names(factors)[j], 
                                              width = 2, 
                                              color = factorVertexColor)
                FactorVertices <- append(FactorVertices, result$FactorVertices)
                if (is.null(PairEdges)) 
                  PairEdges <- result$PairEdges
                else for (i in 1:nrow(result$PairEdges)) {
                  x <- result$PairEdges[i, ]
                  if (!any(apply(PairEdges, 1, function(i) all(i == x)))) 
                    PairEdges <- rbind(PairEdges, x)
                }
                FactorEdges <- append(FactorEdges, 
                  returnFactorEdgeList(two.to.pairs(result$FactorEdges[, 1], 
                                                    result$FactorEdges[, 2]), 
                  Vertices, FactorVertices, color = factorEdgeColor))
            }
        }
    }
    class(FactorVertices) <- "dg.FactorVertexList"
    class(FactorEdges) <- "dg.FactorEdgeList"
    return(list(FactorVertices = FactorVertices, 
                FactorEdges = FactorEdges, 
                PairEdges = PairEdges))
}
