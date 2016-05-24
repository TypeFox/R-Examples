components_ashape3d <-
function (as3d, indexAlpha = 1) 
{
    if (class(as3d) != "ashape3d") {
        cat("Argument is not of class ashape3d.\n")
        return(invisible())
    }
    components = NULL
    if (class(indexAlpha) == "character") 
        if (indexAlpha == "ALL" | indexAlpha == "all") 
            indexAlpha = 1:length(as3d$alpha)
    if (any(indexAlpha > length(as3d$alpha)) | any(indexAlpha <= 
        0)) {
        if (max(indexAlpha) > length(as3d$alpha)) 
            error = max(indexAlpha)
        else error = min(indexAlpha)
        stop(paste("indexAlpha out of bound : valid range = 1:", 
            length(as3d$alpha), ", problematic value = ", error, 
            sep = ""), call. = TRUE)
    }
    for (iAlpha in indexAlpha) {
        nbPoints = dim(as3d$x)[1]
        edgeSelectResult = .C("edgeselect", n = as.integer(nbPoints), 
            x = as.double(as3d$x), ed = as.integer(as3d$edge[, 
                1:2]), nb = as.integer(dim(as3d$edge)[1]), alpha = as.double(as3d$alpha[iAlpha]), 
            nbfinal = as.integer(0), PACKAGE = "alphashape3d")
        edges = as3d$edge[, 1:2][edgeSelectResult$ed[1:edgeSelectResult$nbfinal], 
            ]
        adjacency = rep(list(NULL), nbPoints)
        if (dim(edges)[1] > 0) 
            for (ij in 1:dim(edges)[1]) {
                adjacency[[edges[ij, 1]]] = append(adjacency[[edges[ij, 
                  1]]], edges[ij, 2])
                adjacency[[edges[ij, 2]]] = append(adjacency[[edges[ij, 
                  2]]], edges[ij, 1])
            }
        visited = rep(0, nbPoints)
        for (ij in 1:nbPoints) {
            if (is.null(adjacency[[ij]])) 
                visited[ij] = -1
        }
        nbPart = 0
        while (!all(as.logical(visited))) {
            nbPart = nbPart + 1
            toVisit = which(visited == 0)[1]
            while (!length(toVisit) == 0) {
                cour = toVisit[1]
                neighboor = adjacency[[cour]]
                visited[cour] = nbPart
                for (ii in neighboor) {
                  if (visited[ii] == 0) {
                    toVisit = c(toVisit, ii)
                    visited[ii] = nbPart
                  }
                }
                toVisit = toVisit[-1]
            }
        }
        label = matrix(as.numeric(as.matrix(as.data.frame(table(visited)))), 
            ncol = 2)
        if (any(label[, 1] == -1)) 
            label = label[-which(label[, 1] == -1), , drop = FALSE]
        if (dim(label)[1] > 1) {
            label = label[order(label[, 2], decreasing = TRUE), 
                ]
            visited[visited == -1] = NA
            for (ii in 1:dim(label)[1]) visited[visited == label[ii, 
                1]] = -ii
            visited[is.na(visited)] = 1
            visited = -visited
        }
        if (length(indexAlpha) > 1) {
            components = c(components, list(visited))
        }
        else {
            components = visited
        }
    }
    return(components)
}
