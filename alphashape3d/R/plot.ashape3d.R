plot.ashape3d <-
function (x, clear = TRUE, col = c(2, 2, 2), byComponents = FALSE, 
    indexAlpha = 1, transparency = 1, walpha = FALSE, triangles = TRUE, 
    edges = TRUE, vertices = TRUE, ...) 
{
    arg.triangles <- triangles
    arg.edges <- edges
    arg.vertices <- vertices
    as3d <- x
    triangles <- as3d$triang
    edges <- as3d$edge
    vertex <- as3d$vertex
    x <- as3d$x
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
    if (clear) {
        rgl.clear()
    }
    if (byComponents) {
        components = components_ashape3d(as3d, indexAlpha)
        if (length(indexAlpha) == 1) 
            components = list(components)
        indexComponents = 0
        for (iAlpha in indexAlpha) {
            if (iAlpha != indexAlpha[1]) 
                rgl.open()
            if (walpha) 
                title3d(main = paste("alpha =", as3d$alpha[iAlpha]))
            cat("Device ", rgl.cur(), " : alpha = ", as3d$alpha[iAlpha], 
                "\n")
            indexComponents = indexComponents + 1
            components[[indexComponents]][components[[indexComponents]] == 
                -1] = 0
            colors = c("#000000", sample(rainbow(max(components[[indexComponents]]))))
            if (arg.triangles) {
                tr <- t(triangles[triangles[, 8 + iAlpha] == 
                  2 | triangles[, 8 + iAlpha] == 3, c("tr1", 
                  "tr2", "tr3")])
                if (length(tr) != 0) 
                  rgl.triangles(x[tr, 1], x[tr, 2], x[tr, 3], 
                    col = colors[1 + components[[indexComponents]][tr]], 
                    alpha = transparency, ...)
            }
            if (arg.edges) {
                ed <- t(edges[edges[, 7 + iAlpha] == 2 | edges[, 
                  7 + iAlpha] == 3, c("ed1", "ed2")])
                if (length(ed) != 0) 
                  rgl.lines(x[ed, 1], x[ed, 2], x[ed, 3], col = colors[1 + 
                    components[[indexComponents]][ed]], alpha = transparency, 
                    ...)
            }
            if (arg.vertices) {
                vt <- t(vertex[vertex[, 4 + iAlpha] == 2 | vertex[, 
                  4 + iAlpha] == 3, "v1"])
                if (length(vt) != 0) 
                  rgl.points(x[vt, 1], x[vt, 2], x[vt, 3], col = colors[1 + 
                    components[[indexComponents]][vt]], alpha = transparency, 
                    ...)
            }
        }
    }
    else {
        for (iAlpha in indexAlpha) {
            if (iAlpha != indexAlpha[1]) 
                rgl.open()
            if (walpha) 
                title3d(main = paste("alpha =", as3d$alpha[iAlpha]))
            cat("Device ", rgl.cur(), " : alpha = ", as3d$alpha[iAlpha], 
                "\n")
            if (arg.triangles) {
                tr <- t(triangles[triangles[, 8 + iAlpha] == 
                  2 | triangles[, 8 + iAlpha] == 3, c("tr1", 
                  "tr2", "tr3")])
                if (length(tr) != 0) 
                  rgl.triangles(x[tr, 1], x[tr, 2], x[tr, 3], 
                    col = col[1], , alpha = transparency, ...)
            }
            if (arg.edges) {
                ed <- t(edges[edges[, 7 + iAlpha] == 2 | edges[, 
                  7 + iAlpha] == 3, c("ed1", "ed2")])
                if (length(ed) != 0) 
                  rgl.lines(x[ed, 1], x[ed, 2], x[ed, 3], col = col[2], 
                    alpha = transparency, ...)
            }
            if (arg.vertices) {
                vt <- t(vertex[vertex[, 4 + iAlpha] == 2 | vertex[, 
                  4 + iAlpha] == 3, "v1"])
                if (length(vt) != 0) 
                  rgl.points(x[vt, 1], x[vt, 2], x[vt, 3], col = col[3], 
                    alpha = transparency, ...)
            }
        }
    }
}
