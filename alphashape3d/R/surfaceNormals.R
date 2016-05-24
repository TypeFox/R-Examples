surfaceNormals <-
function (x, indexAlpha = 1, display = FALSE, col = 3, scale = 1, 
    ...) 
{
    as3d <- x
    tetra <- as3d$tetra
    triangles <- as3d$triang
    edges <- as3d$edge
    vertex <- as3d$vertex
    x <- as3d$x
    if (class(indexAlpha) == "character" & (indexAlpha == "ALL" | 
        indexAlpha == "all")) 
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
    normals.obj = NULL
    for (iAlpha in indexAlpha) {
        tr <- triangles[triangles[, 8 + iAlpha] == 2, c("tr1", 
            "tr2", "tr3")]
        te <- tetra[tetra[, 5 + iAlpha] == 1, 1:4]
        normMat <- numeric(length(tr))
        middlePoint <- numeric(length(tr))
        retour <- .C("triangleNormals", as.integer(tr), dim(tr)[1], 
            as.integer(te), dim(te)[1], as.numeric(x), dim(x)[1], 
            normMat, middlePoint)
        normMat = matrix(retour[[7]], ncol = 3)
        middlePoint = matrix(retour[[8]], ncol = 3)
        if (display) {
            segment = matrix(rep(0, dim(normMat)[1] * 6), ncol = 3)
            for (ii in 1:dim(normMat)[1]) {
                segment[2 * ii - 1, ] = middlePoint[ii, ]
                segment[2 * ii, ] = middlePoint[ii, ] + scale * 
                  normMat[ii, ]
            }
            rgl.open()
            plot(as3d, indexAlpha = iAlpha)
            segments3d(segment, col = col, alpha = 1, ...)
        }
        normals.obj <- c(normals.obj, list(list(normals = normMat, 
            centers = middlePoint)))
        class(normals.obj[[length(normals.obj)]]) = "normals"
    }
    if (length(indexAlpha) == 1) {
        normals.obj <- normals.obj[[1]]
    }
    else {
        class(normals.obj) <- "normals-List"
    }
    invisible(normals.obj)
}
