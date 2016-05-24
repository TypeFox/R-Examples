canonicalizeAndMergeScene <- function (scene, ...){
    which <- list(...)
    if (is.Triangles3D(scene)) {
        n.tri <- nrow(scene$v1)
        for (n in which) if (length(scene[[n]]) != n.tri) 
            scene[[n]] <- rep(scene[[n]], length = n.tri)
        scene
    }
    else {
        scene <- lapply(scene, canonicalizeAndMergeScene, ...)
        x <- scene[[1]]
        x$v1 <- do.call(rbind, lapply(scene, function(x) x$v1))
        x$v2 <- do.call(rbind, lapply(scene, function(x) x$v2))
        x$v3 <- do.call(rbind, lapply(scene, function(x) x$v3))
        for (n in which) x[[n]] <- do.call(c, lapply(scene, function(x) x[[n]]))
        x
    }
}


colorScene <- function (scene) {
    if (is.Triangles3D(scene)) 
        colorTriangles(scene)
    else lapply(scene, colorTriangles)
}


t2ve <- function (triangles)  {
    vb <- rbind(triangles$v1, triangles$v2, triangles$v3)
    vbmin <- min(vb)
    vbmax <- max(vb)
    S <- 10^5
    score <- function(v, d) floor(as.vector(v %*% d))
    scale <- function(v) (1 - 1/S) * (v - vbmin)/(vbmax - vbmin)
    d <- c(S, S^2, S^3)
    scores <- score(scale(vb), d)
    vb <- vb[!duplicated(scores), ]
    scores <- score(scale(vb), d)
    ib <- rbind(match(score(scale(triangles$v1), d), scores), 
        match(score(scale(triangles$v2), d), scores), match(score(scale(triangles$v3), 
            d), scores))
    list(vb = t(vb), ib = ib)
}

is.Triangles3D <- function (x) identical(class(x), "Triangles3D")


colorTriangles <- function (triangles){
    if (is.function(triangles$color) || is.function(triangles$color2)) {
        v <- (triangles$v1 + triangles$v2 + triangles$v3)/3
        if (is.function(triangles$color)) 
            triangles$color <- triangles$color(v[, 1], v[, 2], 
                v[, 3])
        if (is.function(triangles$color2)) 
            triangles$color2 <- triangles$color2(v[, 1], v[, 
                2], v[, 3])
        if (is.function(triangles$col.mesh)) 
            triangles$col.mesh <- triangles$col.mesh(v[, 1], 
                v[, 2], v[, 3])
    }
    triangles
}
