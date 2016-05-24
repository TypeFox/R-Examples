drawScene.rgl <- function(scene, add = FALSE, ...) {
    loadRGL()
    if (! rgl.cur())
        open3d()
    if (!add)
        clear3d()

    scene <- colorScene(scene)
    triangles <- canonicalizeAndMergeScene(scene, "color", "color2", "alpha",
                                           "col.mesh", "fill", "smooth")
    use.col2 <- ! all(is.na(triangles$color2))
    if (use.col2 && any(is.na(triangles$color2)))
        warning(paste("mixing surfaces with and without color2 may not",
                      "work properly in the rgl engine"))

    col <- rep(triangles$color, each = 3)
    if (use.col2)
    	col2 <- rep(triangles$color2, each=3)
    alpha <- rep(triangles$alpha, each = 3)
    fill <- rep(triangles$fill, each = 3)
    col.mesh <- rep(triangles$col.mesh, each = 3)
    data <- zipTriangles(triangles)
    if (all(fill)) {
        front <- "filled"
        back <- "filled"
    }
    else if (any(fill))
        ##**** handle these by splitting; OK if no alpha < 1
        stop(paste("for now rgl engine cannot handle mixed fill/wire",
                   "frame surfaces"))
    else {
        front <- "lines"
        back <- "lines"
        col <- col.mesh
    }
    oldstyle = getOption("old.misc3d.orientation")
    if (! is.null(oldstyle) && oldstyle) {
        data <- data[,c(1, 3, 2)]
        data[,3] <- -data[,3]
    }
    if (any(triangles$smooth > 0)) {
        if (any(triangles$smooth == 0))
            stop(paste("for now for the rgl engine cannot handle mixed",
                       "smooth/non-smooth surfaces"))
        normals <- zipTriangles(triangleVertexNormals(triangles))
    }
    else normals <- NULL
    if (nrow(data) > 0) # to avoid a segfault in rgl
    {
        if (! use.col2)
    	    triangles3d(data[,1], data[,2], data[,3],
                      col = col, alpha = alpha, normals = normals,
                      front = front, back = back, ...)
        else {
            triangles3d(data[,1], data[,2], data[,3],
	              col = col, alpha = alpha, normals = normals,
                      front = front, back = "cull", ...)
            triangles3d(data[,1], data[,2], data[,3],
	              col = col2, alpha = alpha, normals = normals,
                      front = "cull", back = back, ...)
        }
    }
        
}

renderScene <- function(scene, box, fill, col.mesh, add, engine, polynum,
                        col.bg, depth, newpage) {
    triangles <- canonicalizeAndMergeScene(scene, "color", "col.light",
                                           "col.mesh", "fill")
    v1 <- triangles$v1
    v2 <- triangles$v2
    v3 <- triangles$v3
    n.tri <- nrow(v1)
    if (fill) {
        fill <- rep(triangles$fill, length = n.tri)
        col.mesh <- rep(triangles$col.mesh, length = n.tri)
    }
    else col.mesh <- rep(col.mesh, length = n.tri)
    col.fill <- ifelse(fill, triangles$col.light, NA)
    z <- (v1[,3] + v2[,3] + v3[,3]) / 3
    if (depth > 0) {
        rgbcol <- col2rgb(col.fill, alpha = TRUE)/255
        rgbcol.bg <- col2rgb(col.bg, alpha = TRUE)/255
        s <- (1 + depth * z) / (1 + depth * max(z))
        col.fill <- rgb(rgbcol[1,] * s + rgbcol.bg[1,] * (1 - s),
                        rgbcol[2,] * s + rgbcol.bg[2,] * (1 - s),
                        rgbcol[3,] * s + rgbcol.bg[3,] * (1 - s),
                        rgbcol[4,])
    }
    col.mesh <- ifelse(is.na(col.mesh), col.fill, col.mesh)
    i <- order(z, na.last = NA)
    if (engine == "grid") 
        render.grid(v1[i,], v2[i,], v3[i,], box, fill[i], col.fill[i],
                    col.mesh[i], add, polynum, newpage)
    else render.standard(v1[i,], v2[i,], v3[i,], box, fill[i], col.fill[i],
                         col.mesh[i], add)
}

render.standard <- function(v1, v2, v3, box, fill, col.fill, col.mesh, add) {
    if (! add) {
        # rr <- screenRange(v1, v2, v3)
        rr <- range(box)
        plot(rr, rr,type="n", axes = FALSE, ann = FALSE)
    }
    xx <- as.vector(rbind(v1[,1], v2[,1], v3[,1], NA))
    yy <- as.vector(rbind(v1[,2], v2[,2], v3[,2], NA))
    polygon(xx, yy, col=col.fill, border=col.mesh)
}

render.grid <- function(v1, v2, v3, box, fill, col.fill, col.mesh,
                        add, polynum, newpage) {
    if (! add) {
        if (newpage) grid::grid.newpage()
        rr <- range(box)
        grid::pushViewport(grid::viewport(w = 0.8, h = 0.8,
                                          xscale = rr, yscale = rr,
                                          name = "misc3dScene"))
        on.exit(grid::upViewport())
    }

    xx <- as.vector(rbind(v1[,1], v2[,1], v3[,1]))
    yy <- as.vector(rbind(v1[,2], v2[,2], v3[,2]))
    n.tri <- nrow(v1)
    idlen <- rep(3, n.tri)
    start <- 1
    end <- start + polynum - 1
    while (start <= n.tri) {
        end <- min(end, n.tri)
        j <- start : end
        j3 <- (3 * start - 2) : (3 * end)
        gp <- grid::gpar(fill = col.fill[j], col = col.mesh[j])
        grid::grid.polygon(x = xx[j3], y = yy[j3], default.units = "native",
                           gp = gp, id.lengths = idlen[j])
        start <- start + polynum
        end <- start + polynum
    }
}

makePerspMatrix <- function(d) {
    rbind(c(1, 0,  0, 0),
          c(0, 1,  0, 0),
          c(0, 0,  1, 0),
          c(0, 0, -1, 1 / d))
}

## drawScene is a simple function for plotting triangles.  The viewer is
## looking down the positive Z axis.
## The returned value is suitable for use with trans3d.

drawScene <- function(scene, light = c(0, 0, 1),
                      screen = list(z = 40, x = -60),
                      scale = TRUE,
                      R.mat = diag(4),
                      perspective = FALSE,
                      distance = if (perspective) 0.2 else 0, 
                      fill = TRUE,
                      xlim = NULL, ylim = NULL, zlim = NULL,
                      aspect = c(1, 1),
                      col.mesh = if (fill) NA else "black",
                      polynum = 100,
                      lighting = phongLighting,
                      add = FALSE,
                      engine = "standard",
                      col.bg = "transparent", depth = 0, newpage = TRUE) {
    scene <- colorScene(scene)
    sr <- sceneRanges(scene, xlim, ylim, zlim)
    if (add)
        rot.mat <- R.mat
    else
        rot.mat <- makeViewTransform(sr, scale, aspect, screen, R.mat)
    scene <- transformScene(scene, rot.mat)
    scene <- lightScene(scene, lighting, light)
    if (distance > 0) {
        scene <- addPerspective(scene, distance)
        rot.mat <- makePerspMatrix(distance) %*% rot.mat
    }
    box <- as.matrix(expand.grid(sr$xlim, sr$ylim, sr$zlim))
    box <- trans3dto3d(box, rot.mat)
    renderScene(scene, box, fill, col.mesh, add, engine, polynum,
                col.bg, depth, newpage)
    invisible(t(rot.mat))
}

