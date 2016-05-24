## These functions work with collections of n triangles.  A collection of
## triangles is a list with components v1, v2, v3 representing the
## coordinates of the three vertices; each of these components is an n by
## 3 matrix.

makeTriangles <- function(v1, v2, v3,
                          color = "red", color2 = NA, alpha = 1,
                          fill = TRUE, col.mesh = if (fill) NA else color,
	                  smooth = 0,  material = "default") {
    if (missing(v2) || missing(v3)) {
        if (missing(v2) && missing(v3))
  	    v <- unzipTriangleMatrix(v1)
        else if (missing(v3))
            v <- ve2t(list(vb = v1, ib = v2))
        else stop("unknown form of triangle specification")
	v1 <- v$v1
	v2 <- v$v2
	v3 <- v$v3
    }
    tris <- structure(list(v1 = v1, v2 = v2, v3 = v3,
                           color = color, color2 = color2, fill = fill,
                           material = material, col.mesh = col.mesh,
                           alpha = alpha, smooth = smooth),
                      class = "Triangles3D")
    colorTriangles(tris)
}

is.Triangles3D <- function(x) identical(class(x), "Triangles3D")

updateTriangles <- function(triangles, color, color2, alpha, fill, col.mesh,
                            material, smooth) {
    if (! missing(color)) triangles$color <- color
    if (! missing(color2)) triangles$color2 <- color2
    if (! missing(fill)) triangles$fill <- fill
    if (! missing(col.mesh)) triangles$col.mesh <- col.mesh
    if (! missing(material)) triangles$material <- material
    if (! missing(alpha)) triangles$alpha <- alpha
    if (! missing(smooth)) triangles$smooth <- smooth
    colorTriangles(triangles)
}

#**** This assumes comparable scaling of dimensions
#**** 5 is the largest exponent for S that will work; smaller is OK
t2ve <- function (triangles)
{
    vb <- rbind(triangles$v1, triangles$v2, triangles$v3)
    vbmin <- min(vb)
    vbmax <- max(vb)
    S <- 10^5
    score <- function(v, d) floor(as.vector(v %*% d))
    scale <- function(v) (1 - 1 / S) * (v - vbmin) / (vbmax - vbmin)
    d <- c(S, S^2, S^3)
    scores <- score(scale(vb), d)
    vb <- vb[! duplicated(scores),]
    scores <- score(scale(vb), d)
    ib <- rbind(match(score(scale(triangles$v1), d), scores),
                match(score(scale(triangles$v2), d), scores),
                match(score(scale(triangles$v3), d), scores))
    list(vb = t(vb), ib = ib)
}

ve2t <- function(ve) {
    list (v1 = t(ve$vb[,ve$ib[1,]]),
          v2 = t(ve$vb[,ve$ib[2,]]),
          v3 = t(ve$vb[,ve$ib[3,]]))
}

unzipTriangleMatrix <- function(tris) {
    if (ncol(tris) != 3)
        stop("triangle matrix must have three columns.")
    if (nrow(tris) %% 3 != 0)
        stop("number of rows in triangle matrix must be divisible by 3")
    n <- nrow(tris) / 3
    list(v1 = tris[3 * (1 : n) - 2,],
         v2 = tris[3 * (1 : n) - 1,],
         v3 = tris[3 * (1 : n),])
}

zipTriangles <- function(tris) {
    n <- nrow(tris$v1)
    if (nrow(tris$v2) != n || nrow(tris$v3) != n)
        stop("vertex arrays must have the same number of rows")
    v <- matrix(0, nrow = 3 * n, ncol = 3)
    v[3 * (1 : n) - 2,] <- tris$v1
    v[3 * (1 : n) - 1,] <- tris$v2
    v[3 * (1 : n),] <- tris$v3
    v
}

colorTriangles <- function(triangles) {
    if (is.function(triangles$color) || is.function(triangles$color2)) {
        v <- (triangles$v1 + triangles$v2 + triangles$v3) / 3
        if (is.function(triangles$color))
            triangles$color <- triangles$color(v[,1], v[,2], v[,3])
        if (is.function(triangles$color2))
            triangles$color2 <- triangles$color2(v[,1], v[,2], v[,3])
        if (is.function(triangles$col.mesh))
            triangles$col.mesh <- triangles$col.mesh(v[,1], v[,2], v[,3])
    }
    triangles
}

colorScene <- function(scene) {
    if (is.Triangles3D(scene))
        colorTriangles(scene)
    else lapply(scene, colorTriangles)
}

## **** better to make new triangles including only requested components?
canonicalizeAndMergeScene <- function(scene, ...) {
    which <- list(...)
    if (is.Triangles3D(scene)) {
        n.tri <- nrow(scene$v1)
        for (n in which)
            if (length(scene[[n]]) != n.tri)
                scene[[n]] <- rep(scene[[n]], length = n.tri)
        scene
    }
    else {
        scene <- lapply(scene, canonicalizeAndMergeScene, ...)
        x <- scene[[1]]
        x$v1 <- do.call(rbind, lapply(scene, function(x) x$v1))
        x$v2 <- do.call(rbind, lapply(scene, function(x) x$v2))
        x$v3 <- do.call(rbind, lapply(scene, function(x) x$v3))
        for (n in which)
            x[[n]] <- do.call(c, lapply(scene, function(x) x[[n]]))
        x
    }
}

expandTriangleGrid <- function(x, y) {
    nx <- length(x) - 1
    ny <- length(y) - 1
    A <- c(0, 0)
    B <- c(1, 0)
    C <- c(1, 1)
    D <- c(0, 1)
    g <- expand.grid(x = 1 : nx, y = 1 : ny)
    even <- (g$x + g$y) %% 2 == 0
    gx11 <- ifelse(even, g$x + A[1], g$x + A[1])
    gy11 <- ifelse(even, g$y + A[2], g$y + A[2])
    gx12 <- ifelse(even, g$x + A[1], g$x + B[1])
    gy12 <- ifelse(even, g$y + A[2], g$y + B[2])
    i1 <- rbind(cbind(gx11, gy11), cbind(gx12, gy12))
    gx21 <- ifelse(even, g$x + B[1], g$x + B[1])
    gy21 <- ifelse(even, g$y + B[2], g$y + B[2])
    gx22 <- ifelse(even, g$x + C[1], g$x + C[1])
    gy22 <- ifelse(even, g$y + C[2], g$y + C[2])
    i2 <- rbind(cbind(gx21, gy21), cbind(gx22, gy22))
    gx31 <- ifelse(even, g$x + C[1], g$x + D[1])
    gy31 <- ifelse(even, g$y + C[2], g$y + D[2])
    gx32 <- ifelse(even, g$x + D[1], g$x + D[1])
    gy32 <- ifelse(even, g$y + D[2], g$y + D[2])
    i3 <- rbind(cbind(gx31, gy31), cbind(gx32, gy32))
    v1 <- cbind(x[i1[,1]], y[i1[,2]])
    v2 <- cbind(x[i2[,1]], y[i2[,2]])
    v3 <- cbind(x[i3[,1]], y[i3[,2]])
    list(v1 = v1, v2 = v2, v3 = v3)
}

## adapted from lattice ltransform3dto3d
trans3dto3d <- function (x, R.mat) {
    if (length(x) == 0)
        return(x)
    val <- R.mat %*% rbind(t(x), 1)
    val[1, ] <- val[1, ]/val[4, ]
    val[2, ] <- val[2, ]/val[4, ]
    val[3, ] <- val[3, ]/val[4, ]
    t(val[1:3, , drop = FALSE])
}

transformTriangles <- function(triangles, R) {
    tr <- function(v) trans3dto3d(v, R)
    triangles$v1 <- tr(triangles$v1)
    triangles$v2 <- tr(triangles$v2)
    triangles$v3 <- tr(triangles$v3)
    triangles
}

transformScene <- function(scene, rot.mat) {
    if (is.Triangles3D(scene))
        transformTriangles(scene, rot.mat)
    else lapply(scene, transformTriangles, rot.mat)
}

translateTriangles <- function(triangles, x = 0, y = 0, z = 0) {
    M <- diag(4)
    M[1:3,4] <- c(x, y, z)
    transformTriangles(triangles, M)
}

scaleTriangles <- function(triangles, x = 1, y = x, z = x) {
    M <- diag(c(x, y, z, 1))
    transformTriangles(triangles, M)
}

## triangleNormals computes the normal vectors to a collection of
## triangles as the vector crossprocuct of the direction from v1 to v2
## and the direction from v2 to v3.  The result is an n by 3 matrix of
## unit representing the n unit normal vectors.

triangleNormals <- function(triangles) {
   x <- triangles$v2 - triangles$v1
   y <- triangles$v3 - triangles$v2
   z <- cbind(x[,2]*y[,3] - x[,3]*y[,2],
              x[,3]*y[,1] - x[,1]*y[,3],
              x[,1]*y[,2] - x[,2]*y[,1])
   z / sqrt(rowSums(z^2))
}

# adapted from lattice ltransform3dMatrix
trans3dMat <- function (screen, P = diag(4)) {
    givens4 <- function(i, j, gamma) {
        T <- diag(4)
        cgamma <- cos(gamma)
        sgamma <- sin(gamma)
        T[c(i,j),c(i,j)] <- matrix(c(cgamma, sgamma, -sgamma, cgamma), 2, 2)
        T
    }
    screen.names <- names(screen)
    for (i in seq(along = screen.names)) {
        if (screen.names[i] == "x")
            P <- givens4(2, 3, screen[[i]] * pi/180) %*% P
        else if (screen.names[i] == "y")
            P <- givens4(1, 3, -screen[[i]] * pi/180) %*% P #**** whi negative?
        else if (screen.names[i] == "z")
            P <- givens4(1, 2, screen[[i]] * pi/180) %*% P
    }
    P
}

makeViewTransform <- function(ranges, scale, aspect, screen, R.mat) {
    m <- c(mean(ranges$xlim), mean(ranges$ylim), mean(ranges$zlim))
    s <- 0.5 * c(diff(ranges$xlim), diff(ranges$ylim), diff(ranges$zlim))
    if (! scale) s <- rep(max(s), 3)
    else s <- s / c(1, aspect)
    A <- diag(1 / c(s, 1))
    A[1:3, 4] <- -m / s
    trans3dMat(screen, R.mat %*% A)
}

trianglesRanges <- function(triangles, xlim, ylim, zlim) {
    v1 <- triangles$v1
    v2 <- triangles$v2
    v3 <- triangles$v3
    if (is.null(xlim)) xlim <- range(v1[,1], v2[,1], v3[,1], na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(v1[,2], v2[,2], v3[,2], na.rm = TRUE)
    if (is.null(zlim)) zlim <- range(v1[,3], v2[,3], v3[,3], na.rm = TRUE)
    list(xlim = xlim, ylim = ylim, zlim = zlim)
}

sceneRanges <- function(scene, xlim, ylim, zlim) {
    if (is.Triangles3D(scene))
        trianglesRanges(scene, xlim, ylim, zlim)
    else {
        ranges <- lapply(scene, trianglesRanges, xlim, ylim, zlim)
        list(xlim = range(sapply(ranges,function(x) x$xlim)),
             ylim = range(sapply(ranges,function(x) x$ylim)),
             zlim = range(sapply(ranges,function(x) x$zlim)))
    }
}

addTrianglesPerspective <- function(triangles, distance) {
    pt <- function(v) {
        v[, 1] <- v[, 1] / (1 / distance - v[, 3])
        v[, 2] <- v[, 2] / (1 / distance - v[, 3])
        v
    }
    triangles$v1 <- pt(triangles$v1)
    triangles$v2 <- pt(triangles$v2)
    triangles$v3 <- pt(triangles$v3)
    triangles
}

addPerspective <- function(scene, distance) {
    if (is.Triangles3D(scene))
        addTrianglesPerspective(scene, distance)
    else lapply(scene, addTrianglesPerspective, distance)
}

screenRange <- function(v1, v2, v3)
    range(v1[,1:2], v2[,1:2], v3[,1:2], na.rm = TRUE)

vertexTriangles <- function(ve) {
    n.vert <- ncol(ve$vb)
    ib <- ve$ib
    vt <- function(i) which(ib[1,] == i | ib[2,] == i | ib[3,] == i)
    lapply(1 : n.vert, vt)
}

# faster version
vertexTriangles <- function(ve) {
    n.vert <- ncol(ve$vb)
    val <- vector("list", n.vert)
    ib <- ve$ib
    for (i in 1 : ncol(ib)) {
        val[[ib[1,i]]] <- c(val[[ib[1,i]]], i)
        val[[ib[2,i]]] <- c(val[[ib[2,i]]], i)
        val[[ib[3,i]]] <- c(val[[ib[3,i]]], i)
    }
    val
}


vertexNormals <- function(vt, N) {
    vn <- function(tris) {
        z <- apply(N[tris,,drop = FALSE], 2, mean, na.rm = TRUE);
        z <- z / sqrt(sum(z^2))
        if (any(is.na(z))) c(1,0,0) else z
    }
    t(sapply(vt, vn))
}

# faster version
vertexNormals <- function(vt, N) {
    val <- matrix(0, nrow = length(vt), ncol = 3)
    for (i in seq(along = vt)) {
        Ni <- N[vt[[i]],,drop = FALSE]
        Ni1 <- Ni[,1]
        Ni2 <- Ni[,2]
        Ni3 <- Ni[,3]
        z1 <- if (any(is.na(Ni1))) mean(Ni1, na.rm = TRUE)
              else sum(Ni1) / length(Ni1)
        z2 <- if (any(is.na(Ni2))) mean(Ni2, na.rm = TRUE)
              else sum(Ni2) / length(Ni2)
        z3 <- if (any(is.na(Ni3))) mean(Ni3, na.rm = TRUE)
              else sum(Ni3) / length(Ni3)
        z <- c(z1, z2, z3)
        z <- z / sqrt(sum(z^2))
        val[i,] <- if (any(is.na(z))) c(1,0,0) else z
    }
    val
}

interpolateVertexNormals <- function(VN, ib) {
    z <- (VN[ib[1,],] + VN[ib[2,],] + VN[ib[3,],]) / 3
    z / sqrt(rowSums(z^2))
}

## triangleVertexNormals computes the normals at the vertices by
## averaging the normals of the incident triangles.  This is used by
## the rgl engine.  The result form is chosen so zipTriangles can be
## used on it.
triangleVertexNormals <- function(v) {
    N <- triangleNormals(v)
    ve <- t2ve(v)
    vt <- vertexTriangles(ve)
    VN <- misc3d:::vertexNormals(vt, N)
    list(v1 = VN[ve$ib[1,],], v2 = VN[ve$ib[2,],], v3 = VN[ve$ib[3,],])
}

vertexColors <- function(vt, col) {
    C <- t(col2rgb(col))
    val <- matrix(0, nrow = length(vt), ncol = 3)
    for (i in seq(along = vt)) {
        vti <- vt[[i]]
        nti <- length(vti)
        Ci <- C[vti,,drop = FALSE]
        Ci1 <- Ci[,1]
        Ci2 <- Ci[,2]
        Ci3 <- Ci[,3]
        val[i,] <- c(sum(Ci1), sum(Ci2), sum(Ci3)) / nti
    }
    val
}

interpolateVertexColors <- function(VC, ib) {
    TC <- (VC[ib[1,],] + VC[ib[2,],] + VC[ib[3,],]) / 3
    rgb(TC[,1], TC[,2], TC[,3], maxColorValue = 255)
}

triangleEdges <- function(vb, ib) {
    edges <- cbind(ib[c(1,2),], ib[c(2,3),], ib[c(3,1),])
    swap <- edges[1,] > edges[2,]
    edges[,swap] <- edges[2:1,swap]
    edges[,! duplicated(edges, MARGIN = 2)]
}

# faster version
triangleEdges <- function(vb, ib) {
    n.vert <- ncol(vb)
    edges <- cbind(ib[c(1,2),], ib[c(2,3),], ib[c(3,1),])
    swap <- edges[1,] > edges[2,]
    edges[,swap] <- edges[2:1,swap]
    score <- as.vector(c(1 + n.vert, 1) %*% edges)
    edges[,! duplicated(score)]
}

triangleMidTriangles <- function(vb, ib, VN) {
    n.vert <- ncol(vb)
    edges <- triangleEdges(vb, ib)
    vb <- (vb[,edges[1,]] + vb[,edges[2,]]) / 2
    d <- c(1 + n.vert, 1)
    scores <- as.vector(d %*% edges)
    mpi <- function(a, b) {
        s <- d[1] * pmin(a, b) + d[2] * pmax(a, b)
        match(s, scores)
    }
    mpi1 <- mpi(ib[1,], ib[2,])
    mpi2 <- mpi(ib[2,], ib[3,])
    mpi3 <- mpi(ib[3,], ib[1,])
    ib <- rbind(mpi1, mpi2, mpi3)
    z <- VN[edges[1,],] + VN[edges[2,],]
    z <- z / sqrt(rowSums(z^2))
    list(vb = vb, ib = ib, VN = z)
}

## surfaceTriangles creates a set of triangles for a grid specified by x,
## y and function falues computed with f if f is a function or taken
## from f if f is a matrix.

surfaceTriangles <- function(x, y, f,
                             color = "red", color2 = NA,  alpha = 1,
                             fill = TRUE, col.mesh = if (fill) NA else color,
                             smooth = 0, material = "default") {
    if (is.function(f))
        ff <- function(ix, iy) f(x[ix], y[iy])
    else
        ff <- function(ix, iy) f[ix + length(x) * (iy - 1)]
    i <- expandTriangleGrid(1 : length(x), 1 : length(y))
    i1 <- i$v1
    i2 <- i$v2
    i3 <- i$v3
    v1 <- cbind(x[i1[,1]], y[i1[,2]], ff(i1[,1], i1[,2]))
    v2 <- cbind(x[i2[,1]], y[i2[,2]], ff(i2[,1], i2[,2]))
    v3 <- cbind(x[i3[,1]], y[i3[,2]], ff(i3[,1], i3[,2]))
    na1 <- is.na(v1[,1]) | is.na(v1[,2]) | is.na(v1[,3])
    na2 <- is.na(v2[,1]) | is.na(v2[,2]) | is.na(v2[,3])
    na3 <- is.na(v3[,1]) | is.na(v3[,2]) | is.na(v3[,3])
    nna <- ! (na1 | na2 | na3)
    makeTriangles(v1[nna,], v2[nna,], v3[nna,],
                  color = color, color2 = color2, fill = fill, smooth = smooth,
                  material = material, col.mesh = col.mesh, alpha = alpha)
}

## pointsTetrahedra computes a collection of tetrahedra centered at
## the specified point locations.  This is useful, for example, for
## displaying raw data along with a density contour in a scene
## rendered with standard or grid graphics. Random orientation might
## be useful to avoid strange results at certain lighting angles.

pointsTetrahedra <- function(x, y, z, size = 0.01, color = "black", ...) {
    n <- length(x)
    if (length(y) != n || length(z) != n)
        stop("coordinate vectors must be the same length.")

    ## Create a basic tetrahedron centered at the origin
    a <- sqrt(3) / 2
    b <- 1 / (2 * sqrt(3))
    h <- sqrt(2 / 3)

    mx <- 1 / 2
    my <- (a + b) / 4
    mz <- h / 4

    A <- c(       -mx,    -my,    -mz)
    B <- c(    1 - mx,    -my,    -mz)
    C <- c(1 / 2 - mx, a - my,    -mz)
    D <- c(1 / 2 - mx, b - my, h - mz)

    v1 <- rbind(B, A, B, C)
    v2 <- rbind(A, B, C, A)
    v3 <- rbind(C, D, D, D)

    ## Scale the tetrahedron
    if (length(size) < 3) size <- rep(size, len = 3)
    if (n == 1) s <- diag(size)
    else s <- diag(size * c(diff(range(x)), diff(range(y)), diff(range(z))))
    sv1 <- v1 %*% s
    sv2 <- v2 %*% s
    sv3 <- v3 %*% s

    ## Compute the tetrahedra for the points, taking advantage of recycling
    x4 <- rep(x, each = 4)
    y4 <- rep(y, each = 4)
    z4 <- rep(z, each = 4)

    V1 <- cbind(x4 + sv1[,1], y4 + sv1[,2], z4 + sv1[,3])
    V2 <- cbind(x4 + sv2[,1], y4 + sv2[,2], z4 + sv2[,3])
    V3 <- cbind(x4 + sv3[,1], y4 + sv3[,2], z4 + sv3[,3])

    makeTriangles(V1, V2, V3, color = color, ...)
}

bresenhamLine <- function(x1, y1, z1, x2, y2, z2, delta){

    if (length(delta) < 3) delta <- rep(delta, len = 3)
   
    vertex <- rep(0,3)
    vertex[1] <- x1
    vertex[2] <- y1
    vertex[3] <- z1
    dx <- x2 - x1
    dy <- y2 - y1
    dz <- z2 - z1
    
    x_inc <- ifelse(dx < 0, -delta, delta)
    l <- abs(dx)/delta[1]
    y_inc <- ifelse(dy < 0, -delta, delta)
    m <- abs(dy)/delta[2]
    z_inc <- ifelse(dz < 0, -delta, delta)
    n <- abs(dz)/delta[3]
    
    dx2 <- 2*l 
    dy2 <- 2*m 
    dz2 <- 2*n

    if ((l >= m) && (l >= n)){
        err_1 <- dy2 - l
        err_2 <- dz2 - l
        Mat <- matrix(0, ncol=3, nrow=l+1)
        ii <- 1
        for (i in 1:l){
            Mat[ii,] <- c(vertex[1],vertex[2],vertex[3])
            if (err_1 > 0){ 
                vertex[2] <- vertex[2] + y_inc
                err_1 <- err_1 - dx2
            }
            if (err_2 > 0){
                vertex[3] <- vertex[3]+ z_inc
                err_2 <- err_2 - dx2
            }
            err_1 <- err_1 + dy2
            err_2 <- err_2 + dz2
            vertex[1] <- vertex[1] + x_inc
            ii <- ii + 1
        }
    }
    else if ((m >= l) && (m >= n)){ 
        err_1 <- dx2 - m
        err_2 <- dz2 - m
        Mat <- matrix(0, ncol=3, nrow=m+1)
        ii <- 1
        for (i in 1:m){ 
            Mat[ii,] <- c(vertex[1],vertex[2],vertex[3])
            if (err_1 > 0){ 
                vertex[1] <- vertex[1] + x_inc
                err_1 <- err_1 - dy2
            }
            if (err_2 > 0){ 
                vertex[3] <- vertex[3] + z_inc
                err_2 <- err_2 - dy2
            }
            err_1 <- err_1 + dx2
            err_2 <- err_2 + dz2
            vertex[2] <- vertex[2] + y_inc
            ii <- ii + 1
        }
    }
    else{ 
        err_1 <- dy2 - n
        err_2 <- dx2 - n
        Mat <- matrix(0, ncol=3, nrow=n+1)
        ii <- 1
        for (i in 1:n){
            Mat[ii,] <- c(vertex[1],vertex[2],vertex[3])
            if (err_1 > 0){ 
                vertex[2] <- vertex[2] + y_inc
                err_1 <- err_1 - dz2
            }
            if (err_2 > 0){ 
                vertex[1] <- vertex[1] + x_inc
                err_2 <- err_2 - dz2
            }
            err_1 <- err_1 + dy2
            err_2 <- err_2 + dx2
            vertex[3] <- vertex[3] + z_inc
            ii <- ii + 1
        }
    }
        
    Mat[ii,] <- c(vertex[1],vertex[2],vertex[3])
    Mat
 }   

linesTetrahedra <- function(x, y, z,
                            delta=c(min(x[,2]-x[,1])/10,
                                    min(y[,2]-y[,1])/10,
                                    min(z[,2]-z[,1])/10),
                            lwd = 0.01, color = "black", ...){

    n <- length(x)
    if (length(y) != n || length(z) != n)
        stop("coordinates must be of the same length.")
    if (is.vector(x)){
        if (!is.vector(y) || !is.vector(z))
            stop("coordinates have to be all vectors or matrices!")
        if (length(x) != 2)
            stop("need to specify the coordinates of starting and ending points.")

        else{
            x <- matrix(x, nrow=1)
            y <- matrix(y, nrow=1)
            z <- matrix(z, nrow=1)
        }
    }
    if (is.matrix(x)){
        if (!is.matrix(y) || !is.matrix(z))
            stop("coordinates have to be all vectors or matrices!")
        if (ncol(x) != 2)
            stop("need to specify the coordinates of starting and ending points.")
    }

    nl <- nrow(x)
    xyz <- do.call(rbind, lapply(1:nl, function(i)
                                 bresenhamLine(x[i,1], y[i,1], z[i,1],
                                               x[i,2], y[i,2], z[i,2],
                                               delta)))
    pointsTetrahedra(xyz[,1], xyz[,2], xyz[,3],
                     size = lwd, color = color, ...) 
}


## Compute for each triangle the indices of triangles that share an
## edge with it.  This could be done more efficiently.
triangleNeighbors <- function(tris) {
   ve <- misc3d:::t2ve(tris)
   vt <- misc3d:::vertexTriangles(ve)
   ib <- ve$ib
   n.tri <- ncol(ib)
   tn <- vector("list", n.tri)
   for (i in 1 : n.tri) {
       v1 <- unique(vt[[ib[1, i]]])
       v2 <- unique(vt[[ib[2, i]]])
       v3 <- unique(vt[[ib[3, i]]])
       i12 <- intersect(v1, v2)
       i23 <- intersect(v2, v3)
       i31 <- intersect(v3, v1)
       u <- union(union(i12, i23), i31)
       tn[[i]] <- u[u != i]
   }
   tn
}
## 'unique' in unique(vt[[ib[1, i]]]) seems to be unnecessary
## unless a triangle has essentially two vertices or one vertex
triangleNeighbors <- function(tris) {
   ve <- misc3d:::t2ve(tris)
   vt <- misc3d:::vertexTriangles(ve)
   ib <- ve$ib
   n.tri <- ncol(ib)
   tn <- vector("list", n.tri)
   for (i in 1 : n.tri) {
       v1 <- vt[[ib[1, i]]]
       v2 <- vt[[ib[2, i]]]
       v3 <- vt[[ib[3, i]]]
       i12 <- intersect(v1, v2)
       i23 <- intersect(v2, v3)
       i31 <- intersect(v3, v1)
       u <- union(union(i12, i23), i31)
       tn[[i]] <- u[u != i]
   }
   tn
}

## Dijkstra's version of Rem's algorithm for computing equivalence
## classes based on a number of vertices 1:nvert and a set of N edges
## provided as an N x 2 matrix.
GetPatches <- function(nvert, edges) {

   f <- 1:nvert

   if (!(is.vector(edges)) && dim(edges)[1] != 0){
       nedge <- nrow(edges)

       for (e in 1:nedge) {
           p0 <- edges[e, 1]
           q0 <- edges[e, 2]
           p1 <- f[p0]
           q1 <- f[q0]
           while (p1 != q1) {
               if (q1 < p1) {
                   f[p0] <- q1
                   p0 <- p1
                   p1 <- f[p1]
               }
               else {
                   f[q0] <- p1
                   q0 <- q1
                   q1 <- f[q1]
               }
           }
       }
   }
   if(is.vector(edges)){
       if(edges[1] < edges[2])
           f[edges[2]] <- edges[1]
       else  f[edges[1]] <- edges[2]
   }

   for (v in 1:nvert)
       f[v] <- f[f[v]]

   split(1:nvert,f)
}

## compute the edges to indicate which triangles share an edge -- this
## needs more error checking
triangleNeighborEdges <- function(tn) {
   edges <- function(i) {
       v <- tn[[i]]
       if (length(v) > 0) cbind(i,v)
       else numeric(0)
   }
   do.call(rbind, lapply(1:length(tn), edges))
}

## separate triangles into disconnected chunks
separateTriangles <- function(contour3dObj){
    tn <- triangleNeighbors(contour3dObj)
    edges <- triangleNeighborEdges(tn)
    edges <- edges[edges[,1] < edges[,2],]
    p <- GetPatches(length(tn), edges)
    newContour3dObj <- vector("list", length(p))
    for(i in 1:length(newContour3dObj)){
        newContour3dObj[[i]] <- contour3dObj
        newContour3dObj[[i]]$v1 <- contour3dObj$v1[p[[i]],]
        newContour3dObj[[i]]$v2 <- contour3dObj$v2[p[[i]],]
        newContour3dObj[[i]]$v3 <- contour3dObj$v3[p[[i]],]
    }
    newContour3dObj

}

