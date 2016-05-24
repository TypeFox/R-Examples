## Lighting functions are functions of the form
##
##     lighting(normals, view, light, color, color2, material)
##
## with
##
##    normals   a matrix of unit normal vectors as computed by triangleNormals
##    view      vector pointing to the viewer
##    light     vector pointing to the light source
##    color     color or vector of n colors for the sides of the triangles
##              in the direction of to the normal vectors
##    color2    color or vector of n colors for the sides of the triangles
##              in the opposite direction of to the normal vectors. NA means
##              same as 'color'
##    alpha     the alpha level to use
##    material  material parameters controlling the lighting calculation.
##
## lighting functions return a vector of n rgb colors corresponding
## to the sides of the trinagles facing the viewer and the lighting
## algorithm.

## phongLighting implements a simple version of the Phong lighting
## model (not shading--that would involve interpolation within the
## triangles). It incorporates ambient and diffuse light, which are
## the same color as the object, and specular light, which is a convex
## combination of the object color and the (white) light color.  This
## is based roughly on the description in Foley and Van Dam.

phongLighting <- function(normals, view, light, color, color2, alpha,
                          material = "default") {
    if (length(light) == 4) {
        LI <- light[4]
        light <- light[1:3]
    }
    else LI <- 1
    if (is.character(material))
        material <- getMaterial(material)
    ambient <- material$ambient
    diffuse <- material$diffuse
    specular <- material$specular
    exponent <- material$exponent
    sr <- material$sr
    V <- view / sqrt(sum(view^2))
    L <- light / sqrt(sum(light^2))
    H <- (L + V) / sqrt(sum((L + V)^2))
    sgn <- as.vector(normals %*% V) > 0
    N <- ifelse(sgn,1, -1) * normals
    Is <- as.vector(specular * abs(N %*% H) ^ exponent)
    Id <-  as.vector(diffuse * pmax(N %*% L,0))
    rgbcol <- t(col2rgb(ifelse(sgn, color, color2)) / 255)
    Lrgbcol <- pmin(LI * ((ambient + Id + sr * Is) * rgbcol + (1 - sr) * Is),
                    1)
    Lrgbcol[is.na(Lrgbcol)] <- 0
    rgb(Lrgbcol[,1], Lrgbcol[,2], Lrgbcol[,3], alpha)
}

## A simple data base is used to register properties of named
## materials.  Some initial entries are based on valuaes for similarly
## named materials in Matlab.

materials.database <- new.env(hash = TRUE)

registerMaterial <- function(name, ambient = 0.3, diffuse = 0.7,
                             specular = 0.1, exponent = 10, sr = 0) {
    value <- list(ambient = ambient, diffuse = diffuse,
                  specular = specular, exponent = exponent, sr = sr)
    assign(name, value, materials.database)
}

getMaterial <- function(name) {
    if (exists(name, materials.database, inherits = FALSE))
        get(name, materials.database)
    else get("default", materials.database, inherits = FALSE)
}

registerMaterial("shiny", ambient = 0.3, diffuse = 0.6, specular = 0.9,
                 exponent = 20, sr = 0)
registerMaterial("dull", ambient = 0.3, diffuse = 0.8, specular = 0.0,
                 exponent = 10, sr = 0)
registerMaterial("metal", ambient = 0.3, diffuse = 0.3, specular = 1.0,
                  exponent = 25, sr = 0.5)
registerMaterial("default", ambient = 0.3, diffuse = 0.7, specular = 0.1,
                 exponent = 10, sr = 0)

# Alternate version of metal, about 50% brighter?
registerMaterial("metal", ambient = 0.45, diffuse = 0.45, specular = 1.5,
                  exponent = 25, sr = 0.5)
# 50% would be 0.45 0.45 1.50?

# Alternate version of shiny, about 20% brighter?
registerMaterial("shiny", ambient = 0.36, diffuse = 0.72, specular = 1.08,
                 exponent = 20, sr = 0)


## perspLighting is an implementation of the lighting algorithm
## described in the help page for persp().  The 'shade' parameter of
## persp is here named is computed from the material's 'exponent'
## component.  To make the "default" material with expone t = 10
## correspond to the shade = 0.75 value of the volcano example from
## the persp help page, 'exponent' is scaled by a factor of 3 / 40.

perspLighting <- function(normals, view, light, color, color2, alpha,
                          material = "default") {
    if (length(light) == 4) {
        LI <- light[4]
        light <- light[1:3]
    }
    else LI <- 1
    if (is.character(material))
        material <- getMaterial(material)
    exponent <- material$exponent
    V <- view / sqrt(sum(view^2))
    L <- light / sqrt(sum(light^2))
    sgn <- as.vector(normals %*% V) > 0
    N <- ifelse(sgn,1, -1) * normals
    I <-  (pmax(1 + as.vector(N %*% L), 0) / 2) ^ (exponent * (3 / 40))
    Lrgbcol <- I * LI * t(col2rgb(ifelse(sgn, color, color2)) / 255)
    rgb(Lrgbcol[,1], Lrgbcol[,2], Lrgbcol[,3], alpha)
}

triangleNormalsPhong <- function(triangles) {
    N <- triangleNormals(triangles)
    ve <- t2ve(triangles)
    vt <- vertexTriangles(ve)
    VN <- vertexNormals(vt, N)
    interpolateVertexNormals(VN, ve$ib)
}

triangleNormalsPhongEX <- function(triangles, reps = 1) {
    N <- triangleNormals(triangles)
    ve <- t2ve(triangles)
    vt <- vertexTriangles(ve)
    VN <- vertexNormals(vt, N)
    vb <- ve$vb
    ib <- ve$ib
    n.tri <- nrow(N)

    while (reps > 0) {
	reps <- reps - 1
	n.ver <- nrow(VN)
	mt <- triangleMidTriangles(vb, ib, VN)
	vb <- cbind(vb, mt$vb)
	VN <- rbind(VN, mt$VN)
	mtib <- mt$ib + n.ver
	ib <- matrix(rbind(ib[1,], mtib[1,],mtib[3,],
			   mtib[1,], ib[2,],mtib[2,],
			   mtib[2,], ib[3,], mtib[3,],
			   mtib),
		     nrow = 3)

	for (i in seq(along = triangles))
	    if (length(triangles[[i]]) == n.tri)
		triangles[[i]] <- rep(triangles[[i]], each = 4)

	n.tri <- 4 * n.tri
    }

    triangles$N <- interpolateVertexNormals(VN, ib)
    triangles$v1 <- t(vb[,ib[1,]])
    triangles$v2 <- t(vb[,ib[2,]])
    triangles$v3 <- t(vb[,ib[3,]])

    triangles
}

# version that handles color interpolation
# **** could lift out and triangleEdges calls
triangleNormalsPhongEX <- function(triangles, reps = 1) {
    N <- triangleNormals(triangles)
    ve <- t2ve(triangles)
    vt <- vertexTriangles(ve)
    VN <- vertexNormals(vt, N)
    vb <- ve$vb
    ib <- ve$ib
    n.tri <- nrow(N)

    color <- rep(triangles$color, length = n.tri)
    color2 <- rep(triangles$color2, length = n.tri)
    col.mesh <- rep(triangles$col.mesh, length = n.tri)
    color2 <- ifelse(is.na(color2), color, color2)
    col.mesh <- ifelse(is.na(col.mesh), color, col.mesh)

    VC <- vertexColors(vt, color)
    VC2 <- vertexColors(vt, color2)
    VCm <- vertexColors(vt, col.mesh)

    while (reps > 0) {
        reps <- reps - 1
        n.ver <- nrow(VN)

        edges <- triangleEdges(vb, ib)
        VC <- rbind(VC, (VC[edges[1,],] + VC[edges[2,],]) / 2)
        VC2 <- rbind(VC2, (VC2[edges[1,],] + VC2[edges[2,],]) / 2)
        VCm <- rbind(VCm, (VCm[edges[1,],] + VCm[edges[2,],]) / 2)

        mt <- triangleMidTriangles(vb, ib, VN)
        vb <- cbind(vb, mt$vb)
        VN <- rbind(VN, mt$VN)
        mtib <- mt$ib + n.ver
        ib <- matrix(rbind(ib[1,], mtib[1,],mtib[3,],
                           mtib[1,], ib[2,],mtib[2,],
                           mtib[2,], ib[3,], mtib[3,],
                           mtib),
                     nrow = 3)

        for (i in seq(along = triangles))
            if (length(triangles[[i]]) == n.tri)
                triangles[[i]] <- rep(triangles[[i]], each = 4)

        n.tri <- 4 * n.tri
    }

    triangles$color <- interpolateVertexColors(VC, ib)
    triangles$color2 <- interpolateVertexColors(VC2, ib)
    triangles$color.mesh <- interpolateVertexColors(VCm, ib)

    triangles$N <- interpolateVertexNormals(VN, ib)
    triangles$v1 <- t(vb[,ib[1,]])
    triangles$v2 <- t(vb[,ib[2,]])
    triangles$v3 <- t(vb[,ib[3,]])

    triangles
}

lightTriangles <- function(triangles, lighting, light) {
    view <- c(0, 0, 1)
    normals <- triangleNormals(triangles)
    smooth <- if (is.null(triangles$smooth)) 0 else triangles$smooth
    if (smooth == 0)
        normals <- triangleNormals(triangles)
    else if (smooth == 1)
        normals <- triangleNormalsPhong(triangles)
    else {
        triangles <- triangleNormalsPhongEX(triangles, reps = smooth - 1)
        normals <- triangles$N
    }
    n.tri <- nrow(normals)
    color <- rep(triangles$color, length = n.tri)
    color2 <- rep(triangles$color2, length = n.tri)
    color2 <- ifelse(is.na(color2), color, color2)
    alpha <- rep(triangles$alpha, length = n.tri)
    mat <- triangles$material
    triangles$col.light <- lighting(normals, view, light, color, color2, 
                                    alpha, mat)
    triangles
}

lightScene <- function(scene, lighting, light) {
    if (is.Triangles3D(scene))
        lightTriangles(scene, lighting, light)
    else lapply(scene, lightTriangles, lighting, light)
}

