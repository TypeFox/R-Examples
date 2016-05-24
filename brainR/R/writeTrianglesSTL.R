#' Write STL triangles (without recalling the ids)
#' 
#' This is code extracted from \code{\link{writeSTL}} in 
#' \code{rgl}.  This allows users to write the triangles in STL
#' without reprinting the rgl (which takes time)
#'
#' 
#' @param scene list of triangles (that have class Triangles3D)
#' @param con filename or connection of stl file to write
#' @param ascii indicator if the file should be written in
#' ascii or binary
#' @export
#' @return filename (invisible) of stl object


writeTrianglesSTL = function (scene, con, ascii = FALSE) 
{
    ascii = FALSE
    writeHeader <- function() {
        ident <- paste(filename, " produced by RGL\n")
        if (ascii) 
            cat("solid ", ident, file = con)
        else {
            padding <- paste(rep(" ", 80), collapse = "")
            ident <- substr(paste("binary", ident, padding), 
                1, 80)
            writeChar(ident, con, nchars = 80, useBytes = TRUE, 
                eos = NULL)
            writeBin(0L, con, size = 4, endian = "little")
        }
    }
    #### from rgl
    normalize = function(v){
        veclen = sqrt(sum(v^2))
        v/veclen
    }
    xprod = function (v, w) {
        c(v[2] * w[3] - v[3] * w[2], 
            v[3] * w[1] - v[1] * w[3], 
            v[1] * w[2] - v[2] * w[1])
    }
    triangles <- 0
    writeTriangles <- function(vertices) {
        if (nrow(vertices)%%3 != 0) 
            stop("Need 3N vertices")
        n <- nrow(vertices)/3
        for (i in seq_len(n)) {
            vec0 <- vertices[3 * i - 2, ]
            vec1 <- vertices[3 * i - 1, ]
            vec2 <- vertices[3 * i, ]
            normal <- normalize(xprod(vec2 - vec0, vec1 - vec0))
            if (ascii) {
                cat("facet normal ", normal, "\n", file = con)
                cat("outer loop\n", file = con)
                cat("vertex ", vec0, "\n", file = con)
                cat("vertex ", vec1, "\n", file = con)
                cat("vertex ", vec2, "\n", file = con)
                cat("endloop\n", file = con)
                cat("endfacet\n", file = con)
            }
            else {
                writeBin(c(normal, vec0, vec1, vec2), con, size = 4, 
                  endian = "little")
                writeBin(0L, con, size = 2, endian = "little")
            }
        }
        triangles <<- triangles + n
    }

    if (!inherits(scene, "list")){
        stop(paste0("scene must be of class 'list', "
            , "if Triangles3D, use list(scene)"))
    }

    if (is.character(con)) {
        con <- file(con, if (ascii) 
            "w"
        else "wb")
        on.exit(close(con))
    }    
    filename <-  summary(con)$description
    alltri = sapply(scene, inherits, "Triangles3D")
    if (!all(alltri)){
        stop("Only implemented for triangles")
    }

    writeHeader()
    for (i in seq_along(scene)) {
        tri = scene[[i]]
        nvert = nrow(tri$v1) + nrow(tri$v2) + nrow(tri$v3)
        vertices = matrix(NA, nrow=nvert, ncol=3)
        vertices[seq(1, nvert, by=3),] = tri$v1
        vertices[seq(2, nvert, by=3),] = tri$v2
        vertices[seq(3, nvert, by=3),] = tri$v3
        colnames(vertices) = c("x", "y", "z")

        # old.vert = rgl.attrib(ids[i], "vertices")
        # diff = abs(old.vert - vertices)
        # any(diff > 1e-5)
        # [1] FALSE
        # colnames(vertices)
        # dim(old.vert)
        writeTriangles(vertices)
        # cat("written triangles\n")
    }
    if (!ascii) {
        seek(con, 80)
        writeBin(as.integer(triangles), con, size = 4, endian = "little")
    }
    invisible(filename)
}
