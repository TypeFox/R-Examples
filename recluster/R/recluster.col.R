recluster.col<-function (mat, st = TRUE, rot = TRUE) 
{
    mat2 <- mat
    if (rot) {
        if ((max(mat[, 1]) - min(mat[, 1])) < (max(mat[, 2]) - 
            min(mat[, 2]))) {
            mat2[, 1] <- mat[, 2]
            mat2[, 2] <- mat[, 1]
        }
    }
    mat3 <- mat2
    if ((max(mat2[, 1]) - min(mat2[, 1])) < (max(mat2[, 2]) - 
        min(mat2[, 2]))) {
        mat3[, 2] <- (mat2[, 2] - min(mat2[, 2]))/(max(mat2[, 
            2]) - min(mat2[, 2]))
        mat3[, 1] <- (mat2[, 1] - min(mat2[, 1]))/(max(mat2[, 
            2]) - min(mat2[, 2]))
        border <- (1 - max(mat3[, 1]))/2
        mat3[, 1] <- mat3[, 1] + border
    }
    if ((max(mat2[, 1]) - min(mat2[, 1])) >= (max(mat2[, 2]) - 
        min(mat2[, 2]))) {
        mat3[, 1] <- (mat2[, 1] - min(mat2[, 1]))/(max(mat2[, 
            1]) - min(mat2[, 1]))
        mat3[, 2] <- (mat2[, 2] - min(mat2[, 2]))/(max(mat2[, 
            1]) - min(mat2[, 1]))
        border <- (1 - max(mat3[, 2]))/2
        mat3[, 2] <- mat3[, 2] + border
    }
    colour <- array(data = 0, dim = c(dim(mat3)[1], (dim(mat3)[2]) + 
        3))
    for (t in 1:dim(mat3)[1]) {
        colour[t, 1] <- mat3[t, 1]
        colour[t, 2] <- mat3[t, 2]
        colour[t, 3] <- round(mat3[t, 2] * 255)
        colour[t, 4] <- round(mat3[t, 1] * 255)
        colour[t, 5] <- round((1 - max(mat3[t, 2],mat3[t, 1])) * 255)
        rownames(colour) <- rownames(mat)
        if (!st) {
            colour[, 1] <- mat2[, 1]
            colour[, 2] <- mat2[, 2]
        }
    }
    return(colour)
}
