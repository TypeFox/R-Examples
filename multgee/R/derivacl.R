derivacl <-
function (fitprob, ncategoriesm1, X_mat) 
{
    n1 <- nrow(X_mat)
    n2 <- ncol(X_mat)
    z1 <- n1/ncategoriesm1
    fitprob1 <- matrix(fitprob, z1, ncategoriesm1,TRUE)
    cumprob1 <- t(apply(fitprob1, 1, cumsum))
    mat1 <- matrix(fitprob, n1, n2)
    mat4 <- matrix(0,n1,n2, TRUE)
    for (i in 1:nrow(fitprob1)) {
        mat2 <- matrix(1 - cumprob1[i, ], ncategoriesm1, ncategoriesm1, TRUE)
        mat3 <- 1 - mat2
        mat2[lower.tri(mat2)] <- -mat3[lower.tri(mat3)]
        mat4[((i - 1) * ncategoriesm1 + 1):(i * ncategoriesm1), 
            1:ncategoriesm1] <- mat2
    }
    if (n2 > ncategoriesm1) {
        dummy2 <- (ncategoriesm1 + 1):n2
        if (length(dummy2) == 1) {
            mat2 <- X_mat[, dummy2] * mat1[, dummy2]
            for (i in 1:nrow(fitprob1)) {
                dummy <- ((i - 1) * ncategoriesm1 + 1):(i * ncategoriesm1)
                mat4[dummy, dummy2] <- X_mat[dummy, dummy2] - 
                  rep(sum(mat2[dummy]), ncategoriesm1)
            }
        }
        else {
            dummy3 <- n2 - ncategoriesm1
            mat2 <- X_mat[, dummy2] * mat1[, dummy2]
            for (i in 1:nrow(fitprob1)) {
                dummy <- ((i - 1) * ncategoriesm1 + 1):(i * ncategoriesm1)
                mat3 <- .colSums(mat2[dummy, ], ncategoriesm1,dummy3,FALSE)
                mat3 <- X_mat[dummy, dummy2] - matrix(mat3, 
                  ncategoriesm1, dummy3, TRUE)
                mat4[dummy, (ncategoriesm1 + 1):n2] <- mat3
            }
        }
    }
    mat4 * mat1
}

