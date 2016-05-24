derivbcl <-
function (fitprob, ncategoriesm1, X_mat) 
{
    n1 <- nrow(X_mat)
    n2 <- ncol(X_mat)
    fitprob1 <- matrix(fitprob, n1/ncategoriesm1, ncategoriesm1, 
        TRUE)
    mat1 <- matrix(0, n1, ncategoriesm1)
    for (i in 1:nrow(fitprob1)) {
        mat1[((i - 1) * ncategoriesm1 + 1):(i * ncategoriesm1), 
            ] <- diagmod(fitprob1[i, ]) - tcrossprod(fitprob1[i, 
            ])
    }
    if (n2 == ncategoriesm1) 
        ans <- mat1
    else {
        mat1 <- t(apply(mat1, 1, function(x) rep(x, each = n2/ncategoriesm1)))
        mat2 <- X_mat[seq(1, n1, by = ncategoriesm1), 1:(n2/ncategoriesm1)]
        mat2 <- t(apply(mat2, 1, function(x) rep(x, ncategoriesm1)))
        mat2 <- apply(mat2, 2, function(x) rep(x, each = ncategoriesm1))
        ans <- mat1 * mat2
    }
    ans
}

