derivmclm <-
function (mueta, ncategoriesm1, X_mat) 
{
    nobs <- nrow(X_mat)/ncategoriesm1
    ans <- diagmod(rep.int(1, ncategoriesm1))
    ans[seq(2, ncategoriesm1^2, ncategoriesm1 + 1)] <- -1
    ans <- apply(ans, 2, function(x) rep.int(x, nobs))
    mat1 <- matrix(mueta, nobs, ncategoriesm1, TRUE)
    mat1 <- apply(mat1, 2, function(x) rep(x, each = ncategoriesm1))
    mat1 <- ans * mat1
    mat2 <- .rowSums(mat1, nrow(mat1), ncol(mat1),FALSE) * X_mat[, 
        -(1:ncategoriesm1)]
    mat2 <- cbind(mat1, mat2)
    mat2
}

