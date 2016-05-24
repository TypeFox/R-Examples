# Convert Sigma matrix to V matrix
sigma2v <- function(sigma_vec)
{   
    k = sqrt(length(as.vector(sigma_vec))) + 1
    sigma_mat = matrix(sigma_vec, ncol = k - 1)
    A = rbind(cbind(diag(1, ncol = k - 1, nrow = k - 1), rep(-1, k - 1)), rep(1, k))
    B = rbind(cbind(sigma_mat, rep(0, k - 1)), c(rep(0, k - 1), k)) 
    return(solve(A) %*% B %*% solve(t(A)))
}
