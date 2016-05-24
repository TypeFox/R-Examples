# Convert Sigma vector to correlation matrix
sigma2corr <- function(sigma_vec)
{   
    k = sqrt(length(sigma_vec)) + 1
    sigma_mat = matrix(sigma_vec, ncol = k - 1)
    corr_mat = cov2cor(sigma_mat)
    return(as.vector(corr_mat))
}