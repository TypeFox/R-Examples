# Convert the triangula vector to symetric matrix
triangula2mat <- function(tri)
{
    k = (sqrt(1 + 8 * length(tri)) - 1) / 2
    tempmat = matrix(rep(0, k^2), ncol = k)
    lower_ind = lower.tri(tempmat, diag = TRUE)
    tempmat[lower_ind] = as.numeric(tri)
    upper_ind = upper.tri(tempmat, diag = FALSE)
    lower_ind = lower.tri(tempmat, diag = FALSE)
    tempmat[upper_ind] = tempmat[lower_ind]
    return(tempmat)
}