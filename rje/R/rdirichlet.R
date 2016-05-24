rdirichlet <-
function (n, alpha) 
{
    if (is.matrix(alpha)) {
        k = dim(alpha)[2]
        alpha = as.vector(t(alpha))
    }
    else k = length(alpha)
    out = matrix(rgamma(n * k, shape = alpha), n, k, byrow = TRUE)
    out = out/rowSums(out)
    out
}
