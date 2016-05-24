permutations <-
function (n, set = 1:n) 
{
    if (!is.vector(n, mode = "numeric") || is.factor(n) || length(n) != 
        1 || n != trunc(n)) 
        stop("'n' must be an integer")
    if (n <= 0) 
        return(NULL)
    set <- as.vector(set)
    if (length(set) != n) 
        stop(paste("'set' must contain", n, "elements"))
    perm.mat <- t(sapply(seq(gamma(n + 1)), perm, n = n))
    if (is.numeric(set) && all.equal(set, 1:n) == TRUE) 
        return(perm.mat)
    else return(matrix(set[perm.mat], factorial(n)))
}
