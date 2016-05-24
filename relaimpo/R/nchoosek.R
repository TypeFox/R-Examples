"nchoosek" <-
function (n, k) 
{
# function for calculating all subsets of size k from n objects
# taken from package vsn (provided by Wolfgang Huber under LGPL)
# slightly modified to also work for n=k=2 by Ulrike Groemping
    if (!is.numeric(n) || !is.numeric(k) || is.na(n) || is.na(k) || 
        length(n) != 1 || length(k) != 1) 
        stop("arguments must be non-NA numeric scalars.")
    if (k > n || k < 0) 
        stop("Arguments must satisfy 0 <= k <= n.")
    nck = choose(n, k)
    res = matrix(NA, nrow = k, ncol = nck)
    res[, 1] = 1:k
    j = 2
    repeat {
        if (j > nck) 
            break
        res[, j] = res[, j - 1]
        i = k
        repeat {
            res[i, j] = res[i, j] + 1
            if (res[i, j] <= n - (k - i)) 
                break
            i = i - 1
            stopifnot(i >= 1)
        }
        if (i < k) 
            res[(i + 1):k, j] = res[i, j] + 1:(k - i)
        j = j + 1
    }
    stopifnot(all(res[, nck] == (n - k + 1):n))
    stopifnot(all(res <= n) && all(res >= 1))
    return(res)
}

