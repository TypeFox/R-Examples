AllHarmonicMean <- function (vectors) 
{
    n = length(vectors)/2
    m = rep(0, n)
    for (i in 1:n) {
        m[i] = sum_distances2 = n / (sum(1 / sqrt((vectors[i, 2] - 
            vectors[i, 1]) * (vectors[i, 2] - vectors[i, 1]) + 
            (vectors[-i, 2] - vectors[-i, 1]) * (vectors[-i, 
                2] - vectors[-i, 1]))))
    }
    return(m)
}
