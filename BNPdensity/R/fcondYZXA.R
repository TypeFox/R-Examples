fcondYZXA <-
function (x, distr, Tauy, Tauz, J) 
{
    K <- matrix(NA, nrow = length(Tauy), ncol = length(x))
    for (i in seq(Tauy)) {
        K[i, ] <- dk(x, distr = distr, mu = Tauy[i], sigma = Tauz[i]) * 
            J[i]
    }
    if (any(is.na(K))) 
        print(K, Tauy, Tauz, J)
    pK <- prop.table(K, margin = 2)
    j <- apply(pK, 2, function(x) sample(length(Tauy), size = 1, 
        prob = x))
    return(matrix(c(y = Tauy[j], z = Tauz[j]), nrow = length(x), 
        ncol = 2))
}
