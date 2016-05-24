fcondXA2 <-
function (x, distr, Tauy, Tauz, J) 
{
    pJ <- J/sum(J)
    K <- matrix(NA, nrow = length(Tauy), ncol = length(x))
    for (i in seq(Tauy)) {
        K[i, ] <- dk(x, distr = distr, mu = Tauy[i], sigma = Tauz[i])
    }
    fcondXA2 <- apply(K, 2, function(x) sum(x * pJ))
    return(fcondXA2)
}
