rest.est <-
function(gamma, z, resvector, resmatrix, S)
{
    repeat
    {
        projection <- lsConstrain.fit(x = z, b = resvector, s = S, a = resmatrix,
                                      0 * {1:dim(resmatrix)[1]},
                                      itmax = 4000, eps = 1e-06, eps2 = 1e-06)
        Pz <- projection$x.final
        estimation <- Pz - gamma * {z - Pz}
        if(gamma == 0 || prod(resmatrix %*% estimation <= resvector) == 1)
            break
        else
            z <- estimation
    }
    return(estimation)
}