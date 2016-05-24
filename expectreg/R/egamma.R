egamma <-
function (asy, shape, rate = 1, scale = 1/rate) 
{
    asy[asy > 1 | asy < 0] = NA
    zz = 0 * asy
    lower = rep(0, length(asy))
    upper = rep(1.5 * qgamma(0.995, shape, rate), length(asy))
    diff = 1
    index = 1
    while (diff > 1e-10 && index < 1000) {
        root = pegamma(zz, shape, rate) - asy
        root[is.na(root)] = 0
        lower[root < 0] = zz[root < 0]
        upper[root > 0] = zz[root > 0]
        zz = (upper + lower)/2
        diff = max(abs(root), na.rm = T)
        index = index + 1
    }
    zz[is.na(asy)] = NA
    return(zz)
}
