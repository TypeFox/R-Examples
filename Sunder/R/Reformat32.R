Reformat32 <-
function (g) 
{
    if (length(dim(g)) != 3) 
        stop("g must be 3d")
    nloc = dim(g)[2]
    nind = dim(g)[1]
    xb = matrix(rep(NA, 2 * nloc * nind), nind)
    w1 = which(g == 1, arr.ind = TRUE)
    w2 = which(g == 2, arr.ind = TRUE)
    for (i in 1:nrow(w2)) {
        xb[w2[i, 1], 2 * w2[i, 2]] = w2[i, 3]
        xb[w2[i, 1], 2 * w2[i, 2] - 1] = w2[i, 3]
    }
    for (i in 1:nrow(w1)) {
        xx = w1[i, 1]
        yy = 2 * w1[i, 2] - 1
        zz = w1[i, 3]
        if (!is.na(xb[xx, yy])) 
            yy = yy + 1
        xb[xx, yy] = zz
    }
    y = tapply(as.vector(xb), rep(1:nloc, each = 2 * nind), function(x) as.numeric(as.factor(x)))
    nal = as.vector(sapply(y, max, na.rm = TRUE))
    return(xb)
}
