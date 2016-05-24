robout <- function (data, nclass, meth = c("mve", "mcd"), rep = 10) 
{
    ncol = dim(data)[2]
    tempo = data[data[, ncol] == nclass, 1:(ncol - 1)]
    namestempo = rownames(tempo)
    nrow = dim(tempo)[1]
    roboutl = NULL
    roboutall = matrix(0, nrow, rep)
    rownames(roboutall) = namestempo
    for (i in 1:rep) {
        mcdc = cov.rob(tempo, method = meth)
        mbc = sqrt(mahalanobis(tempo, mcdc$center, mcdc$cov, 
                   to = 1e-14))
        roboutl = c(roboutl, boxplot(mbc, plot = FALSE)$out)
        roboutall[, i] = mbc
    }
    a = as.matrix(roboutl)
    b = apply(roboutall, 1, mean)
    outme = rev(sort(b))
    topo = rev(sort(b))[1:10]
    top = table(as.numeric(rownames(a)))
    top1 = top[top > ceiling(rep / 2)]
    topout = as.numeric(names(top1))
    ntops = length(topout)
    outly = rep(0, ntops)
    for (i in 1:ntops) {
        outly[i] = mean(a[as.numeric(rownames(a)) == topout[i]])
    }
    return(list(outliers = topout, depth.total = outme))
}
