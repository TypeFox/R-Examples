robout <-
function (data, nclass=0, meth = c("mve", "mcd"), rep = 10, plot = TRUE) 
{
#    require(MASS)
    if (sum(is.na(data))> 0) 
        stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    p=dim(data)[2]
    if(nclass==0)
        data = as.data.frame(data[,-p])
    else        data=data[as.numeric(factor(data[,p]))==nclass,-p ]
    #ncol = dim(data)[2]
    #tempo = data[data[, ncol] == nclass, 1:(ncol - 1)]
    tempo=data
    namestempo = rownames(tempo)
    nrow = dim(tempo)[1]
    roboutl = NULL
    roboutall = matrix(0, nrow, rep)
    rownames(roboutall) = namestempo
    for (i in 1:rep) {
        mcdc = MASS::cov.rob(tempo, method = meth)
        mbc = sqrt(mahalanobis(tempo, mcdc$center, mcdc$cov, 
            to = 1e-14))
        roboutl = c(roboutl, boxplot(mbc, plot = FALSE)$out)
        roboutall[, i] = mbc
    }
    a = as.matrix(roboutl)
    b = apply(roboutall, 1, mean)
    outme = rev(sort(b))
    topo = rev(sort(b))[1:10]
    if (plot) {
        plot(rev(sort(b)), ylab = paste("Mahalabobis distance(", 
            meth, ")"))
        text(1:10, topo, names(topo), cex = 0.6, pos = 4)
    }
    top = table(as.numeric(rownames(a)))
    top1 = top[top > ceiling(rep/2)]
    cat("\nTop outliers by frequency\n")
    print(top1)
    topout = as.numeric(names(top1))
    ntops = length(topout)
    outly = rep(0, ntops)
    for (i in 1:ntops) {
        outly[i] = mean(a[as.numeric(rownames(a)) == topout[i]])
    }
    cat("\nTop outliers by outlyngness measure\n")
    print(outme[1:ntops])
    cat("\n")
    list(outme = outme)
}
