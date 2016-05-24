mahaout <-
function (data, nclass=0, plot = TRUE) 
{
#    require(MASS)
    if (sum(is.na(data))> 0) 
        stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    p=dim(data)[2]
    if(nclass==0)
        data = as.data.frame(data[,-p])
    else        data=data[as.numeric(factor(data[,p]))==nclass,-p ]
#    ncol = dim(data)[2]
    tempo = data
    namestempo = rownames(tempo)
    nrow = dim(tempo)[1]
    roboutl = NULL
    mcdc = MASS::cov.rob(tempo, method = "classical")
    mbc = sqrt(mahalanobis(tempo, mcdc$center, mcdc$cov, to = 1e-14))
    roboutl = c(roboutl, boxplot(mbc, plot = FALSE)$out)
    cat("Ouliers given by the boxplot of the  Mahalanobis distance\n")
    print(rev(sort(roboutl)))
    outme = rev(sort(mbc))
    topo = rev(sort(mbc))[1:10]
    if (plot) {
        plot(rev(sort(mbc)), ylab = "Mahalabobis distance")
        text(1:10, topo, names(topo), cex = 0.6, pos = 4)
    }
    cat("\n")
    list(outme = outme)
}
