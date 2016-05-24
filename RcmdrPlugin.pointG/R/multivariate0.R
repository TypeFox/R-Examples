multivariate0 <-function (XX, YY = NULL) 
{
    X <- XX
    n <- nrow(X)
    X1 <- data.frame(1:n, X)
    X2 <- na.omit(X1)
    X.NA <- X2[, -1]
    selection <- X2[, 1]
    dudi.X <- dudi.mix(X.NA,scannf=FALSE)
dev.new()
barplot(dudi.X$eig,col=rev(brewer.pal(3,name="PuRd"))[1])
    plotNum(X.NA, dudi.X$li)
    plotCat2(X.NA, dudi.X$li)
    if (!is.null(YY)) {
        Y <- YY
        Y.NA <- Y[selection, ]
        plotNum(Y.NA, dudi.X$li)
        plotCat2(Y.NA, dudi.X$li)
    }
    dudi.X
}

