Yule04<-
function (X,YY,levX,nameYY,tri=0,alpha=0.05) 
{
    T4 <- Yule03(X, YY, levX,nameYY,tri=tri)
T4<-na.omit(T4)
    couleur<-rep(1,nrow(T4))
    couleur[T4[,4]<alpha]<-2

    L <- nrow(T4)
    MIN <- min(-1, T4[, 2] - 1.96 * T4[, 3])
    MAX <- max(1, T4[, 2] + 1.96 * T4[, 3])
    dotchart(T4[, 2], labels = rownames(T4), xlim = c(MIN, MAX), 
        main = levX,color=couleur,xlab="Yule's Q")
    abline(v = 0, lty = 2)
    segments(T4[, 2] - 1.96 * T4[, 3], 1:L, T4[, 2] + 1.96 * 
        T4[, 3], 1:L,col=couleur)
}
