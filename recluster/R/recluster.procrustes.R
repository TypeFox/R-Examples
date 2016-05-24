recluster.procrustes<- function (X, Y, Yv=FALSE, num=nrow(X), scale = TRUE, ...) 
{
    Yvrot<-NULL    
    X <- scores(X, display = "sites")
    Y <- scores(Y, display = "sites")
    if (abs(sum(Yv))>0) {
		Yv<-rbind(c(0,0),Yv)
    	Yv <- scores(Yv, display = "sites")
	}
    ctrace <- function(MAT) sum(MAT^2)
    c <- 1
    xmean <- apply(X[1:num,], 2, mean)
    ymean <- apply(Y[1:num,], 2, mean)
    X[,1] <- X[,1]-xmean[1]
    X[,2] <- X[,2]-xmean[2]
    Y[,1] <- Y[,1]-ymean[1]
    Y[,2] <- Y[,2]-ymean[2]
    if (abs(sum(Yv))>0) {
		Yv[,1] <- Yv[,1]-ymean[1]
		Yv[,2] <- Yv[,2]-ymean[2]
	}
    XY <- crossprod(X[1:num,], Y[1:num,])
    sol <- svd(XY)
    A <- sol$v %*% t(sol$u)
    if (scale) {
        c <- sum(sol$d)/ctrace(Y[1:num,])
    }
    Yrot <- c * Y %*% A
    if (abs(sum(Yv))>0)	{
		Yvrot <- c * Yv %*% A
	}
    b <- xmean - c * ymean %*% A
    R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
    reslt <- list(Yrot = Yrot, X = X, Yvrot=Yvrot, ss = R2, rotation = A, translation = b, scale = c, xmean = xmean, call = match.call())
    reslt$svd <- sol
    class(reslt) <- "procrustes"
    reslt
}
