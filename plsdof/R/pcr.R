pcr<-function (X, y, scale = TRUE, m = min(ncol(X), nrow(X) - 1), 
    eps = 1e-06,supervised=FALSE) 
{
    p <- ncol(X)
    n <- nrow(X)
    Beta <- matrix(, p, m)
    mean.y <- mean(y)
    y <- scale(y, scale = FALSE)
    mean.X <- apply(X, 2, mean)
    if (scale == FALSE) {
        sd.X <- rep(1, p)
    }
    if (scale == TRUE) {
        sd.X <- apply(X, 2, sd)
        sd.X[sd.X == 0] = 1
    }
    X <- X - rep(1, nrow(X)) %*% t(mean.X)
    X <- X/(rep(1, nrow(X)) %*% t(sd.X))
    my.svd = svd(X)
    sigma<- (my.svd$d)[1:m]
    U <- my.svd$v[,1:m]
    V<-my.svd$u[,1:m]
    
    if (supervised==TRUE){
	cor2<-as.vector(cor(V,y))^2
	my.sort=order(cor2,decreasing=TRUE)
	sigma=sigma[my.sort]
	U=U[,my.sort]
	V=V[,my.sort]
    }
    
    sigmainv <- 1/sigma
    sigmainv[sigma < eps] = 0
    my.d<-as.vector(t(V) %*% y)
    #cat(paste("length of my.d: ",length(my.d),"\n"))
    #cat(paste("length of sinv: ",length(sigmainv),"\n"))
    dummy <- sigmainv*my.d
    for (i in 1:m) {
        Beta[, i] <- U[, 1:i, drop = FALSE] %*% dummy[1:i]
    }
    coefficients <- matrix(0, p, m + 1)
    coefficients[, 2:(m + 1)] = Beta/(sd.X %*% t(rep(1, m)))
    intercept <- rep(mean.y, m + 1) - t(coefficients) %*% mean.X
    return(list(intercept = intercept, coefficients = coefficients))
}
