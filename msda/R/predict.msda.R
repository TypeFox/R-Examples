predict.msda <- function(object, newx, ...) {
    theta <- object$theta
    mu <- object$mu
    prior <- object$prior
    mubar <- sweep(mu[, -1], 1, mu[, 1], "+")/2
    n <- nrow(newx)
    p <- ncol(newx)
    x.train <- object$x
    y.train <- object$y
    nclass <- length(prior)
    nlambda <- length(theta)
    pred <- matrix(0, n, nlambda)
    pred[1] <- which.max(prior)
    for (i in 1:nlambda) {
        nz <- sum(theta[[i]][, 1] != 0)
        if (nz == 0) {
            pred[,i] <- which.max(prior)
        } else {
            xfit <- x.train %*% theta[[i]][, 1:(min(nclass - 1, nz)),drop=FALSE]
            xfit.sd<-matrix(0,nclass,ncol(xfit))
            for(j in 1:nclass){
            xfit.sd[j,]<-apply(xfit[y.train==j,,drop=FALSE],2,sd)}
            xfit.sd<-apply(xfit.sd,2,min)
            if(min(xfit.sd)<1e-4){pred[,i]<-which.max(prior)}else{
            l <- lda(xfit, y.train)
            pred[, i] <- predict(l, newx %*% theta[[i]][, 1:(min(nclass - 
                1, nz))])$class}
        }
    }
    pred
}
