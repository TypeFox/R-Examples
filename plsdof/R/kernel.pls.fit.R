kernel.pls.fit=function (X, y, m=ncol(X), compute.jacobian=FALSE,DoF.max=min(ncol(X)+1,nrow(X)-1)) {
    p <- ncol(X)
    n <- nrow(X)
    m<-min(m,DoF.max)
    #Beta <- matrix(, p, m) # matrix of regression coefficients
    #W <- V <- Beta
    #dW <- dBeta <- dV <- array(dim = c(m, p, n))
    X0<-X
    y0<-y
    # scaling of the data
    mean.y <- mean(y)
    y <- scale(y, scale = FALSE)
    mean.X <- apply(X, 2, mean)
    sd.X <- apply(X, 2, sd)
    sd.X[sd.X == 0] = 1
    X <- X - rep(1, nrow(X)) %*% t(mean.X)
    X <- X/(rep(1, nrow(X)) %*% t(sd.X))
    K<-X%*%t(X)
    #dcoefficients=NULL
    yhat<-vector(length=m+1)
    TT <- matrix(, n, m)
    Yhat <- matrix(, n, m)
    Alpha <- matrix(0, n, m)
    Gamma <- matrix(, n, m)
    DoF = NULL
    dYhat = NULL
    if (compute.jacobian == TRUE) {
        dYhat = array(dim = c(m, n, n))
        dtildeT <- array(dim = c(m, n, n))
        dT <- dtildeT
        DoF = vector(length = m)
    }
    for (i in 1:m) {
        if (i == 1) {
            ri <- y
            gi <- ri
            ti <- tildeti <- K %*% y
            if (compute.jacobian == TRUE) {
                dT[i, , ] <- dtildeT[i, , ] <- K
                dT[i, , ] <- dnormalize(ti, dT[i, , ])
                dtildeT[i, , ] <- dT[i, , ]
            }
            dummy <- normalize(ti, gi)
            ti <- dummy$v
            gi <- dummy$w
            tildeti <- ti
            TT[, i] <- ti
            Gamma[, i] <- gi
            Alpha[, i] <- Gamma[, i] * sum(ti * y)
            Yhat[, i] = vvtz(ti, y)
            if (compute.jacobian == TRUE) {
                ti <- as.vector(ti)
                dYhat[i, , ] = dvvtz(ti, y, dT[i, , ], diag(n))
                DoF[i] <- sum(diag(dYhat[i, , ]))
            }
        }
        if (i > 1) {
            ri <- y - Yhat[, i - 1]
            tildeti <- K %*% ri
            if (compute.jacobian == TRUE) {
                dtildeT[i, , ] <- K %*% (diag(n) - dYhat[i - 
                  1, , ])
            }
            TTi <- TT[, 1:(i - 1), drop = FALSE]
            ti <- tildeti - vvtz(TTi, tildeti)
            dummy <- rep(0, n)
            for (r in 1:(i - 1)) {
                dummy <- dummy + sum(tildeti * TT[, r]) * Gamma[,r]
            }
            gi <- ri - dummy
            if (compute.jacobian == TRUE) {
                dT[i, , ] = dtildeT[i, , ] - dvvtz(TTi, tildeti, 
                  dT[1:(i - 1), , , drop = FALSE], dtildeT[i, 
                    , ])
                dT[i, , ] = dnormalize(ti, dT[i, , ])
            }
            dummy <- normalize(ti, gi)
            ti <- dummy$v
            gi <- dummy$w
            TT[, i] <- ti
            Gamma[, i] <- gi
            Yhat[, i] <- Yhat[, i - 1] + vvtz(ti, y)
            #yhat[i]<-sum(Yhat[,i]^2)
            Alpha[, i] <- Alpha[, i - 1] + Gamma[, i] * sum(TT[, 
                i] * y)
            if (compute.jacobian == TRUE) {
                ti <- as.vector(ti)
                dYhat[i, , ] <- dYhat[i - 1, , ] + dvvtz(ti, 
                  y, dT[i, , ], diag(n))
                DoF[i] <- sum(diag(dYhat[i, , ]))
            }
        }
    }
    Yhat<-cbind(rep(0,n),Yhat)
    Yhat<-Yhat+mean.y
    yhat<-apply(Yhat^2,2,sum)
    Alpha<-cbind(rep(0,n),Alpha)
    if (compute.jacobian==TRUE){
        DoF<-c(0,DoF)+1
    }
    if (compute.jacobian==FALSE){
        DoF<-1:(m+1)
    }
        RSS<-vector(length=m+1)
        for (i in 1:(m+1)) {
            res <- Yhat[, i] - y0
            RSS[i]<-sum(res^2)
        }
        if ((compute.jacobian) == TRUE) {
        dummy<-array(1/n,dim=c(m+1,n,n))
        dummy[2:(m+1),,]<-dummy[2:(m+1),,] + dYhat
        dYhat<-dummy
        denominator <- vector(length = m+1)
        for (i in 1:(m+1)) {
            denominator[i] <- sum(diag((diag(n) - dYhat[i, , 
                ]) %*% ((diag(n) - t(dYhat[i, , ])))))
        }
        sigmahat <- sqrt(RSS/denominator)
    }
    if (compute.jacobian==FALSE){
        sigmahat<-sqrt(RSS/(n-DoF))
    }
    coefficients<-t(X)%*%Alpha
    coefficients = coefficients/(sd.X %*% t(rep(1, m+1)))
    intercept <- rep(mean.y, m+1) - t(coefficients) %*% mean.X
    intercept<-as.vector(intercept)
    DoF[DoF>DoF.max]=DoF.max
    covariance=NULL
    return(list(coefficients=coefficients, intercept=intercept,DoF = DoF, Yhat=Yhat,yhat=yhat,
        RSS = RSS, sigmahat = sigmahat,covariance=covariance, TT = TT))
}
