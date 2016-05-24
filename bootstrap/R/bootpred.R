"bootpred" <- function(x,y,nboot,theta.fit,theta.predict,err.meas,...) {
    call <- match.call()
    x <- as.matrix(x)
    n <- length(y)
    saveii <- NULL
    fit0 <- theta.fit(x,y,...)
    yhat0 <- theta.predict(fit0,x)
    app.err <- mean(err.meas(y,yhat0))
    err1 <- matrix(0,nrow=nboot,ncol=n)
    err2 <- rep(0,nboot)
    for(b in 1:nboot){
        ii <- sample(1:n,replace = TRUE)
        saveii <- cbind(saveii,ii)
        fit <- theta.fit(x[ii,],y[ii],...)
        yhat1 <- theta.predict(fit,x[ii,])
        yhat2 <- theta.predict(fit,x)  
        err1[b,] <- err.meas(y,yhat2)
        err2[b] <- mean(err.meas(y[ii],yhat1))
    }

    optim <- mean( apply(err1,1,mean)-err2)

    junk <- function(x,i){sum(x==i)}
    e0 <- 0
    for(i in 1:n){
        o <- apply(saveii,2,junk,i)
        if( sum(o==0)==0)
            cat("increase nboot for computation of the .632 estimator",
                fill = TRUE)
        e0 <- e0+ (1/n)*sum(err1[o==0,i])/sum(o==0)
    }
    err.632 <- .368*app.err + .632*e0
    return(list(app.err, 
                optim, 
                err.632, 
                call=call))
}
