monreg.wrapper <- function(x, y, hr.vals=seq(0.01, 0.40, l=40), x0, folds=5){
    ##`monreg.wrapper` <- function(x,y,hr.vals=seq(0.01,0.40,l=40),folds=5){
    ##require(bootstrap)
    ##require(monreg)
    mae <- numeric(length(hr.vals))
    theta.fit <- function(x,y) list(x=x,y=y)
    theta.predict <- function(fit,x) monreg(fit$x,fit$y,hd=hd,hr=hr,t=x)$estimation
    ii <- 0
    for(hr in hr.vals){
        ii <- ii+1
        hd <- hr^2
        y.cv <- crossval(x,y,theta.fit,theta.predict,ngroup=folds)$cv.fit
        mae[ii] <- mean(abs(y-y.cv))
    }
    hr <- hr.vals[which.min(mae)]
    hd <- hr^2
    if (missing(x0)) x0 <- x
    xstd <- (x0-min(x0))/diff(range(x0))
    out <- monreg(x,y,hd=hd,hr=hr,t=xstd)
    return(out)
}

