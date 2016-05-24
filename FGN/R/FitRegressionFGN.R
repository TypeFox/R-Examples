`FitRegressionFGN` <-
function(X,y){
    MaxIter<-10
    n<-length(y)
    if (nrow(X)!=n) 
            stop("error: number of rows in X != length of y")
    res<-lsfit(X,y,intercept=FALSE)$residuals
    out<-GetFitFGN(res)
    LLprev<-out$loglikelihood
    H0<-out$H
    etol<-10
    iter<-0
    while(iter<MaxIter&&etol>0.001){
            r<-acvfFGN(H0, n-1)
            Ginv<-TrenchInverse(toeplitz(r))
            XpGinv<-crossprod(X,Ginv)
            XpXinv<-solve(XpGinv%*%X)
            alpha<-crossprod(XpXinv,XpGinv)%*%y
            res<-y-X%*%alpha
            out<-GetFitFGN(res, MeanZeroQ=TRUE)
            LL<-out$loglikelihood
            etol<-LL-LLprev
            H0<-out$H
            iter<-iter+1
            }
    list(loglikelihood=LL,H=H0,coef=alpha)
}

