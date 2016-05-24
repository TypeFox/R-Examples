CVDH <-
function(X, y, K=10, REP=1){
#K-fold adjusted CV. sd also computed. 
#Algorithm 6.5 page 295
#Reference: Davison and Hinkley, 
#           Boostrap Methods and Their Application
n <- length(y)
if (n != nrow(X))
    stop(paste("error: X must have ", n, "rows"))
p <- ncol(X) #if zero, handle separately
Xy<- data.frame(as.data.frame(X), y=y)
DFF<-mean(resid(lm(y~., data=Xy))^2)
CVAdj<-0
varCV<-0
for (iRep in 1:REP){
    fold <- sample(rep(1:K,length=n))
    f<-tabulate(fold)/n
    SumSqErr <- 0
    Errs <- numeric(K)
    DFFk<-numeric(K)
    for (j in 1:K) {
        iTest <- fold==j
        if (p == 0) {#no covariates
            yHat<-mean(y[!iTest])
            DFFk[j]<-mean((y-yHat)^2)
            }
        else {
            Xj <- X[!iTest,,drop=FALSE]
            yj <- y[!iTest]
            Xyj<- data.frame(as.data.frame(Xj), y=yj)
            ansj<- lm(y~., data=Xyj)
            yHat <- predict(ansj, newdata=X[iTest,,drop=FALSE])
            DFFk[j]<-mean((y-predict(ansj, newdata=X))^2)
        }
        Errs[j]<-mean((y[iTest]-yHat)^2)
        SumSqErr <- SumSqErr + sum((y[iTest]-yHat)^2)
        }
    varCV <- varCV+var(Errs)
    CVErr <- SumSqErr/n
    CVAdj <- CVAdj + CVErr + DFF - sum(f*DFFk)
    }
CVErr<-CVAdj/REP
sdCV<-sqrt(varCV/REP)
c(CVErr, sdCV)
}

