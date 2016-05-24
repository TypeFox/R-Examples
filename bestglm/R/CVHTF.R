CVHTF <-
function(X, y, K=10, REP=1, family=gaussian, ...){
#CV K-fold with replication
#compute standard devations as per HTF
n <- length(y)
gaussianQ <- "gaussian"==deparse(substitute(family))
if (n != nrow(X))
    stop(paste("error: X must have ", n, "rows"))
p <- ncol(X) #if zero, handle separately
CVErr <- 0
sumVarCV <- 0
for (iREP in 1:REP){
    fold <- sample(rep(1:K,length=n)) #depends on RNG
    SumSqErr <- 0
    Errs <- numeric(K)
    for (k in 1:K) {
        iTest <- fold==k
        if (p == 0) #no covariates
            yHat<-mean(y[!iTest])
        else {
            Xk <- X[!iTest,,drop=FALSE]
            yk <- y[!iTest]
            Xyk<- data.frame(as.data.frame(Xk), y=yk)
            if (gaussianQ) {
                ansj<- lm(y~., data=Xyk, ...)
                yHat <- predict(ansj, newdata=X[iTest,,drop=FALSE])
                }
            else {
                ansj<- glm(y~., data=Xyk,family=family, ...)
                yHat <- predict(ansj, newdata=X[iTest,,drop=FALSE],type = "response")
                }
        }
        Errs[k]<-mean((y[iTest]-yHat)^2)
        SumSqErr <- SumSqErr + sum((y[iTest]-yHat)^2)
        }
#CVerr: average EPE overall all K test sets
        CVErr <- CVErr+SumSqErr/n
#Changjiang. divide by K. Oct 29, 09. 
#Brillant!! Yes, we need the variance of the average!
        sumVarCV <- sumVarCV + var(Errs)/K
    }
#average over number of replications
    CVErr<-CVErr/REP
    sdCV<-sqrt(sumVarCV/REP)
    c(CVErr, sdCV)
}

