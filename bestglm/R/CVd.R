CVd <-
function(X, y, d=ceiling(n*(1-1/(log(n) - 1))), REP=100, family=gaussian, ...){
#CV delete-d method. With sd computation.
#Jun Shao (1993, 1997)
MAXIT <- min(20, REP*0.1) #stop if too many rank degenerate regressions
gaussianQ <- "gaussian"==deparse(substitute(family))
n <- length(y)
if (n != nrow(X))
    stop(paste("error: X must have ", n, "rows"))
p <- ncol(X) #if zero, handle separately
CVErr <- 0
Errs <- numeric(REP)
MSErr <- 0
iREP <- ITER <- 0
while (iREP < REP) {
    iTest <- (1:n)%in%sample(1:n, size=d)
    if (p == 0) #no covariates
        yHat<-mean(y[!iTest])
    else {
        Xj <- X[!iTest,,drop=FALSE]
        yj <- y[!iTest]
        Xyj<- data.frame(as.data.frame(Xj), y=yj)
        if (gaussianQ) {
                ansj<- lm(y~., data=Xyj, ...)
                yHat <- predict(ansj, newdata=X[iTest,,drop=FALSE])
        }
        else {
                ansj<- glm(y~., data=Xyj, family=family, ...)
                yHat <- predict(ansj, newdata=X[iTest,,drop=FALSE],type = "response")
        }        
        if (gaussianQ && ncol(Xyj) > ansj$rank){
            ITER <- ITER+1
             if (ITER > MAXIT)
                stop("Too many rank deficient iterations. Try decreasing d.\nCurrent d =",d)
            next
        }
    }
    iREP <- iREP+1
    Errs[iREP]<-mean((y[iTest]-yHat)^2)
    MSErr <- MSErr + Errs[iREP]       
}
CVErr <- CVErr+MSErr/REP
sdCV <- sd(Errs)
#turned off since Sweave does not accept.
#if (ITER > 0) 
#    warning("Decreasing d or eliminating some variable would improve efficiency.\nThere were ", ITER, " iterations rejected\n due to rank deficient regression matricies.")
c(CVErr, sdCV)
}

