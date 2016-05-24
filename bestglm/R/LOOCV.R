`LOOCV` <-
function(X, y){
#LOOCV (equivalent to PRESS)
n <- length(y)
if (n != nrow(X))
    stop(paste("error: X must have ", n, "rows"))
p <- ncol(X) #if zero, handle separately
Xy<- data.frame(as.data.frame(X), y=y)
ans<-lm(y~., data=Xy)
Errs <- (resid(ans)^2) / ((1 - lm.influence(ans)$hat)^2)
sdCV <- sd(Errs)
CVErr <- mean(Errs)
c(CVErr, sdCV)
}

