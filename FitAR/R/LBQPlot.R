`LBQPlot` <-
function(res, lag.max=30, StartLag=k+1, k=0, SquaredQ=FALSE){
stopifnot(k>=0, lag.max-StartLag>0, length(res)>lag.max)
ans<-LjungBoxTest(res, k=k, lag.max=lag.max, StartLag=StartLag, SquaredQ=FALSE)
m <- ans[,"m"]
pv <- ans[,"pvalue"]
plot(m, pv, xlab="lag", ylab="p-value", ylim=c(0,1), main="Ljung-Box Test")
abline(h=0.05, col="red", lty=2)
}

