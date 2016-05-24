loess.surf<-function(Y, X, span = 0.75, degree = 1, family = "gaussian", phi = 20, theta = 50, xlab = "X", ylab = "Y", 
zlab = "Fit", line.col = 1, line.type = 1, scale = TRUE, duplicate = "error", expand = 0.5, ...){
lom<-loess(Y ~ X, span = span, degree = degree, family = family)
r<-ncol(lom$x)
yhat<-lom$fitted
if(r>2)stop("Number of predictors must be < 2")
if(r==1){
plot(X,Y,xlab=xlab,ylab=ylab,...)
d<-10^round(log(min(X),10),0)
xv<-seq(min(X),max(X),d/100)
yv<-predict(lom, data.frame(X=xv))
lines(xv,yv,col=line.col,lty=line.type)
}
if(r==2){##Thanks to Wilcox(2005)
iout<-c(1:length(yhat))
nm1<-length(yhat)-1
for(i in 1:nm1){
ip1<-i+1
for(k in ip1:length(yhat))if(sum(X[i,]==X[k,])==2)iout[k]<-0
}
yhat1 <- yhat[iout >= 1] 
X1 <- X[iout >= 1,]
interp <- NA
rm(interp)
fit <- akima::interp(X1[,1], X1[,2], yhat1, duplicate = duplicate)
persp(fit, theta = theta, phi = phi, xlab = xlab, ylab = ylab, zlab = zlab, scale = scale, expand = expand)
}
}