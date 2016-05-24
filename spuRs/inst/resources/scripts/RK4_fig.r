par(mfrow=c(2,2))
for (h in c(10, 5, 2)) {
  y <- RK4(LV, 0, c(100, 50), h, 1000/h,
           a=0.05, K=Inf, b=0.0002, g=0.8, e=0.03)
  plot(y[,1], y[,2], type='p',
       xlab='prey', ylab='pred', main=paste('RK4, h =',h))
}
plot(seq(0, 1000, h), y[,1], type='l', xlab='time',
     ylab='prey solid pred dashed', main=paste('RK4, h =',h))
lines(seq(0, 1000, h), y[,2], lty=2)
par(mfrow=c(1,1))
