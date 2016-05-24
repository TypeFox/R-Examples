plotpi <-
function(k, global=0.5){
  x = seq(0, 1, 0.01)
  y = sapply(x, function(pr) piFunction(pr, kparam=k, qparam=global))
  plot(x, y, type='l', main=paste('k =',round(k,3)), 
       xlab='Fitted probability', ylab='Treatment probability', 
       xaxt='n', yaxt='n', bty='n')
  axis(side=1, lwd=0, lwd.ticks=1, at=seq(0,1,.25))
  axis(side=2, lwd=0, lwd.ticks=1, at=seq(0,1,.25))
  abline(h=c(0,global), lty=1:2, col=c('black','grey'))
  abline(v=c(0,global), lty=1:2, col=c('black','grey'))
  lines(x, y)
}
