library('skellam')
lambda1 = 3
lambda2 = 2
par(mfrow=c(1,2))

x = rskellam(100,lambda1, lambda2)
hist(x,prob=TRUE, main='skellam')

xseq = seq(min(x)-1, max(x)+1)
lines(xseq, dskellam(xseq, lambda1, lambda2),col='red')
legend('topright', fill=c('black','red'), legend=c('hist','dens'))


pseq = seq(0,1,len=1000)
plot(qskellam(pseq, lambda1, lambda2), pseq, type='l', xlab='x', ylab='quant')
lines(xseq, pskellam(xseq, lambda1, lambda2),col='red')
legend('topright', fill=c('black','red'), legend=c('qskellam','pskellam'))
