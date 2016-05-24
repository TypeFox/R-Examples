# seasrescheck.R
# seasonal residual checks
# April 2009

seasrescheck<-function(res){
op <- par(no.readonly = TRUE) # the whole list of settable par's.
par(mfrow=c(2,2),lwd=1)
# histogram
par(mai=c(0.7,0.7,0.1,0.1)) # c(bottom, left, top, right)
hist(res,col='gray',main='',xlab='')
# scatter plot
plot(res,type='p',main='',xlab='')
lines(c(1,length(res)),c(0,0),lty=2)
# autocovariance
acf(res,type='correlation',main='',ylab='Autocorrelation')
# cumulative periodogram
cpgram(res,main='')
par(mfrow=c(1,1))
# box plot?
par(op) # restore graphic settings
}