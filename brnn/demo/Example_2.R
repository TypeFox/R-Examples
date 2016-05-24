#Example 2
#sin wave function, example in the Matlab 2010b demo.

x = seq(-1,0.5,length.out=100)
y = sin(2*pi*x)+rnorm(length(x),sd=0.1)

#With the formula interface
out=brnn(y~x,neurons=3)

#With the default method the call is
#out=brnn(y=y,x=as.matrix(x),neurons=3)

plot(x,y)
lines(x,predict(out),col="blue",lty=2)

legend("bottomright",legend="Fitted model",col="blue",lty=2,bty="n")