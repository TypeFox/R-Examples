#Example 1 
#Noise triangle wave function, similar to example 1 in Foresee and Hagan (1997)

#Generating the data
x1=seq(0,0.23,length.out=25)
y1=4*x1+rnorm(25,sd=0.1)
x2=seq(0.25,0.75,length.out=50)
y2=2-4*x2+rnorm(50,sd=0.1)
x3=seq(0.77,1,length.out=25)
y3=4*x3-4+rnorm(25,sd=0.1)
x=c(x1,x2,x3)
y=c(y1,y2,y3)


#With the formula interface
out=brnn(y~x,neurons=2)

#With the default S3 method the call is
#out=brnn(y=y,x=as.matrix(x),neurons=2)

plot(x,y,xlim=c(0,1),ylim=c(-1.5,1.5),
     main="Bayesian Regularization for ANN 1-2-1")
lines(x,predict(out),col="blue",lty=2)
legend("topright",legend="Fitted model",col="blue",lty=2,bty="n")

