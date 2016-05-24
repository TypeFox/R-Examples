library(exact2x2)

## check error found by Hamilton Gross
## prior to version 1.1-1.0 the following gave an error:
x<-factor(c(0,0,1,1,1))
y<-factor(c(0,0,0,0,0),levels=c(0,1))
mcnemar.exact(x,y)