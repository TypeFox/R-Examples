
library(marqLevAlg)

## first test

f1 <- function(b){	
	return(4*(b[1]-5)^2+(b[2]-6)^2)	
}

## derivate
gr <- function(b){
	return(c(8*(b[1]-5),2*(b[2]-6)))
}

## hessian
hes <- function(b){
	return(c(-8,0,-2))
}

### without parameters
test.marq <- marqLevAlg(m=2,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1)
test.marq

### without length of parameters
test.marq <- marqLevAlg(b=c(8,9),maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1)
test.marq

### with parameters and length of parameters
test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1)

### with gradiant
test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1,gr=gr)

### with hessian
test.marq <- marqLevAlg(b=c(8,9),m=2,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f1,gr=gr,hess=hes)
test.marq


## second test

f2 <- function(b){(b[1]+10*b[2])^2+5*(b[3]-b[4])^2+(b[2]-2*b[3])^4+10*(b[1]-b[4])^4}


### without parameters
test.marq2 <- marqLevAlg(m=4,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f2)
test.marq2

### without length of parameters
test.marq2 <- marqLevAlg(b=c(3,-1,0,1),maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f2)
test.marq2

### with parameters and length of parameters
test.marq2 <- marqLevAlg(b=c(3,-1,0,1),m=4,maxiter=100,epsa=0.001,epsb=0.001,epsd=0.001,fn=f2)
test.marq2


## third test
f<-function(x){(x[1]^42*(x[2]-3)^22)} 
test.marq3 <- marqLevAlg(b=c(10,10),maxiter=10000,epsa=0.1^5,epsb=0.1^5,epsd=0.1^5,fn=f)
test.marq3


## test 4
f<-function(x){x^120} 
test.marq4 <- marqLevAlg(b=10,maxiter=10000,epsa=0.1^5,epsb=0.1^5,epsd=0.1^5,fn=f)
test.marq4



