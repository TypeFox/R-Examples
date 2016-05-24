### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: First Example
###################################################
library(smallarea)
response=c(1,2,3,4,5)  # response vector
designmatrix=cbind(c(1,1,1,1,1),c(1,2,4,4,1),c(2,1,3,1,5)) 
# designmatrix with 5 rows and 3 columns, 
# the first column has all entries equal to one 
sampling.var=c(0.5,0.7,0.8,0.4,0.5) 
# This is the vector of sampling variances
answer=prasadraoest(response,designmatrix,sampling.var)
answer
answer$estimate


###################################################
### code chunk number 2: Second Example
###################################################
response=c(1,2,3,4,5)
designmatrix=cbind(c(1,1,1,1,1),c(1,2,4,4,1),c(2,1,3,1,5))
sampling.var=c(0.5,0.7,0.8,0.4,0.5)
fayherriot(response,designmatrix,sampling.var)


###################################################
### code chunk number 3: Third Example
###################################################
set.seed(55)
response=rnorm(5,3,1.5)
designmatrix=cbind(c(1,1,1,1,1),c(1,2,4,4,1),c(2,1,3,1,5))
sampling.var=c(0.5,0.7,0.8,0.4,0.5)
maximlikelihood(response,designmatrix,sampling.var)


###################################################
### code chunk number 4: Fourth Example
###################################################
set.seed(55)
response=rnorm(5,3,1.5)
designmatrix=cbind(c(1,1,1,1,1),c(1,2,4,4,1),c(2,1,3,1,5))
sampling.var=c(0.5,0.7,0.8,0.4,0.5)
resimaxilikelihood(response,designmatrix,sampling.var,maxiter=100)


###################################################
### code chunk number 5: Fifth Example
###################################################
data=data.frame(response=rnorm(5,3,1.5),
x1=c(1,2,4,4,1),x2=c(2,1,3,1,5),D=c(0.5,0.7,0.8,0.4,0.5))
data
ans=smallareafit(response~D+x1+x2,data,method="FH")
ans1=smallareafit(response~D+x1+x2,data,method="REML")
ans2=smallareafit(response~D+x1+x2,data,method="PR")
ans3=smallareafit(response~D+x1+x2,data,method="ML")
ans # FH method
ans1 # REML method
ans2 # PR method
ans3 # ML method


###################################################
### code chunk number 6: Fifth Example
###################################################
data=data.frame(response=rnorm(5,3,1.5),
x1=c(1,2,4,4,1),D=c(0.5,0.7,0.8,0.4,0.5))
attach(data)
ans=smallareafit(response~D+x1,method="FH")
ans1=smallareafit(response~D+x1,method="REML")
ans2=smallareafit(response~D+x1,method="PR")
ans3=smallareafit(response~D+x1,method="ML")
ans
ans1
ans2
ans3


###################################################
### code chunk number 7: Sixth example
###################################################
set.seed(55)
# the sampling variances
D=c(rep(0.7,3),rep(0.6,3),rep(0.5,3),rep(0.4,3),rep(0.3,3))
# generating the errors
e1=rnorm(3,0,sqrt(D[1]))
e2=rnorm(3,0,sqrt(D[4]))
e3=rnorm(3,0,sqrt(D[7]))
e4=rnorm(3,0,sqrt(D[10]))
e5=rnorm(3,0,sqrt(D[13]))
e=c(e1,e2,e3,e4,e5)
psi=1
# generating the random small area effects
v=rnorm(15,0,psi)
# response  
y=v+e
data1=data.frame(response=y,D=D)
head(data1)
fit1.pr=smallareafit(response~D,data1,method="PR")
fit1.pr


