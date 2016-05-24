#
#
#    Script for adjustment of lambda
#
#
library(aws)
############################################################
#
#                     univariate   
#
############################################################

#  Gaussian models
y <- rnorm(10000)

# local const    alpha < 0.01
yhat <- aws(y,hmax=10000,ladjust=1,u=0,testprop=TRUE,graph=TRUE)#11.3

# local polynomial degree 1   alpha=0.05

yhatInf <-  lpaws(y,degree=1,hmax=1000,ladjust=1e10,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)
yhat <-  lpaws(y,degree=1,hmax=1000,ladjust=1,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)#6
max(yhat@mae/yhatInf@mae-1)

# local polynomial degree 0   alpha=0.05

yhatInf <-  lpaws(y,degree=0,hmax=1000,ladjust=1e10,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)
yhat <-  lpaws(y,degree=0,hmax=1000,ladjust=1,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)#6
max(yhat@mae/yhatInf@mae-1)

# local polynomial degree 2   alpha=0.05

yhatInf <-  lpaws(y,degree=2,hmax=1000,ladjust=1e10,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)
yhat <-  lpaws(y,degree=2,hmax=1000,ladjust=1,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)#10.4
max(yhat@mae/yhatInf@mae-1)

# irregular design alpha=0.08

x <- sort(runif(10000,0,10000))
yhatInf <-  aws.irreg(y,x,hmax=1000,nbins=1000,ladjust=1e10,graph=TRUE)
yhat <-  aws.irreg(y,x,hmax=1000,nbins=1000,ladjust=1,graph=TRUE)#14.2
sd(yhat@theta)/sd(yhatInf@theta)-1

# local const with variance model    alpha < 0.01

yhatInf <-  aws.gaussian(y,hmax=1000,ladjust=1e10,u=0,graph=TRUE,homogen=FALSE)
yhat <-  aws.gaussian(y,hmax=1000,ladjust=1,u=0,graph=TRUE,homogen=FALSE)#11.3
max(yhat@mae/yhatInf@mae-1)

# local const  Bernoulli  alpha = 0.07

y <- rbinom(10000,1,.05) # extreme values are more critical
yhat <- aws(y,hmax=10000,ladjust=1,u=0.05,family="Bernoulli",testprop=TRUE,graph=TRUE)#9.6


# local const  Poisson  alpha = 0.01

y <- rpois(10000,1) 
yhat <- aws(y,hmax=10000,ladjust=1,u=1,family="Poisson",testprop=TRUE,graph=TRUE)#9.6


# local const  Exponential alpha < 0.01
y <- rexp(10000,1)
yhat <- aws(y,hmax=4000,ladjust=1,u=1,family="Exponential",testprop=TRUE,graph=TRUE)#14.2


# local const  Volatility alpha =0.02
y <- rnorm(10000)
yhat <- aws(y,hmax=4000,ladjust=1,u=1,family="Volatility",testprop=TRUE,graph=TRUE)#10

# local const  Variance alpha = 0.05
y <- rnorm(10000)^2
yhat <- aws(y,hmax=4000,ladjust=1,u=1,family="Variance",testprop=TRUE,graph=TRUE,shape=1)#12.8

############################################################
#
#                     bivariate   
#
############################################################

#  Gaussian models
y <- matrix(rnorm(512^2),512,512)

# local const    alpha = 0.07
yhat <- aws(y,hmax=10,ladjust=1,u=0,testprop=TRUE,graph=TRUE)#6.1

# local polynomial degree 1   alpha=0.07

yhatInf <-  lpaws(y,degree=1,hmax=15,ladjust=1e10,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)
yhat <-  lpaws(y,degree=1,hmax=15,ladjust=1,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)#11.3
max(yhat@mae/yhatInf@mae-1)

# local polynomial degree 2   alpha=0.09

yhatInf <-  lpaws(y,degree=2,hmax=20,ladjust=1e10,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)
yhat <-  lpaws(y,degree=2,hmax=20,ladjust=1,u=0,graph=TRUE,homogen=FALSE,earlystop=FALSE)# 27
max(yhat@mae/yhatInf@mae-1)

# irregular design alpha=0.05

y <- rnorm(10000)
x <- cbind(runif(10000,0,1),runif(10000,0,1))
yhatInf <-  aws.irreg(y,x,hmax=25,nbins=200,ladjust=1e10,graph=TRUE)
yhat <-  aws.irreg(y,x,hmax=25,nbins=200,ladjust=1,graph=TRUE)#9.1
sd(yhat@theta)/sd(yhatInf@theta)-1

# local const with variance model    alpha = 0.025

y <- matrix(rnorm(512^2),512,512)
yhatInf <-  aws.gaussian(y,hmax=10,ladjust=1e10,u=0,graph=TRUE,homogen=FALSE)
yhat <-  aws.gaussian(y,hmax=10,ladjust=1,u=0,graph=TRUE,homogen=FALSE)#7.6
max(yhat@mae/yhatInf@mae-1)

# local const  Bernoulli  alpha = 0.025

y <- matrix(rbinom(512^2,1,.05),512,512) # extreme values are more critical
yhat <- aws(y,hmax=10,ladjust=1,u=0.05,family="Bernoulli",testprop=TRUE,graph=TRUE)#7.6


# local const  Poisson  alpha = 0.02

y <- matrix(rpois(512^2,1),512,512) 
yhat <- aws(y,hmax=10,ladjust=1,u=1,family="Poisson",testprop=TRUE,graph=TRUE)#7.6


# local const  Exponential alpha = 0.03
y <- matrix(rexp(512^2,1),512,512) 
yhat <- aws(y,hmax=10,ladjust=1,u=1,family="Exponential",testprop=TRUE,graph=TRUE)#6.8


# local const  Volatility alpha =0.05
y <- matrix(rnorm(512^2),512,512) 
yhat <- aws(y,hmax=10,ladjust=1,u=1,family="Volatility",testprop=TRUE,graph=TRUE)#6.1

# local const  Variance alpha = 0.04
y <- matrix(rnorm(512^2)^2,512,512)
yhat <- aws(y,hmax=10,ladjust=1,u=1,family="Variance",testprop=TRUE,graph=TRUE,shape=1)#6.1

############################################################
#
#                     3D
#
############################################################

#  Gaussian models
y <- array(rnorm(64^3),c(64,64,64))

# local const    alpha = 0.05
yhat <- aws(y,hmax=5,ladjust=1,u=0,testprop=TRUE,graph=TRUE)#6.2


# local const  Bernoulli  alpha = 0.02

y <- array(rbinom(64^3,1,.05),c(64,64,64)) # extreme values are more critical
yhat <- aws(y,hmax=5,ladjust=1,u=0.05,family="Bernoulli",testprop=TRUE,graph=TRUE)#6.9


# local const  Poisson  alpha = 0.02

y <- array(rpois(64^3,1),c(64,64,64)) 
yhat <- aws(y,hmax=5,ladjust=1,u=1,family="Poisson",testprop=TRUE,graph=TRUE)#6.9


# local const  Exponential alpha = 0.04
y <- array(rexp(64^3,1),c(64,64,64)) 
yhat <- aws(y,hmax=5,ladjust=1,u=1,family="Exponential",testprop=TRUE,graph=TRUE)#6.1


# local const  Volatility alpha =0.045
y <- array(rnorm(64^3),c(64,64,64)) 
yhat <- aws(y,hmax=5,ladjust=1,u=1,family="Volatility",testprop=TRUE,graph=TRUE)#6.1

# local const  Variance alpha = 0.04
y <- array(rnorm(64^3)^2,c(64,64,64))
yhat <- aws(y,hmax=5,ladjust=1,u=1,family="Variance",testprop=TRUE,graph=TRUE,shape=1)#6.1

