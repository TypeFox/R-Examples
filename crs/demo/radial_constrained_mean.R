## Code to conduct restricted regression splines on evaluation
## data. Presumes continuous regressors, accepts an arbitrary number
## of regressors, and accepts arbitrary derivative restrictions.

## IMPORTANT NOTE - the code that follows is only valid for the tensor
## basis (basis="tensor")

## Load libraries

require(crs)
require(quadprog)

## Parameters to be set.

set.seed(42)

n <- 1000
n.eval <- 50

x.min <- -5
x.max <- 5

## These will need to be modified if/when you modify Amat and bvec

lower <- 0
upper <- 0.5

## IMPORTANT - you must be careful to NOT read data from environment -
## this appears to work - create a data frame.

## IMPORTANT - code that follows presumes y is the first variable in
## the data frame and all remaining variables are regressors used for
## the estimation.

## Generate a DGP, or read in your own data and create y, x1,
## etc. When you change this by adding or removing variables you need
## to change `data', `rm(...)', and `bw <- ...'. After that all code
## will need no modification.

x1 <- runif(n,x.min,x.max)
x2 <- runif(n,x.min,x.max)

y <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2) + rnorm(n,sd=.1)

data.train <- data.frame(y,x1,x2)

x1.seq <- seq(min(x1),max(x1),length=n.eval)
x2.seq <- seq(min(x2),max(x2),length=n.eval)

rm(y,x1,x2)

data.eval <- data.frame(y=0,expand.grid(x1=x1.seq,x2=x2.seq))

model.unres <- crs(y~x1+x2,
                   data=data.train,
                   nmulti=5,
                   basis="tensor")

summary(model.unres)

## If you wish to alter the constraints, you need to modify Amat and
## bvec.

## Generate the estimated model computed for the training data. Note -
## we need to premultiply the weights by n and each column must be
## multiplied by y

B <- model.matrix(model.unres$model.lm)
Aymat.res <- t(B%*%solve(t(B)%*%B)%*%t(B))*data.train$y

## Here is Amat

Amat <- cbind(Aymat.res,
              -Aymat.res)

rm(Aymat.res)

## Here is bvec

bvec <- c(rep(lower,n),
          -rep(upper,n))

## Solve the quadratic programming problem

QP.output <- solve.QP(Dmat=diag(n),dvec=rep(1,n),Amat=Amat,bvec=bvec)

if(is.nan(QP.output$value)) stop(" solve.QP failed. Try smoother curve (larger bandwidths or polynomial order)")

## No longer needed...

rm(Amat,bvec)

## Get the solution

p.hat <- QP.output$solution

## Now estimate the restricted model

data.trans <- data.frame(y=p.hat*data.train$y,data.train[,2:ncol(data.train),drop=FALSE])

model.res <- crs(y~x1+x2,cv="none",
                 degree=model.unres$degree,
                 segments=model.unres$segments,
                 basis=model.unres$basis,                                  
                 data=data.trans,
                 deriv=1)

## That's it!

## Create a 3D perspective plot of the constrained and unconstrained
## surfaces

fitted.unres <- matrix(predict(model.unres,newdata=data.eval), n.eval, n.eval)
fitted.res <- matrix(predict(model.res,newdata=data.eval), n.eval, n.eval)

zlim <- c(min(fitted.unres,fitted.res),max(fitted.unres,fitted.res))

par(mfrow=c(1,2))

persp(x1.seq, x2.seq,
      fitted.unres,
      main="Unconstrained Regression Spline",
      col="lightblue",
      ticktype="detailed", 
      ylab="X2",
      xlab="X1",
      zlim=zlim,
      zlab="Conditional Expectation",
      theta=300,
      phi=30)

persp(x1.seq, x2.seq,
      fitted.res,
      main="Constrained Regression Spline",
      sub="0 <= g(x1,x2) <= 1/2",
      col="lightblue",
      ticktype="detailed", 
      ylab="X2",
      xlab="X1",
      zlim=zlim,
      zlab="Conditional Expectation",
      theta=300,
      phi=30)

par(mfrow=c(1,1))


