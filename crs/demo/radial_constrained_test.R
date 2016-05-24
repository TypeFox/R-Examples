## Here we conduct a hypothesis test using a nonparametric bootstrap
## and the constrained model. num.boot is the number of bootstrap
## replications. D is the test statistic. This code to conduct
## restricted regression splines on evaluation data. Presumes
## continuous regressors, accepts an arbitrary number of regressors,
## and accepts arbitrary derivative restrictions.

## Load libraries

require(crs)
options(crs.messages=FALSE)
require(quadprog)

num.boot <- 999
set.seed(42)

n <- 250

## Test statistic

D <- function(p) {
  n <- length(p)
  ## Smallest weights returned are 6.938894e-18
  if(isTRUE(all.equal(rep(1/n,n),p))) {
    return(0)
  } else {
    return(mean((p-rep(1/n,n))^2))
  }
}

D.boot <- numeric(num.boot)

## Parameters to be set.

x.min <- -5
x.max <- 5

## These will need to be modified if/when you modify Amat and
## bvec. Note the true function lies in the range -0.217 to 1.00. So
## restricting the function to lie strictly within this range is
## imposing an incorrect constraint (for instance, set lower to 0 and
## upper to 0.5 and you would expect to reject the null, while setting
## lower and upper to those below and you would expect to fail to
## reject the null). Note also that imposing nonbinding constraints
## can be checked by inspecting p.hat prior to conducting the
## bootstrap... if they are all equal to 1/n then you have imposed a
## non-binding constraint and in this case the null distribution is
## degenerate.

lower <- -0.217
upper <- 1.00

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

rm(y,x1,x2)

model.unres <- crs(y~x1+x2,
                   basis="auto",
                   data=data.train,
                   nmulti=5)

## If you wish to alter the constraints, you need to modify Amat and
## bvec.

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

D.stat <- D(p.hat)

if(D.stat > 0) {
  
  ## Generate fitted values and data from constrained model
  
  data.trans <- data.frame(y=p.hat*data.train$y,data.train[,2:ncol(data.train)])
  names(data.trans) <- names(data.train) ## Necessary when there is only 1 regressor
  model.res <- crs(y~x1+x2,cv="none",
                   degree=model.unres$degree,
                   segments=model.unres$segments,
                   basis=model.unres$basis,                                  
                   data=data.trans)
  
  yhat <- fitted(model.res)
  model.resid <- data.train$y - fitted(model.unres)
  
  for(b in 1:num.boot) {

    ## Draw a bootstrap resample from the constrained model
    
    y.star <- yhat + sample(model.resid,replace=TRUE)
    
    ## Now compute model for bootstrap resample...
    
    data.boot <- data.frame(y=y.star,data.train[,2:ncol(data.train)])
    names(data.boot) <- names(data.train) ## Necessary when there is only 1 regressor
    model.boot <- crs(y~x1+x2,cv="none",
                      degree=model.unres$degree,
                      segments=model.unres$segments,
                      basis=model.unres$basis,                                  
                      data=data.boot)

    B <- model.matrix(model.boot$model.lm)
    Aymat.res <- t(B%*%solve(t(B)%*%B)%*%t(B))*data.boot$y

    ## Here is Amat
    
    Amat <- cbind(Aymat.res,
                  -Aymat.res)
    
    ## Here is bvec
    
    bvec <- c(rep(lower,n),
              -rep(upper,n))
    
    ## Solve the quadratic programming problem
    
    QP.output <- solve.QP(Dmat=diag(n),dvec=rep(1,n),Amat=Amat,bvec=bvec)

    if(is.nan(QP.output$value)) stop(" solve.QP failed. Try smoother curve (larger bandwidths or polynomial order)")
    
    ## Get the solution
    
    p.hat <- QP.output$solution
    
    D.boot[b] <- D(p.hat)
    
  }
  
  D.boot <- sort(D.boot)
  
}

if(D.stat > 0 ) {
  
  P <- mean(ifelse(D.boot > D.stat, 1, 0))
  ## Check for degenerate case (all bootstrap D statistics identical)
  if(length(unique(D.boot))==1) P <- runif(1) 

  plot(density(D.boot),
       main="Null Distribution",
       sub=paste("Test statistic: ", formatC(D.stat,digits=3,format="g"),
         ", P-value = ", formatC(P,digits=2,format="f"),sep=""))

} else {

  ## If the D statistic is zero then the constraints are non-binding
  ## and the test is degenerate, so issue a warning to this effect
  
  warning("The test statistic is 0, the bootstrap is degenerate, and the constraints are non-binding")
  
}

