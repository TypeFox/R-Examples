## This illustration considers the `radial function' and plots the
## results in a 3D perspective plot.

require(crs)

set.seed(42)

## Interactively request number of observations, whether to do NOMAD
## or exhaustive search, and if NOMAD the number of multistarts

n <- as.numeric(readline(prompt="Input the number of observations desired: "))
cv <- as.numeric(readline(prompt="Input the cv method (0=nomad, 1=exhaustive): "))
cv <- ifelse(cv==0,"nomad","exhaustive")
if(cv=="nomad") nmulti <- as.numeric(readline(prompt="Input the number of multistarts desired (e.g. 10): "))
num.eval <- as.numeric(readline(prompt="Input the number of evaluation observations desired (e.g. 50): "))

x1 <- runif(n,-5,5)
x2 <- runif(n,-5,5)

dgp <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2)

y <- dgp + rnorm(n,sd=.1)

model <- crs(y~x1+x2,
             cv=cv,
             complexity="degree-knots",
             knots="uniform",
             deriv=1,
             cv.func="cv.aic",
             nmulti=nmulti)

summary(model)

# Perspective plot
x1.seq <- seq(min(x1),max(x1),length=num.eval)
x2.seq <- seq(min(x2),max(x2),length=num.eval)
x.grid <- expand.grid(x1.seq,x2.seq)
newdata <- data.frame(x1=x.grid[,1],x2=x.grid[,2])
z <- matrix(predict(model,newdata=newdata),num.eval,num.eval)
persp(x=x1.seq,y=x2.seq,z=z,
      xlab="X1",ylab="X2",zlab="Y",
      ticktype="detailed",      
      border="red",
      main="Conditional Mean",
      theta=45,phi=45)

## Perspective plot - derivative wrt x1
z <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,1],num.eval,num.eval)

persp(x=x1.seq,y=x2.seq,z=z,
      xlab="X1",ylab="X2",zlab="Y",
      ticktype="detailed",      
      border="red",
      main="d g(x1,x2)/d x1 (x2=med(x2))",
      theta=45,phi=45)

## Perspective plot - derivative wrt x2
z <- matrix(attr(predict(model,newdata=newdata),"deriv.mat")[,2],num.eval,num.eval)

persp(x=x1.seq,y=x2.seq,z=z,
      xlab="X1",ylab="X2",zlab="Y",
      ticktype="detailed",      
      border="red",
      main="d g(x1,x2)/d x2 (x1=med(x1))",
      theta=45,phi=45)
