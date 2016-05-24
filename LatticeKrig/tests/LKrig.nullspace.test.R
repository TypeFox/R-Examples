# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library(LatticeKrig)
 options( echo=FALSE)
 test.for.zero.flag<- 1
  set.seed(123)
data(ozone2)

 x<-ozone2$lon.lat
 y<- ozone2$y[16,]
good<- !is.na( y)
x<- x[good,]
y<- y[good]
xnew<- ozone2$lon.lat[!good,]

obj<- LKrig( x,y, NC=4, nlevel=2, a.wght=5, nu=1, lambda=.1,
             fixedFunction=NULL)
hold<- predict( obj, xnew=x)

test.for.zero( hold, obj$fitted.values, tag="predict and fitted values")

Sigma<- obj$rho.MLE*LKrig.cov( x,x, LKinfo=obj$LKinfo) +
                  (obj$sigma.MLE^2) * diag( 1, length(y))
S0<- obj$rho.MLE*LKrig.cov( x,xnew, LKinfo=obj$LKinfo)

hold<- predict( obj, xnew=xnew)
hold2<- t(S0)%*%solve( Sigma)%*%y
test.for.zero( hold, hold2, tag="predicts at new values")

#hold<- predictSE( obj, xnew=xnew)
#hold2<-  obj$rho.MLE  - diag(  t(S0)%*%solve( Sigma) %*% S0 )
#test.for.zero( hold, hold2, tag=" SE at new values")

