# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

#library( fields, lib.loc="lib.test")

library( fields)
options(echo=FALSE)
test.for.zero.flag<- 1

DD<- cbind( seq(.01,2,,50))
look2<- RadialBasis(DD,  dimension=2,M=3,derivative=1) 

look1<- ( RadialBasis(DD+1e-5, dimension=2,M=3,derivative=0 )
-  RadialBasis(DD-1e-5, dimension=2,M=3,derivative=0))/2e-5

test.for.zero( look1, look2,tol=1e-6, tag="radial basis function exact" )


  set.seed( 234)
  x<- matrix( runif(10), ncol=2)
  ctest<- rep(0,5)
  ctest[3]<- 1
  stationary.cov( x,x, Covariance="RadialBasis", dimension=2,M=3,derivative=0)-> look0
  RadialBasis( rdist(x,x), dimension=2,M=3,derivative=0)-> sanity.look
  test.for.zero( look0, sanity.look, tag="sanity test of stationary.cov with RadialBasis")

  Rad.cov(x,x,p= (2*3 -2))-> look1
  test.for.zero( sanity.look, look1, tag="sanity test of Rad.cov")

  sanity.look%*% ctest->look0
  stationary.cov( x,x, Covariance="RadialBasis", dimension=2,M=3,
                                  derivative=0, C=ctest)-> look
  test.for.zero( look0, look, tag="stat.cov Radbas C multiply")
  Rad.cov(x,x,p= (2*3 -2), C=ctest)-> look1
  test.for.zero( look0, look1, tag="Rad.cov C multiply")


############################ end of radial basis


DD<- cbind( seq(.01,2,,50))
look2<- Wendland(DD, theta=1.0, dimension=2,k=3,derivative=1) 

look1<- (Wendland(DD+1e-5, theta=1.0, dimension=2,k=3)
- Wendland(DD-1e-5, theta=1.0, dimension=2,k=3))/2e-5

test.for.zero( look1, look2,tol=1e-6)



look2<- Wendland(DD, theta=1.5, dimension=2,k=3,derivative=1) 

look1<- (Wendland(DD+1e-5, theta=1.5, dimension=2,k=3)
- Wendland(DD-1e-5, theta=1.5, dimension=2,k=3))/2e-5

test.for.zero( look1, look2,tol=1e-6, tag="Wendland exact")

x<- seq( -1,1,,5)

ctest<- rep(0,5)
ctest[3]<- 1

wendland.cov( x,x, k=2, theta=.75)-> look0
Wendland( rdist(x,x)/.75, k=2, dimension=1)-> sanity.look
test.for.zero( look0, sanity.look)

look0%*% ctest->look0

wendland.cov( x,x, k=2, theta=.75, C=ctest, derivative=0)-> look

test.for.zero( look0, look, tag="Wendland C multiply")


wendland.cov( x,x, k=2, theta=1.0, C=ctest, derivative=1)-> look

wendland.cov( x+1e-5, x, k=2, theta=1.0, C=ctest)-
wendland.cov( x-1e-5, x, k=2, theta=1.0, C=ctest)-> look2
look2<- look2/2e-5
 
test.for.zero( look, look2,tol=1e-7, tag="Wendland.cov theta=1.0")


wendland.cov( x,x, k=2, theta=.75, C=ctest, derivative=1)-> look
wendland.cov( x+1e-5, x, k=2, theta=.75, C=ctest)-
wendland.cov( x-1e-5, x, k=2, theta=.75, C=ctest)-> look2
look2<- look2/2e-5
test.for.zero( look, look2,tol=1e-7, tag="Wendland.cov theta=.75")


stationary.cov( x,x, Covariance="Wendland", dimension=1,
                k=2, theta=1.0, C=ctest, derivative=0)-> look
look0<- Wendland( rdist(x,x), k=2, dimension=1)%*%ctest
test.for.zero( look0, look, tag="stationary.cov and exact C multiply for Wendland")

wendland.cov( x,x, k=2,C=ctest, theta=1.0)-> look
look0<- Wendland( rdist(x,x), k=2, dimension=1)%*%ctest
test.for.zero( look0, look, tag="  Wendland C multiply")

####### 2 -d quadratic surface 

x<- make.surface.grid( list(x=seq( -1,1,,20), y=seq( -1,1,,20)))
y<- (.123*x[,1] + .234*x[,2])
obj<- mKrig( x,y, lambda=0, cov.function="wendland.cov", k=3, theta=.4)

xp<- make.surface.grid( list(x=seq(-.5,.5,,24),y= seq( -.5,.5,,24)) )
predict( obj, xp, derivative=1)-> outd
test.for.zero( outd[,1],.123, tag="2-d derivs from wend.cov/mKrig")
test.for.zero( outd[,2],.234)


#%%%%%%%% repeat to check derivatives in stationary.cov

x<- make.surface.grid( list(x=seq( -1,1,,20), y=seq( -1,1,,20)))
y<- (.123*x[,1] + .234*x[,2])
obj<- mKrig( x,y, lambda=0, cov.function="stationary.cov",
            cov.args=list(k=3, theta=.2, dimension=2, Covariance="Wendland"))

xp<- make.surface.grid( list(x=seq(-.5,.5,,24),y= seq( -.5,.5,,24)) )
predict( obj, xp, derivative=1)-> outd
test.for.zero( outd[,1],.123, tag="2-d derivs from stationary-wend/mKrig")
test.for.zero( outd[,2],.234)


############## quadratic surface
x<- make.surface.grid( list(x=seq( -1,1,,20), y=seq( -1,1,,20)))
y<- (x[,1]**2 - 2* x[,1]*x[,2] +  x[,2]**2)/2

############## wendland.cov
obj<- mKrig( x,y, lambda=0, cov.function="wendland.cov", k=3, theta=.8)
xp<- make.surface.grid( list(x=seq(-.5,.5,,24),y= seq( -.5,.5,,24)) )
true<- cbind( xp[,1] -  xp[,2], xp[,2]- xp[,1])
############## wendland.cov
predict( obj, xp, derivative=1)-> outd
rmse<-sqrt(mean((true[,1] - outd[,1])**2))/sqrt(mean(true[,1]**2))
test.for.zero( rmse,0, tol=5e-3,relative=FALSE, tag="wendland.cov quad 2-d")

############## stationary cov
x<- make.surface.grid( list(x=seq( -1,1,,20), y=seq( -1,1,,20)))
y<- (x[,1]**3 +  x[,2]**3)
obj<- mKrig( x,y, lambda=0, cov.function="stationary.cov",
            cov.args=list(k=3, theta=.8, dimension=2, Covariance="Wendland"))

xp<- make.surface.grid( list(x=seq(-.5,.5,,24),y= seq( -.5,.5,,24)) )
true<- cbind( 3*xp[,1]**2 , 3*xp[,2]**2)
predict( obj, xp, derivative=1)-> outd2
rmse<-sqrt(mean((true[,1] - outd2[,1])**2))/sqrt(mean(true[,1]**2))
test.for.zero( rmse,0, tol=1e-2,relative=FALSE,
                        tag="stationary.cov/Wendland cubic 2-d")

############## stationary cov  with radial basis
x<- make.surface.grid( list(x=seq( -1,1,,20), y=seq( -1,1,,20)))
y<- (x[,1]**3 +  x[,2]**3)
obj<- Krig( x,y,  cov.function="stationary.cov", m=3,
            cov.args=list(M=3, dimension=2, Covariance="RadialBasis"))

xp<- make.surface.grid( list(x=seq(-.5,.5,,24),y= seq( -.5,.5,,24)) )
true<- cbind( 3*xp[,1]**2 , 3*xp[,2]**2)
predictDerivative.Krig( obj, xp)-> outd2
look<- as.surface( xp, outd2[,1])
rmse<-sqrt(mean((true[,1] - outd2[,1])**2))/sqrt(mean(true[,1]**2))
test.for.zero( rmse,0, tol=5e-3,relative=FALSE,
                        tag="stationary.cov/Wendland cubic 2-d")


#########################
  x<- make.surface.grid( list(x=seq( -1,1,,20), y=seq( -1,1,,20)))
  y<- (x[,1]**3 +  x[,2]**3)

  obj<- mKrig( x,y, lambda=0, cov.function="wendland.cov", k=3,  V=diag(c( 1.1,1.1) ))
  xp<- make.surface.grid( list(x=seq(-.5,.5,,24),y= seq( -.5,.5,,24)) )
  predict( obj, xp, derivative=1)-> outd
  true<- cbind( 3*xp[,1]**2 , 3*xp[,2]**2)
  rmse<-sqrt(mean((true[,1] - outd[,1])**2)/mean(true[,1]**2))
  test.for.zero( rmse,0, tol=5e-3,relative=FALSE)

  obj<- Tps( x,y,lambda=0)
  predictDerivative.Krig( obj, xp, derivative=1)-> outd
  look<- as.surface( xp, outd[,1])
  rmse<-sqrt(mean((true[,1] - outd[,1])**2)/mean(true[,1]**2))
  test.for.zero( rmse,0, tol=2e-4,relative=FALSE, tag="Tps derivative x1")
  rmse<-sqrt(mean((true[,2] - outd[,2])**2)/mean(true[,2]**2))
  test.for.zero( rmse,0, tol=2e-4,relative=FALSE, tag="Tps derivative x2")

cat("done with dervative tests", fill=TRUE)
options( echo=TRUE)

