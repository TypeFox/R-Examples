# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( LatticeKrig)
options( echo=FALSE)


# tests for computing the determinant and quad form
  test.for.zero.flag<- 1
  alpha<- c(1,.5,.5)
  nlevel<-3
  a.wght<-  c(5,5,10)
  lnDet<- function(A){
  sum( log( eigen( A, symmetric=TRUE)$values))}

  data( ozone2)

  x<-ozone2$lon.lat[1:20,]
  y<- ozone2$y[16,1:20]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
#x<- transformx(x, "range")
  N<- length( y)
  lambda <- .8
 set.seed(123)
  weights<- runif(N)
  W<- diag(weights)
# a micro sized lattice so determinant is not too big or small
  obj<- LKrig( x,y,NC=5, weights= weights, lambda=lambda,nlevel=nlevel,alpha=alpha,a.wght=a.wght,
                              NtrA=5,iseed=122)
  LKinfo<- obj$LKinfo
# now check these formulas as implemented in LatticeKrig
    obj0<- mKrig( x,y, weights= weights, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=5, iseed=122)
 
###### check of formula with weights
  PHI<- as.matrix(LKrig.basis( x,LKinfo))
  Q <- as.matrix(LKrig.precision(LKinfo))
  M1<- PHI%*%solve( Q)%*%t(PHI) +  lambda*solve( W) 
  B1<- (t(PHI)%*%(W/lambda)%*%PHI + Q)
  B2<- (1/lambda) * ( t(PHI)%*%(W)%*%PHI + lambda*Q)
  B3<-  t(PHI)%*%(W)%*%PHI + lambda*Q
  N2<- nrow(Q)
#  lnDet( M1)
#  lnDet( B1) - lnDet(Q) - lnDet( W/lambda)
#  lnDet( B2) - lnDet(Q) - lnDet( W/lambda)
  test.for.zero( lnDet( B3) - lnDet(Q) - sum( log( weights))  + (N-N2)*log(lambda),
                         lnDet( M1), tag="Direct formula")
  test.for.zero( obj$lnDetCov,  obj0$lnDetCov, tag= "lnDetCov for mKrig and LatticeKrig")
  test.for.zero( obj$quad.form,  obj0$quad.form, tag= "quadratic forms for rho MLE")
  test.for.zero( obj0$lnProfileLike, obj$lnProfileLike,
                                tag="Profile Likelihood concentrated on lambda" )
  Y<- cbind( y,2*y, y*3)
  obj<- LKrig( x,Y,NC=5, weights= weights, lambda=lambda,nlevel=nlevel,alpha=alpha,a.wght=a.wght,
                              NtrA=5,iseed=122)
  LKinfo<- obj$LKinfo
# now check these formulas as implemented in LatticeKrig
  obj0<- mKrig( x,Y, weights= weights, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
  test.for.zero(  obj0$rho.MLE, obj$rho.MLE,
                                tag="MLEs for rho with replicate fields" )

  test.for.zero(  obj0$lnProfileLike, obj$lnProfileLike,
                                tag="Profile Likelihood concentrated on lambda with replicate fields" )
# test of full likelihood
#    t(y- T%*%d.coef) %*% solve( sigma^2/w + rho*K)%*% y- T%*%d.coef)
#  this is the quad form found from shortcut formula
#
   sigma<- .1
   rho<- 2.2
   lambda<- sigma^2/rho
   obj.test<- LKrig( x,Y,NC=5, weights= weights, nlevel=nlevel,alpha=alpha,a.wght=a.wght,
                              NtrA=5,iseed=122, sigma=sigma, rho=rho)
   Tmatrix <- cbind(fields.mkpoly(x, 2))
   res<- Y - Tmatrix%*%obj.test$d.coef
   M<- (diag(lambda/weights) + LKrig.cov(x,x,LKinfo=obj.test$LKinfo))
   qtest<- diag( t(res)%*% solve( M)%*% res)
   test.for.zero( qtest, obj.test$quad.form, tag="quad form for arbitrary sigma and rho")

   n<- nrow( Y)
   lnDetCov.test<- sum( log( eigen(M)$values))
   test.for.zero( obj.test$lnDetCov, lnDetCov.test, tag="lnDetCov with arbitrary rho sigma")

   lnLike.test<-   (-qtest/(2*rho) - log(2*pi)*(n/2)
                       -(n/2)*log(rho) - (1/2) * lnDetCov.test) 
   test.for.zero( lnLike.test, obj.test$lnLike, tag="lnLike arbitrary sigma and rho")



cat("all done with lnPLike with weights", fill=TRUE)
options( echo=FALSE)
