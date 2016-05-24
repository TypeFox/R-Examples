# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( LatticeKrig)
options( echo=FALSE)

##########################################
  test.for.zero.flag<- 1
  data( ozone2)
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]

  N<- length( y)
  a.wght<- 5
  lambda <-  1.5
  obj<- LKrig( x,y,NC=16, lambda=lambda, a.wght=a.wght, alpha=1, nlevel=1, NtrA=20,iseed=122)
  LKinfo<- obj$LKinfo
  K<- LKrig.cov( x,x,LKinfo)
  tempM<-  K
  diag(tempM) <- (lambda) + diag(tempM)
# Mi is proportional to the  inverse of the covariance matrix for the observations 
  Mi<- solve( tempM)
  T.matrix<- cbind( rep(1,N),x)
# this is estimating the fixed part using generalized least squares
  d.coef0 <-  solve( t(T.matrix)%*%Mi%*%T.matrix, t(T.matrix)%*%Mi%*%y)
  test.for.zero( obj$d.coef, d.coef0, tag="d from LKrig and by hand")
#### this is "c" coefficients for standard Kriging equations as done in mKrig
  temp2<- chol( tempM)
  c.coef0 <- forwardsolve(temp2, transpose = TRUE,
                        (y- T.matrix%*%d.coef0), upper.tri = TRUE)
  c.coef0 <- backsolve(temp2, c.coef0)
### find these using mKrig (still standard Kriging) but using the the LatticeKrig covariance function:
  obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
  test.for.zero( obj0$c, c.coef0, tag="c from mKrig and by hand" )
# we also know that for standard Kriging
# residuals = lambda* c.coef0
# use this to check the initial LatticeKrig result
 test.for.zero( obj0$fitted.values, obj$fitted.values)
# here is a nontrivial test: compaaring the same Kriging estimate
# using mKrig (via covariance function) and LatticeKrig (via precision and S-M-W identities)
 test.for.zero( lambda*obj0$c, (y-obj$fitted.values),
               tag="c from mKrig and from residuals of LatticeKrig (this is big!)" )
# compare Monte Carlo estimates of trace
 test.for.zero( obj$trA.info, obj0$trA.info, tag="Monte Carlo traces")
#
# test more complex covariance model:
#
  alpha<- c(1,.5,.2)
  nlevel<-3
  a.wght<-  c(5,5,10)
  lambda<- .1
  obj<- LKrig( x,y,NC=5, lambda=lambda,
                        nlevel=nlevel, alpha=alpha,a.wght=a.wght, NtrA=20,iseed=122)
  LKinfo<- obj$LKinfo

  obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
  test.for.zero( obj0$fitted.values, obj$fitted.values)
  test.for.zero( obj$d.coef, obj0$d, tag= "d from Lattice Krig and mKrig")

#### tests with  spatially varying alpha's

# first a sanity check that the marginalization and alpha option is working
# even when alpha's constant (obviously this is useful for debugging!)
# nu weighted alpha that constant by level to get the basis function counts.
   LKinfo<- LKrigSetup( x, NC=5, nlevel=3, nu=1, a.wght=5)
   N<- LKinfo$latticeInfo$mx[,1] * LKinfo$latticeInfo$mx[,2]
   glist<- fields.x.to.grid( x,10, 10)
   xg<-  make.surface.grid(glist)
   alpha.list<-  list( rep( 1,N[1]), rep(.5,N[2]), rep( .2,N[3]))
   LKinfo2<- LKrigSetup( x, NC=5, nlevel=3, alpha= alpha.list, a.wght=5)
   PHI1<- LKrig.basis(x, LKinfo2)
   PHI2<- LKrig.basis(xg, LKinfo2)                  
   Q<- LKrig.precision( LKinfo2)
# find covariane matrix "by hand" for this case.
   Ktest1<- PHI1%*%solve(Q)%*%t(PHI2)
   test.for.zero( Ktest1, LKrig.cov( x,xg, LKinfo=LKinfo2), tag="spatial alpha cov")
# check marginal variance 
   Ktest2<- PHI2%*%solve(Q)%*%t(PHI2)
   test.for.zero( diag(Ktest2), LKrig.cov( xg, LKinfo=LKinfo2, marginal =TRUE), tag="spatial alpha mar var")

# now varying alphas using same kind of tests
   set.seed( 678)
# here alpha weights are all different and random.
   alpha.list<- list(  runif(N[1]), runif(N[2]), runif(N[3]) )
 LKinfo2<- LKrigSetup( x, NC=5, nlevel=3, alpha= alpha.list, a.wght=5)
   PHI1<- LKrig.basis(x, LKinfo2)
   PHI2<- LKrig.basis(xg, LKinfo2)                  
   Q<- LKrig.precision( LKinfo2)
   Ktest1<- PHI1%*%solve(Q)%*%t(PHI2)
   test.for.zero( Ktest1, LKrig.cov( x,xg, LKinfo=LKinfo2), tag="spatial alpha cov w/ tricksy alpha")
   Ktest2<- PHI2%*%solve(Q)%*%t(PHI2)
   test.for.zero( diag(Ktest2), LKrig.cov( xg, LKinfo=LKinfo2, marginal =TRUE), tag="spatial alpha mar var  tricksy alpha")

  lambda<- .1
  obj<- LKrig( x,y,NC=5, lambda=lambda,
                        nlevel=nlevel, alpha=alpha,a.wght=a.wght, NtrA=20,iseed=122)
  LKinfo<- obj$LKinfo
  obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=LKinfo),
                                 NtrA=20, iseed=122)
#  check predicted values
  glist<- fields.x.to.grid( x,10, 10)
  xg<-  make.surface.grid(glist)
  grid.info<- obj$LKinfo$grid.info
  LKinfo<- obj$LKinfo
# first "by hand"
  Tmatrix<- cbind( rep(1,nrow(xg)), xg)
  yhat0<- Tmatrix%*%obj0$d +
            LKrig.cov( xg,x, LKinfo)%*%obj0$c
  PHIg<- LKrig.basis( xg, LKinfo)
  yhat1<- Tmatrix%*%obj$d.coef + PHIg%*%obj$c.coef
  test.for.zero( yhat1, predict(obj,xg), tag="predicted values from LatticeKrig and by hand, spatial alpha"  )
  test.for.zero( predict(obj,xg), predict(obj0,xg),
                           tag="predicted values LatticeKrig and mKrig, spatial alpha")
########## done with spatially varying alpha

# tests for computing the determinant and quad form from log likelihood
# see LatticeKrig tech report to sort these opaque computations!
  rm( obj, obj0) # remove previous objects
  test.for.zero.flag<- 1
  alpha<- c(1,.5,.5)
  nlevel<-3
  a.wght<-  c(5,5,10)
  lnDet<- function(A){
  sum( log( eigen( A, symmetric=TRUE)$values))}

  data( ozone2)
  x<-ozone2$lon.lat
  y<- ozone2$y[,16]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  x<- x[1:6,]
  y<- y[1:6]
#x<- transformx(x, "range")
  N<- length( y)
  lambda <- .8
# a micro sized lattice so determinant is not too big or small
  obj<- LKrig( x,y,NC=3, NC.buffer=1, lambda=lambda,nlevel=nlevel,alpha=alpha,a.wght=a.wght,
                              NtrA=5,iseed=122)
  LKinfo<- obj$LKinfo
  grid.info<- LKinfo$grid.info
  PHI<- LKrig.basis( x,LKinfo)
  Q <- LKrig.precision(LKinfo)
# coerce to full matrix
  Q<- as.matrix(Q)
  Mtest<- PHI%*% (solve( Q)) %*% t( PHI) + diag(lambda, N)
  temp<- t(PHI)%*%PHI + lambda*Q
  A<- Q*lambda
  B1<-  PHI%*% (solve( A)) %*% t( PHI) + diag(1, N)
  B2<-  t(PHI)%*%PHI + A
# the bullet proof application of identity 
  test.for.zero(lnDet( B1),lnDet( B2)- lnDet(A))
  test.for.zero(
                 lnDet( PHI%*% (solve( Q*lambda)) %*% t( PHI) + diag(1, N)),
                 lnDet( t(PHI)%*%PHI + Q*lambda) - lnDet(Q*lambda) )

# now adjusting for lambda factor 
  test.for.zero( lambda*B1, Mtest)
  test.for.zero(lnDet( Mtest), lnDet(B2) - lnDet(lambda*Q) + N*log(lambda) )
  test.for.zero(lnDet( Mtest), lnDet(B2) - lnDet(Q) + (-LKinfo$latticeInfo$m + N)*log(lambda) )

# find log determinant of temp using cholesky factors
# applying det identity
   temp<- t(PHI)%*%PHI + lambda*Q
   chol( temp)-> Mc

   lnDetReg <- 2 * sum(log(diag(Mc)))
   lnDetQ<-  2* sum( log( diag( chol(Q))))
   lnDetCov<- lnDetReg - lnDetQ + (-LKinfo$latticeInfo$m + N)*log(lambda)
   test.for.zero( lnDetCov, lnDet( Mtest))
   test.for.zero( obj$lnDetCov, lnDet( Mtest), tag="LatticeKrig and direct test of lnDetCov")
#
###### check of formula with weights
  set.seed(123)
  weights<- runif(N)
  W<- diag(weights)
  lambda<- .5
  PHI<- as.matrix(LKrig.basis( x,LKinfo))
  Q <- as.matrix(LKrig.precision(LKinfo))
  M1<- PHI%*%solve( Q)%*%t(PHI) +  lambda*solve( W) 

  
  B1<- (t(PHI)%*%(W/lambda)%*%PHI + Q)
  B2<- (1/lambda) * ( t(PHI)%*%(W)%*%PHI + lambda*Q)
  B3<-  t(PHI)%*%(W)%*%PHI + lambda*Q
  N2<- nrow(Q)
  hold<- lnDet( M1)
test.for.zero(   lnDet( B1) - lnDet(Q) - lnDet( W/lambda), hold, tag="Det with weights1")
test.for.zero(   lnDet( B2) - lnDet(Q) - lnDet( W/lambda), hold,  tag="Det with weights1=2")
test.for.zero(  lnDet( B3) - lnDet(Q) - sum( log( weights))  + (N-N2)*log(lambda), hold, tag="Det with weights3")

         
# now check these formulas as implemented in LatticeKrig
 rm( obj) # remove previous objects
 data( ozone2)
  x<-ozone2$lon.lat[1:10,]
  y<- ozone2$y[16,1:10]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
#x<- transformx(x, "range")
  N<- length( y)
  lambda <- .8
# a micro sized lattice so determinant is not too big or small
  obj<- LKrig( x,y,NC=5, lambda=lambda,nlevel=nlevel,alpha=alpha,a.wght=a.wght,
                              NtrA=5,iseed=122)
    obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=obj$LKinfo),
                                 NtrA=20, iseed=122)
 
 test.for.zero( obj$lnDetCov,obj0$lnDetCov, tag= "lnDetCov for mKrig and LatticeKrig")
 test.for.zero( obj$quad.form,  obj0$quad.form, tag= "quadratic forms for rho hat")
 test.for.zero(  obj0$lnProfileLike, obj$lnProfileLike,
                                tag="Profile Likelihood concentrated on lambda" )

# repeat tests for weighted measurement errors.
# recopy data to make reading easier
  rm( obj, obj0) # remove previous objects
  data( ozone2)
  x<-ozone2$lon.lat[1:10,]
  y<- ozone2$y[16,1:10]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
#x<- transformx(x, "range")
  N<- length( y)
  alpha<- c(1,.5,.25)
  nlevel<-3
  a.wght<-  c(5, 5, 4.5)
  lambda <- .5
  N<- length(y)
  set.seed(243)
  weights<- runif(N)*10 + 30
#  weights<- rep( 1, N)
  test.for.zero.flag<- 1
  obj<- LKrig( x,y,weights,NC=5,
                    lambda=lambda,alpha=alpha,nlevel=nlevel, a.wght=a.wght, NtrA=5,iseed=122)
# compare mKrig and Krig with weights and LatticeKrig
  obj0<- mKrig( x,y,weights, lambda=lambda, m=2, cov.function="LKrig.cov",
                                 cov.args=list(LKinfo=obj$LKinfo),
                                 NtrA=20, iseed=122)
 
  obj1<- Krig( x,y,weights=weights, lambda=lambda,GCV=TRUE, m=2,
               cov.function="LKrig.cov", cov.args=list(LKinfo=obj$LKinfo))
            
 test.for.zero( obj0$fitted.values, obj1$fitted.values)
 test.for.zero( predict(obj0), predict(obj1), tag="predicted  values mKrig/Krig  w/weights")
 test.for.zero( obj0$rhohat, obj1$rhohat,tag="compare rhohat for mKrig and Krig with weights")

############ now tests for LatticeKrig

 test.for.zero( obj$fitted.values, obj0$fitted.values)
 test.for.zero( obj$rho.MLE, obj0$rho.MLE)
 test.for.zero( obj$lnDetCov, obj0$lnDetCov)
############# tests using reuse Mc options
 data( ozone2)
  x<-ozone2$lon.lat[1:20,]
  y<- ozone2$y[16,1:20]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
  N<- length(y)
  set.seed(243)
  weights<- runif(N)*10
#x<- transformx(x, "range")
  N<- length( y)
  alpha<- c(1,.5,.5)
  nlevel<-3
  a.wght<-  c(5,5,10)
  lambda <- .8

 obj<- LKrig( x,y,weights=weights,NC=15, lambda=lambda,alpha=alpha,
                    nlevel=nlevel,a.wght=a.wght, return.cholesky=TRUE)
 obj2<- LKrig( x,y,weights=weights,NC=15, lambda=2*lambda,alpha=alpha,
                    nlevel=nlevel,a.wght=a.wght, use.cholesky=obj$Mc)
 obj3<-  LKrig( x,y,weights=weights,NC=15, lambda=2*lambda,alpha=alpha,
                    nlevel=nlevel,a.wght=a.wght, return.cholesky=TRUE)

 test.for.zero( obj3$c.coef, obj2$c.coef, tag="test of LatticeKrig.coef c")
 test.for.zero( obj3$d.coef, obj2$d.coef, tag="test of LatticeKrig.coef d")

 Q<- LKrig.precision(obj3$LKinfo)
 look2<-LKrig.lnPlike(obj3$Mc,Q,sqrt(weights)*y, obj3$residuals, weights,obj3$LKinfo)
 test.for.zero( look2$lnProfileLike, obj3$lnProfileLike)

# all done!
 cat("Done testing LatticeKrig",fill=TRUE)
 options( echo=FALSE)


# SUPPLEMENT: commented out sanity checks for  weighted/unweighted versions of mKrig and Krig
#hold0<-Krig ( x,y,weights=weights,method="user",GCV=TRUE,lambda=1e-3,
#             cov.function="Exp.simple.cov", cov.args=list( theta=300) )
#hold1<-mKrig(x,y,weights, lambda=1e-3,cov.function="Exp.simple.cov", cov.args=list( theta=300))
#test.for.zero( predict(hold0), predict(hold1))









