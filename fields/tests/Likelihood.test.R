# this is a test script to verify the likelihood computations are 
# correct with the eigen decomposition format used in Krig
# see Krig.flplike for the concise computation.
#  

library(fields)

#options( echo=FALSE)
test.for.zero.flag<- 1


data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]
is.good <- !is.na( y)
x<- x[is.good,]
y<- y[is.good]

theta<- 2.0

# check log likelihood calculation
  nu<- 1.5
  lambda<- .2
  out<- mKrig( x,y, theta=theta,Covariance="Matern", smoothness=nu, lambda=lambda) 

# peg rho and sigma as MLEs from mKrig
  rho <- out$rho.MLE
  sigma2<- rho*lambda
  N<- length( y)
  dd<- rdist( x,x)
  M<-  rho* Matern( dd, range= theta, smoothness=nu) + sigma2* diag( 1, N)
  X<- fields.mkpoly( x, 2)
  Mi<- solve( M)
  betahat<-  solve(t(X)%*%Mi%*%X)%*% t(X)%*% Mi%*% y
  res<- y - X%*%betahat
  ccoef<- ( Mi%*% ( res))*rho

# sanity check that estimates are the same
  test.for.zero( ccoef, out$c, tag="check ccoef")

# find full log likelihood
  chol(M)-> cM
  lLike<-  -(N/2)*log(2*pi) - (1/2)* (2*sum( log( diag(cM)))) - (1/2)* t(res)%*% Mi %*% res
# formula for full likelihood using peices from mKrig
  lLike.test<- -(N/2)*log(2*pi) - (1/2)* out$lnDetCov - (1/2)*(N)*log( rho) - (1/2)*out$quad.form/rho

  test.for.zero( lLike, lLike.test, tag="llike full verses rhohat")
  test.for.zero( lLike, out$lnProfileLike, tag="llike profile from mKrig")

# REML check
  nu<- 1.5
  theta<- .6
  obj<- Krig( x,y, theta=theta,Covariance="Matern", smoothness=nu )

# sanity check that c coefficients agree with Krig
  rho<- 500
  lambda<- .2
  sigma2<- lambda*rho

  hold<- REML.test( x,y,rho, sigma2, theta, nu=1.5)
  ccoef2<- Krig.coef( obj, lambda)$c
  test.for.zero( hold$ccoef, ccoef2, tag="ccoefs")

# check RSS with Krig decomposition.
  RSS1<- sum( (lambda*ccoef2)**2)
  lD <- obj$matrices$D * lambda
  RSS2 <- sum(((obj$matrices$u * lD)/(1 + lD))^2) 
  test.for.zero( RSS2, RSS1, tag=" RSS using matrices")

# check quadratic form with Krig
  D.temp<- obj$matrices$D[  obj$matrices$D>0]
  A3test<- (1/lambda)* obj$matrices$V %*% diag((D.temp*lambda)/ (1 +D.temp*lambda) )%*% t( obj$matrices$V)
  test.for.zero(solve(A3test), hold$A/rho, tol=5e-8)
  Quad3<-   sum( D.temp*(obj$matrices$u[obj$matrices$D>0])^2/(1+lambda*D.temp))

  test.for.zero( hold$quad.form, Quad3/rho, tag="quad form")

# test determinants
  N2<- length( D.temp)
  det4<- -sum( log(D.temp/(1 + D.temp*lambda)) )
  det1<- sum( log(eigen( hold$A/rho)$values)) 
  test.for.zero( det1, det4, tag="det" )

# test REML Likelihood
  lLikeREML.test<--1*( (N2/2)*log(2*pi) - (1/2)*(sum( log(D.temp/(1 + D.temp*lambda)) ) - N2*log(rho)) +
                                      (1/2)*sum( lD*(obj$matrices$u)^2/(1+lD)) /(lambda*rho) )

test.for.zero( hold$REML.like, lLikeREML.test, tag="REML using matrices")


# profile likelihood

# lnProfileLike <- (-np/2 - log(2*pi)*(np/2)
#                      - (np/2)*log(rho.MLE) - (1/2) * lnDetCov)
#  test using full REML likelihood.
  nu<- 1.5
  rho<- 7000
  lambda<- .02
  sigma2<- lambda*rho
  theta<- 2.0
  obj<- Krig( x,y, theta=theta,Covariance="Matern", smoothness=nu )
  hold<- REML.test(x,y,rho, sigma2, theta, nu=1.5)
  np<- hold$N2
  rho.MLE<- c(hold$rhohat)
  lnDetCov<-sum( log(eigen( hold$A/rho)$values))  

  l0<- REML.test(x,y,rho.MLE, rho.MLE*lambda, theta, nu=1.5)$REML.like
  l1<-   (-np/2 - log(2*pi)*(np/2)- (np/2)*log(rho.MLE) - (1/2) * lnDetCov)
  l2<-  (-1)*Krig.flplike( lambda, obj)

  test.for.zero( l0,l2, tag="REML profile flplike")
  test.for.zero( l1,l2, tag="REML profile flplike")


cat("all done with likelihood  tests", fill=TRUE)
options( echo=TRUE)
      
