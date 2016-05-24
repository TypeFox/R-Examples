# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( fields)
options( echo=FALSE)
test.for.zero.flag<- 1
data(ozone2)
y<- ozone2$y[16,]
x<- ozone2$lon.lat
#
# Omit the NAs
good<- !is.na( y)
x<- x[good,]
y<- y[good]
#source("~/Home/Src/fields/R/mKrig.family.R")

# now look at mKrig w/o sparse matrix 
mKrig( x,y, cov.function="stationary.cov", theta=10, lambda=.3,
                      chol.args=list( pivot=FALSE))-> look


Krig( x,y, cov.function="stationary.cov",
               theta=10) -> look2

test.df<-Krig.ftrace(look$lambda,look2$matrices$D)

test<- Krig.coef( look2, lambda=look$lambda)

test.for.zero( look$d, test$d, tag="Krig mKrig d coef")
test.for.zero( look$c, test$c, tag="Krig mKrig c coef")

# test of trace calculation

mKrig( x,y, cov.function="stationary.cov", theta=10, lambda=.3,
         
          find.trA=TRUE, NtrA= 1000, iseed=243)-> look

test.for.zero( look$eff.df, test.df,tol=.01, tag="Monte Carlo eff.df")


# 
Krig( x,y, cov.function="stationary.cov",
               theta=350, Distance="rdist.earth",Covariance="Wendland", 
               cov.args=list( k=2, dimension=2) ) -> look2

mKrig( x,y, cov.function="stationary.cov", 
        theta=350, 
        Distance="rdist.earth",Covariance="Wendland",  
        cov.args=list( k=2, dimension=2),
        lambda=look2$lambda,
        find.trA=TRUE, NtrA= 1000, iseed=243)-> look

test.for.zero( look$c, look2$c, tag="Test of wendland and great circle")

test.for.zero(look$eff.df, Krig.ftrace( look2$lambda, look2$matrices$D)
              ,tol=.01, tag="eff.df")

# same calculation using sparse matrices.

mKrig( x,y, cov.function="wendland.cov", 
        theta=350, 
        Dist.args=list( method="greatcircle"),  
        cov.args=list( k=2),
        lambda=look2$lambda,
        find.trA=TRUE, NtrA=100, iseed=243)-> look

test.for.zero( look$c, look2$c,tol=8e-7, 
           tag="Test of sparse wendland and great circle")
test.for.zero(look$eff.df, Krig.ftrace( look2$lambda, look2$matrices$D),
                        tol=.01, tag="sparse eff.df")

# great circle distance switch has been a  big bug -- test some options

mKrig( x,y, cov.function="wendland.cov", 
 theta=350, Dist.args=list( method="greatcircle"),  
 cov.args=list( k=2),lambda=look2$lambda,
 find.trA=TRUE, NtrA=200, iseed=243)-> look

test.for.zero(look$eff.df, Krig.ftrace( look2$lambda, look2$matrices$D),
                   tol=1e-6, tag="exact sparse eff.df")

# compare to fast Tps 
 fastTps( x,y,theta=350,lambda=look2$lambda, NtrA=200, iseed=243, 
                lon.lat=TRUE)-> look3

test.for.zero(look3$eff.df, Krig.ftrace( look2$lambda, look2$matrices$D),
                   tol=1e-6, tag="exact sparse eff.df -- fastTps")

# calculations of likelihood, rho and sigma

lam<-.2

out<- mKrig( x,y, cov.function =Exp.cov, theta=4, lambda=lam)
out2<- Krig( x,y, cov.function =Exp.cov, theta=4, lambda=lam)

            
Sigma<- Exp.cov( x,x,theta=4)
X<-  cbind( rep(1, nrow(x)), x)

Sinv<- solve( Sigma + lam* diag( 1, nrow( x)))

#checks on  likelihoods            

# quadratic form:
dhat<- c(solve( t(X)%*%Sinv%*%(X) ) %*% t(X) %*%Sinv%*%y)
test.for.zero( dhat, out$d, tag="initial check on d for likelihood")
r<- y -X%*%dhat
N<- nrow(x)
look<-  t( r)%*%(Sinv)%*%r/N



test.for.zero( look, out$rho.MLE, tag="rho hat from likelihood")

test.for.zero( look, out2$rhohat, tag="rho hat from likelihood compared to Krig")



# check determinant
lam<- .2
Sigma<- Exp.cov( x,x,theta=4)
M<- Sigma + lam * diag( 1, nrow(x))
chol( M)-> Mc
look2<- sum( log(diag( Mc)))*2

out<-mKrig( x,y,cov.function =Exp.cov, theta=4, lambda=lam)

test.for.zero( out$lnDetCov, look2)
test.for.zero( out$lnDetCov, determinant(M, log=TRUE)$modulus)

# weighted version 
lam<- .2
Sigma<- Exp.cov( x,x,theta=4)
set.seed( 123)
weights<- runif(nrow( x))
M<- Sigma +  diag(lam/ weights)
chol( M)-> Mc
look2<- sum( log(diag( Mc)))*2

out<-mKrig( x,y,weights=weights, cov.function =Exp.cov, theta=4, lambda=lam)

test.for.zero( out$lnDetCov, look2)
test.for.zero(  look2, determinant(M, log=TRUE)$modulus)
test.for.zero( out$lnDetCov, determinant(M, log=TRUE)$modulus)



# check profile likelihood by estimating MLE
lam.true<- .2
N<- nrow( x)
Sigma<- Exp.cov( x,x,theta=4)
M<- Sigma + lam.true * diag( 1, nrow(x))
chol( M)-> Mc
t(Mc)%*%Mc -> test




##D set.seed( 234)
##D NSIM<- 100
##D hold2<-rep( NA, NSIM)
##D temp.fun<- function(lglam){
##D             out<-mKrig( x,ytemp,
##D                         cov.function =Exp.cov, theta=4, lambda=exp(lglam))
##D             return(-1* out$lnProfileLike)}

##D hold1<-rep( NA, NSIM)
##D yt<- rep( 1, N) 
##D obj<- Krig( x,yt, theta=4)


##D E<- matrix( rnorm( NSIM*N), ncol=NSIM)

##D for ( j in 1:NSIM){
##D cat( j, " ")
##D ytemp <- x%*%c(1,2) +  t(Mc)%*%E[,j] 
##D out<- optim( log(.2), temp.fun, method="BFGS")
##D hold2[j]<- exp(out$par)
##D hold1[j]<-  gcv.Krig(obj, y=ytemp)$lambda.est[6,1]

##D }
##D test.for.zero( median( hold1), .2, tol=.08)
##D test.for.zero( median( hold2), .2, tol=.12)



            










