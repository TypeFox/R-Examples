# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html



# tests of predictSE
# against direct linear algebra 

library(fields)
options( echo=FALSE)

test.for.zero.flag<- TRUE

x0<- cbind( 0,4)

Krig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", theta=50,
      lambda=.06, GCV=FALSE)-> out

# direct calculation
Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%ChicagoO3$y, predict( out, x0),tag="Amatrix vs. predict")

Sigma0<- out$rhohat*Exp.cov( ChicagoO3$x, ChicagoO3$x, theta=50)
S0<- out$rhohat*c(Exp.cov( x0, x0, theta=50))
S1<- out$rhohat*Exp.cov( out$x, x0, theta=50)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

look<- S0 - t(S1)%*% t(A) - A%*%S1 +  
       A%*% ( Sigma0 + diag(out$shat.MLE**2/out$weightsM))%*% t(A)
#
#compare to 
# diagonal elements


test2<- predictSE( out, x= x0) 
test.for.zero( sqrt(diag(  look)), test2,tag="Marginal predictSE")


# now test shortcut formula that leverages the prediction step for Kriging
#

Sigma<-  Exp.cov( ChicagoO3$x, ChicagoO3$x, theta=50) +
          diag(out$lambda/out$weightsM)

#Sigma<-  ( Sigma0 + diag(out$shat.MLE**2/out$weightsM))

Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM)))

Omega<-  solve( t(Tmatrix)%*% solve( Sigma)%*% Tmatrix)
Id<- diag( 1, nrow( Tmatrix))

 qr.R( out$matrices$qr.VT) -> Rmat

Omega.test<- solve(t(Rmat)%*% (Rmat))

# Omega is the GLS covariance matrix for estimated parameters in fixed part of
# spatial model (the d coefficients). These are usually the "spatial drift" -- a
# low order polynomial

test.for.zero( Omega, Omega.test, tag="comparing Omega")

# M1 and M2 are matrices that go from obs to the estimated coefficients (d,c)

M1<- Omega%*% t(Tmatrix)%*% solve( Sigma)
M2<- solve( Sigma)%*% ( Id - Tmatrix%*% M1)


x0<- cbind( 0,4)

k0<-  Exp.cov( out$x, x0, theta=50)

#k0<- S1

t0<- c( 1, c(x0))

hold<- t( t0)%*%M1 + t(k0)%*% M2
test.for.zero( hold, A)
test.for.zero( M2%*%Sigma%*%t( M2), M2)

# benchmark using standard predictSE function

SE0<- predictSE( out, x=x0)

# shortcut formula explicitly
MSE<- S0  + out$rhohat*t(t0)%*% Omega %*%t0 -
            out$rhohat*(t(k0)%*% M2 %*% k0  + t(t0)%*% M1%*% k0) -
            out$rhohat*t(t0)%*% M1%*% k0

# collecting terms to make this look like  two predict steps. 
MSE2<- S0 + out$rhohat*t(t0)%*% Omega %*%t0 -
            out$rhohat* predict( out, yM= k0, x=x0) -
            out$rhohat* predict( out, yM= k0, x=x0, just.fixed=TRUE)
            
hold<- Krig.coef(out, y=k0)

tempc<-  t(k0)%*% hold$c 
tempd<-  t(t0)%*%hold$d 

MSE4<- S0 + out$rhohat*t(t0)%*% Omega %*%t0 -
                 out$rhohat * (tempc +2*tempd)

test.for.zero(SE0, sqrt( MSE4), tag="test of formula with explicit d and c")


# test of new function

Krig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", theta=50,lambda=.06)-> out0
SE0<- predictSE.Krig( out0, x=x0)
mKrig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", theta=50, lambda=.06)-> out2
SE3<- predictSE.mKrig( out2, xnew=x0)

test.for.zero(SE0, sqrt( MSE), tag="Krig function and direct formula")


test.for.zero(sqrt(MSE), sqrt( MSE2), 
             tag="new predict formula and direct formula")

test.for.zero( SE3, SE0,  tag="New se _function_ and old Krig _function_")
#
# test of vectors of locations.


# receate object
Krig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", theta=50, lambda=.06)-> out0
mKrig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", theta=50, lambda=.06)-> out


x0<-rep( c( -20, -10,10,20),4)

x0 <- cbind( x0 , sort( x0))
x0<- rbind( c(0,4), x0)

k0<-  Exp.cov(  ChicagoO3$x,x0, theta=50)
t0<- t(fields.mkpoly(x0, m=out$m))
hold<- Krig.coef(out0, y=k0)

MSE5<- (rep( S0,nrow(x0)) +
                  out0$rhohat * colSums( t0 *(out0$matrices$Omega%*%t0)) 
                  -out0$rhohat* colSums((k0)*hold$c) - 
                   2*out0$rhohat*colSums(t0*hold$d))

hold<- mKrig.coef(out, y=k0)
MSE6<- (rep( S0,nrow(x0)) +
                  out$rhohat * colSums( t0 *(out$Omega%*%t0))
                  -out$rhohat* colSums((k0)*hold$c) -
                   2*out$rhohat*colSums(t0*hold$d))

test.for.zero( predictSE( out0, x0), sqrt(MSE5), 
                   tag="Benchmark of formula")

test.for.zero( predictSE( out0, x0), sqrt(MSE6), 
                   tag="Benchmark of formula mKrig coefs")

test.for.zero( predictSE( out, x0), predictSE.mKrig(out, x0),
                   tag="test function with several locations Krig mKrig functions" )


