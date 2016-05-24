# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


library( fields)

# tests of predictSE
# against direct linear algebra 

options( echo=FALSE)

test.for.zero.flag<- 1

x0<- expand.grid( c(-8,-4,0,20,30), c(10,8,4,0))


Krig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", theta=50)-> out


# direct calculation
Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%ChicagoO3$y, predict( out, x0),tag="Amatrix vs. predict")

Sigma<- out$rhohat*Exp.cov( ChicagoO3$x, ChicagoO3$x, theta=50)
S0<- out$rhohat*c(Exp.cov( x0, x0, theta=50))
S1<- out$rhohat*Exp.cov( out$x, x0, theta=50)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

look<- S0 - t(S1)%*% t(A) - A%*%S1 +  
       A%*% ( Sigma + diag(out$shat.MLE**2/out$weightsM))%*% t(A)
#
#compare to 
# diagonal elements


test2<- predictSE( out, x= x0) 
test.for.zero( sqrt(diag(  look)), test2,tag="Marginal predictSE")

out2<- Krig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", theta=50,
            lambda=out$lambda)

test2<- predictSE( out2, x= x0) 
test.for.zero( sqrt(diag(  look)), test2,tag="Marginal predictSE fixed ")

test<- predictSE( out, x= x0, cov=TRUE)
test.for.zero( look, test,tag="Full covariance predictSE")


# simulation based.

set.seed( 333)

sim.Krig( out, x0,M=4e3)-> test

var(test)-> look

predictSE( out, x=x0)-> test2
mean( diag( look)/ test2**2)-> look2
test.for.zero(look2, 1.0, tol=1.5e-2, tag="Marginal standard Cond. Sim.")

predictSE( out, x=x0, cov=TRUE)-> test2

# multiply simulated values by inverse square root of covariance
# to make them white

eigen( test2, symmetric=TRUE)-> hold
hold$vectors%*% diag( 1/sqrt( hold$values))%*% t( hold$vectors)-> hold
cor(test%*% hold)-> hold2
# off diagonal elements of correlations -- expected values are zero. 

abs(hold2[ col(hold2)> row( hold2)])-> hold3

test.for.zero(   mean(hold3), 0, relative=FALSE, tol=.02,
          tag="Full covariance standard Cond. Sim.")


# test of sim.Krig.approx.R
#
# first create and check a gridded test case. 


data( ozone2)
as.image(ozone2$y[16,], x= ozone2$lon.lat, ny=24, nx=20, 
          na.rm=TRUE)-> dtemp
#
# A useful disctrtized version of ozone2 data
 
x<- dtemp$xd
y<- dtemp$z[ dtemp$ind]
weights<- dtemp$weights[ dtemp$ind]

Krig( x, y, Covariance="Matern", 
   theta=1.0, smoothness=1.0, weights=weights) -> out



  set.seed(234)
  ind0<- cbind( sample( 1:20, 5), sample( 1:24, 5))

  x0<- cbind( dtemp$x[ind0[,1]], dtemp$y[ind0[,2]]) 

# an  inline check plot(out$x, cex=2); points( x0, col="red", pch="+",cex=2)

# direct calculation as backup ( also checks weighted case)

Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%out$yM, predict( out, x0),tag="Amatrix vs. predict")

Sigma<- out$rhohat*stationary.cov( 
out$xM, out$xM, theta=1.0,smoothness=1.0, Covariance="Matern")

S0<- out$rhohat*stationary.cov( 
x0, x0, theta=1.0,smoothness=1.0, Covariance="Matern")

S1<- out$rhohat*stationary.cov(
out$xM, x0, theta=1.0,smoothness=1.0, Covariance="Matern")



#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)
 
look<- S0 - t(S1)%*% t(A) - A%*%S1 +
       A%*% ( Sigma + diag(out$shat.MLE**2/out$weightsM) )%*% t(A)

test<- predictSE( out, x0, cov=TRUE)

test.for.zero( c( look), c( test), tag="Weighted case and exact for ozone2 full 
cov", tol=1e-8)

########################################################################
######### redo test with smaller grid to speed things up
#cat("Conditional simulation test -- this takes some time", fill=TRUE)

# redo data set to smaller grid size
##D N1<-4
##D N2<-5
##D as.image(ozone2$y[16,], x= ozone2$lon.lat, ny=N2, nx=N1, 
##D          na.rm=TRUE)-> dtemp
#
# A useful discretized version of ozone2 data
 
##D xd<- dtemp$xd
##D y<- dtemp$z[ dtemp$ind]
##D weights<- dtemp$weights[ dtemp$ind]

##D Krig( xd, y, Covariance="Matern", 
##D    theta=1.0, smoothness=1.0, weights=weights) -> out


##D xr<- range( dtemp$x)
##D yr<- range( dtemp$y)
##D M1<-N1
##D M2<- N2
##D glist<- list( x=seq( xr[1], xr[2],,M1) , y=seq( yr[1], yr[2],,M2))

##D set.seed( 233)
# with extrap TRUE this finesses problems with
# how NAs are handled in var below

##D sim.Krig.approx( out, grid= glist, M=3000, extrap=TRUE)-> look

##D predictSE( out, make.surface.grid( glist))-> test


##D look2<- matrix( NA, M1,M2)

##D for(  k in 1:M2){
##D     for ( j in 1:M1){
##D       look2[j,k] <- sqrt(var( look$z[j,k,], na.rm=TRUE)) }
##D }


##D test.for.zero(  1-mean(c(look2/test), na.rm=TRUE), 0, relative=FALSE, 
##D tol=.001, tag="Conditional simulation marginal se for grid")

cat("all done testing predictSE ", fill=TRUE)
options( echo=TRUE)
