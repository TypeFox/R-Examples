# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( fields)
# tests of predict.se
# using approximations for conditional simulation on a grid. 
#
options( echo=FALSE)
test.for.zero.flag<-1

long.test.flag<- FALSE

data( ozone2)
as.image(ozone2$y[16,], x= ozone2$lon.lat, ny=24, nx=20, 
          na.rm=TRUE)-> dtemp
#
# A useful discretized version of ozone2 data
 
x<- cbind(dtemp$x[dtemp$ind[,1]], dtemp$y[dtemp$ind[,2]])
y<- dtemp$z[ dtemp$ind]
weights<- dtemp$weights[ dtemp$ind]

Krig( x, y, Covariance="Matern", 
   theta=1.0, smoothness=1.0, weights=weights) -> out


if(long.test.flag){

# the grid ...

glist<- list( x= dtemp$x, y=dtemp$y)

set.seed( 233)
sim.Krig.approx( out, grid= glist, M=200, extrap=TRUE)-> look

predict.surface.se( out, grid=glist, extrap=TRUE)-> test

look2<- matrix( NA, 20,24)

for(  k in 1:24){
for ( j in 1:20){
look2[j,k] <- sqrt(var( look$z[j,k,], na.rm=TRUE))
}
}


test.for.zero(  mean( abs(look2- test$z)/test$z), 0, relative=FALSE,
tol=.05, tag="Conditional simulation marginal se for grid")

#
# test for covariances


ind0<- expand.grid( c(1,4,5,10), c(3,4,5,10, 15))

x0<- cbind( glist$x[ind0[,1]], glist$y[ind0[,2]])
look2<- matrix( NA, 200,20)
for(  k in 1:20){
look2[,k] <- look$z[ ind0[k,1], ind0[k,2],]}

predict.se( out, x0, cov=TRUE)-> test2
ds<- 1/sqrt(diag(test2))
test3<-  diag(ds)%*% test2 %*% diag(ds)

#check plot( diag( test2), diag( var( look2)))

# Another plot to look at plot( c(test3), c(cor(look2)))

hold<-cor(look2)
upper<- col(hold)> row( hold)
dd<- (c(hold)- c(test3))[upper]

test.for.zero(   mean( abs(dd)) ,0, relative=FALSE,
tol=.05, tag="Conditional simulation correlations for grid (RMSE) ")

} # end long test block

cat( "all done with grid based se tests", fill=TRUE)
options( echo=TRUE)
