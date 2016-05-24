# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( fields)
options( echo=FALSE)
test.for.zero.flag<- 1

# test data
data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]

#first test addToDiagC

I3 = diag(nrow=3)
twoI3 = I3*2
.Call("addToDiagC", I3, rep(1.0, 3), as.integer(3))
test.for.zero(twoI3, I3, tag="addToDiag")

# turning spam on and off
Krig(x,y, cov.function = "stationary.taper.cov", theta=1.5,
     cov.args= list( spam.format=FALSE,
                     Taper.args= list( theta=2.0,k=2, dimension=2) )
) -> out1

Krig(x,y, cov.function = "stationary.taper.cov", lambda=2.0, theta=1.5,
     cov.args= list( spam.format=TRUE,
                     Taper.args= list( theta=2.0,k=2, dimension=2) )
) -> out2

temp1<- predict( out1,lambda=2.0)
temp2<- predict( out2)
test.for.zero( temp1, temp2, tag="spam vs no spam")

#
# Omit the NAs
good<- !is.na( y)
x<- x[good,]
y<- y[good]

# now look at mKrig w/o sparse matrix 
mKrig( x,y, cov.function="stationary.cov", theta=10, lambda=.3,
       chol.args=list( pivot=FALSE))-> look

Krig( x,y, cov.function="stationary.cov", theta=10, lambda=.3) -> look2

test.for.zero( look$d, look2$d, tag="Krig mKrig d coef")
test.for.zero( look$c, look2$c, tag="Krig mKrig c coef")


set.seed(123)
xnew<- cbind( (runif(20)-.5)*5, (runif(20)-.5)*5)
temp<- predict( look, xnew)
temp2<- predict( look2, xnew)
test.for.zero( temp, temp2, tag="test of predict at new locations")

# test of matrix of obs
N<- length( y)
Y<- cbind( runif(N), y,runif(N), y)

mKrig( x,Y, cov.function="stationary.cov", 
       theta=10, lambda=.3)-> lookY
temp3<-  predict( lookY, xnew)[,4]

test.for.zero( temp, temp3, tag="test of matrix Y predicts" )
predictSurface( look)-> temp
predictSurface( look2)-> temp2

good<- !is.na( temp2$z)
test.for.zero( temp$z[good], temp2$z[good])

# testing stationary taper covariance 
# and also surface prediction

N<- length( y)
mKrig( x,y, cov.function="stationary.taper.cov", theta=2, 
       spam.format=FALSE, lambda=.35 )-> look

Krig( x,y, cov.function="stationary.taper.cov", theta=2, 
      spam.format=FALSE, lambda=.35)-> look2

predictSurface( look, nx=50, ny=45)-> temp
predictSurface( look2, nx=50, ny=45)-> temp2

good<- !is.na( temp2$z)
test.for.zero( temp$z[good], temp2$z[good], tag="predictSurface with mKrig")

# 
# Use Wendland with sparse off and on.
Krig( x,y, cov.function="wendland.cov", 
      cov.args=list( k=2, theta=2.8), 
      lambda=.3, spam.format=FALSE)-> look

mKrig( x,y, cov.function="wendland.cov",k=2, theta=2.8,
       spam.format=FALSE, lambda=.3)-> look2

mKrig( x,y, cov.function="wendland.cov",k=2, theta=2.8,
       spam.format=TRUE, lambda=.3)-> look3

# final tests for  predict.
set.seed(223)
xnew<- cbind(runif( N)*.5 + x[,1], runif(N)*.5 + x[,2])
temp<- predict( look, xnew)
temp2<- predict( look2, xnew)
temp3<- predict( look3, xnew)
test.for.zero( temp, temp2, tag="Wendland/no spam")
test.for.zero( temp2, temp3, tag="Wendland/spam")


### testing coefficients for new data 
mKrig.coef( look2, cbind(y+1,y+2))-> newc
test.for.zero( look2$c, newc$c[,2], tag="new coef c no spam")

test.for.zero( look2$d,
               c(newc$d[1,2] -2, newc$d[2:3,2]), tag="new d coef no spam")

mKrig.coef( look3, cbind(y+1,y+2))-> newc
test.for.zero( look3$c, newc$c[,2], tag="new coef c spam")

test.for.zero( look3$d,
               c(newc$d[1,2] -2, newc$d[2:3,2]), tag="new d coef spam")

###


### bigger sample size
set.seed( 334)
N<- 1000
x<- matrix( runif(2*N),ncol=2)
y<- rnorm( N)
nzero <- length( wendland.cov(x,x, k=2,theta=.1)@entries)


mKrig( x,y, cov.function="wendland.cov",k=2,
       theta=.1, lambda=.3)-> look2


test.for.zero( look2$non.zero.entires, nzero, tag="nzero in call to mKrig")

###### 
### test out passing to chol

data( ozone2)
y<- ozone2$y[16,]
good<- !is.na( y)
y<-y[good]
x<- ozone2$lon.lat[good,]

# interpolate using defaults (Exponential)
# stationary covariance
mKrig( x,y, theta = 1.5, lambda=.2)-> out
#
# NOTE this should be identical to 
Krig( x,y, theta=1.5, lambda=.2) -> out2
temp<- predict( out)
temp2<- predict( out2)

test.for.zero( temp, temp2, tag="mKrig vs. Krig for ozone2")

# test passing arguments for chol 

set.seed( 334)
N<- 300
x<- matrix( 2*(runif(2*N)-.5),ncol=2)
y<- sin( 3*pi*x[,1])*sin( 3.5*pi*x[,2]) + rnorm( N)*.01


Krig( x,y, Covariance="Wendland",
      cov.args= list(k=2, theta=.8, dimension=2),                   , 
      give.warnings=FALSE,
      lambda=1e2) -> out

mKrig( x,y, 
       cov.function="wendland.cov",k=2, theta=.8, 
       lambda=1e2, 
       chol.args=list( memory=list( nnzR=1e5)), 
)-> out2

temp<- predict( out)
temp2<- predict( out2)

test.for.zero( temp, temp2, tag="predict Wendland  mKrig vs Krig")




# test of fastTps
nx<- 50
ny<- 60
x<- seq( 0,1,,nx)
y<- seq( 0,1,,ny)
gl<- list( x=x, y=y)
xg<- make.surface.grid(gl)
ztrue<- sin( xg[,1]*pi*3)* cos(xg[,2]*pi*2.5)
#image.plot(x,y,matriz( ztrue, nx,ny)) 
set.seed( 222)
ind<- sample( 1:(nx*ny), 600)
xdat<- xg[ind,]
ydat <- ztrue[ind]
out<- fastTps(xdat, ydat, theta=.3)
out.p<-predictSurface( out, grid=gl, extrap=TRUE)
# perfect agreement at data
test.for.zero( ydat, c( out.p$z)[ind], tag="fastTps interp1")
#image.plot(x,y,matrix( ztrue, nx,ny)- out.p$z) 
rmse<- sqrt(mean( (ztrue- c( out.p$z))^2)/ mean( (ztrue)^2))
test.for.zero( rmse,0,tol=.01, relative=FALSE,tag="fastTps interp2")


##### test precomputing distance matrices:
#

set.seed(1)

# test data
data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]

#
# Omit the NAs
good<- !is.na( y)
x<- x[good,]
y<- y[good]
compactDistMat = rdist(x, compact=TRUE)
distMat = rdist(x)

##### test using distance matrix
print("testing using distance matrix")

mKrig(x,y, cov.function = "stationary.cov", lambda=2.0, theta=1.5) -> out1

mKrig(x,y, cov.args= list(Covariance="Exponential", Distance="rdist", Dist.args=list(compact=TRUE)), 
      lambda=2.0, theta=1.5) -> out2

#NOTE: compact distance matrix should not be used by user for fields compatibility reasons
mKrig(x,y, cov.args= list(Covariance="Exponential", Dist.args=list(compact=TRUE)), 
      lambda=2.0, theta=1.5, distMat=compactDistMat) -> out3

mKrig(x,y, cov.args= list(Covariance="Exponential"), 
      lambda=2.0, theta=1.5, distMat=distMat) -> out4

temp1<- predict( out1)
temp2<- predict( out2)
temp3 = predict( out3)
temp4 = predict( out4)
test.for.zero( temp1, temp2, tag="predict: stationary.cov versus Exp.cov")
test.for.zero( temp2, temp3, tag="predict: no distance matrix versus compact distance matrix")
test.for.zero( temp2, temp4, tag="predict: no distance matrix versus distance matrix")

##### test SE
print("testing using predictSE")

temp1 = predictSE(out1)
temp2 = predictSE(out2)
temp3 = predictSE(out3)
temp4 = predictSE(out4)

test.for.zero( temp1, temp2, tag="predictSE: stationary.cov with exponential versus Exp.cov")
test.for.zero( temp2, temp3, tag="predictSE: no distance matrix versus compact distance matrix")
test.for.zero( temp2, temp4, tag="predictSE: no distance matrix versus distance matrix")





cat("all done with mKrig tests", fill=TRUE)
options( echo=TRUE)



