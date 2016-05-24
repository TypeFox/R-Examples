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
x1<- x[1:20,]
x2<- x[1:10,]

look<- exp(-1*rdist(x1,x2)/4)
look2<- stationary.cov( x1,x2, theta=4)
test.for.zero( look, look2)

V<- matrix( c(2,1,0,4), 2,2)
Vi<- solve( V)

u1<- t(Vi%*% t(x1))
u2<- t(Vi%*% t(x2))


look<- exp(-1*rdist(u1,u2))
look2<- stationary.cov( x1,x2, V= V)
test.for.zero( look, look2)

look<- Wendland(rdist(u1,u2), k=3, dimension=2)
look2<- stationary.cov( x1,x2, V= V, Covariance = "Wendland",
                        k=3, dimension=2)


test.for.zero( look, look2)


x1<- x[1:5,]
x2<- x[2:6,]
V<- matrix( c(2,1,0,4), 2,2)
Vi<- solve( V)

u1<- x1
u2<- x2

look1a<- exp(-1*rdist(u1,u2))
look1b<-  Wendland(rdist(u1,u2),
                   k=3, dimension=2, theta= 1)
look1<- look1a*look1b
look2<- stationary.taper.cov( x1,x2, theta=1,
                              Taper.args=list( theta=1,k=3, dimension=2), verbose=FALSE)
test.for.zero( look1, as.matrix(look2))


u1<- t(Vi%*% t(x1))
u2<- t(Vi%*% t(x2))


look1a<- exp(-1*rdist(u1,u2))
look1b<-  Wendland(rdist(u1,u2),
                   k=3, dimension=2, theta= 1.5)
look1<- look1a*look1b
look2<- stationary.taper.cov( x1,x2,V=V,
                              Taper.args=list( theta=1.5,k=3, dimension=2), verbose=FALSE)
test.for.zero( look1, as.matrix(look2))


u1<- t(Vi%*% t(x1))
u2<- t(Vi%*% t(x2))


look1a<- Matern(rdist(u1,u2), smoothness=1.5)
look1b<-  Wendland(rdist(u1,u2),
                   k=3, dimension=2, theta= 1.5)
look1<- look1a*look1b
look2<- stationary.taper.cov( x1,x2,V=V,Covariance=Matern, smoothness=1.5,
                              Taper.args=list( theta=1.5,k=3, dimension=2), verbose=FALSE)
test.for.zero( look1, as.matrix(look2))


# some tests of great circle distance


stationary.taper.cov( x[1:3,],x[1:10,] , theta=200, Taper.args= 
                        list(k=2,theta=300, dimension=2),
                      Dist.args=list( method="greatcircle") )-> temp

# temp is now a tapered 3X10 cross covariance matrix in sparse format. 
# should be identical to
# the direct matrix product

temp2<- Exponential( rdist.earth(x[1:3,],x[1:10,]), range=200) * 
  Wendland(rdist.earth(x[1:3,],x[1:10,]), theta= 300, k=2, dimension=2)

test.for.zero(  as.matrix(temp), temp2, tol=1e-6, tag="taper with great circle")

# example of calling the taper version directly 
# Note that default covariance is exponential and default taper is 
# Wendland (k=2).

stationary.taper.cov( x[1:3,],x[1:10,] , theta=1.5, Taper.args= 
                        list(k=2,theta=2.0, dimension=2) )-> temp
# temp is now a tapered 5X10 cross covariance matrix in sparse format. 
# should be identical to
# the direct matrix product

temp2<- Exp.cov( x[1:3,],x[1:10,], theta=1.5) * 
  Wendland(rdist(x[1:3,],x[1:10,]),
           theta= 2.0, k=2, dimension=2)

test.for.zero(  as.matrix(temp), temp2, tag= "high level test of taper cov")

stationary.taper.cov( x[1:3,],x[1:10,] , range=1.5,
                      Taper.args= list(k=2,theta=2.0,
                                       dimension=2) )-> temp

test.for.zero(  as.matrix(temp), temp2, tag= "high level test of taper cov")



##### Test precomputing distance matrix
#

y<- ozone2$y[16,]
x<- ozone2$lon.lat

#
# Omit the NAs

good<- !is.na( y)
x<- x[good,]
y<- y[good]

#####test that stationary.cov returns the same result when passed distance matrix:

#with x1 == x2:

x1<- x[1:20,]
compactDistMat = rdist(x1, compact=TRUE)
distMat = rdist(x1)
look<- stationary.cov(x1, theta=4)
look2 <- stationary.cov(x1, theta=4, distMat = compactDistMat)
look3 <- stationary.cov(x1, theta=4, distMat = distMat)
test.for.zero( look, look2, tag="stationary.cov versus stationary.cov compact distMat")
test.for.zero( look, look3, tag="stationary.cov versus stationary.cov matrix distMat")

#with x1 != x2:

x2=x[1:10,]
distMat = rdist(x1, x2)
look<- stationary.cov(x1, x2, theta=4)
look2 <- stationary.cov(x1, x2, theta=4, distMat = distMat)
test.for.zero( look, look2, tag="stationary.cov versus stationary.cov asymmetric distMat")

#####test that stationary.cov returns the same result when passed distance matrix:

#with x1 == x2:
distMat = rdist(x1, x1)
compactDistMat = rdist(x1, compact=TRUE)

look<- Exp.cov(x1, theta=4)
look2 <- Exp.cov(x1, theta=4, distMat = compactDistMat)
look3 <- Exp.cov(x1, theta=4, distMat = distMat)
test.for.zero( look, look2, tag="Exp.cov versus Exp.cov compact distMat")
test.for.zero( look, look3, tag="Exp.cov versus Exp.cov matrix distMat")

#with x1 != x2:

x1<- x[1:20,]
x2=x[1:10,]
distMat = rdist(x1, x2)
look<- Exp.cov(x1, x2, theta=4)
look2 <- Exp.cov(x1, x2, theta=4, distMat = distMat)
test.for.zero( look, look2, tag="Exp.cov versus Exp.cov asymmetric distMat")

##### test for correct value when using C argument:

Ctest<- rnorm(10)

#with x1 == x2:

x1 = x[1:10,]
compactDistMat = rdist(x1, compact=TRUE)
distMat = rdist(x1, x1)

temp1<- stationary.cov( x1, C= Ctest, theta=4 )
temp2 = stationary.cov( x1, C= Ctest, theta=4, distMat=compactDistMat )
temp3 = stationary.cov( x1, C= Ctest, theta=4, distMat=distMat )

exp1<- Exp.cov( x1, C= Ctest, theta=4 )
exp2 = Exp.cov( x1, C= Ctest, theta=4, distMat=compactDistMat )
exp3 = Exp.cov( x1, C= Ctest, theta=4, distMat=distMat )

test.for.zero(temp1, temp2, tag="stationary.cov vs stationary.cov with C set, compact distMat")
test.for.zero(temp1, temp3, tag="stationary.cov vs stationary.cov with C set, matrix distMat")
test.for.zero(temp1, exp1, tag="stationary.cov vs Exp.cov with C set, no distMat")
test.for.zero(temp2, exp2, tag="stationary.cov vs Exp.cov with C set, compact distMat")
test.for.zero(temp3, temp3, tag="stationary.cov vs Exp.cov with C set, matrix distMat")

#with x1 != x2:

x1 = x
x2 = x[1:10,]

distMat = rdist(x1, x1)

temp1<- stationary.cov( x1, x2, C= Ctest, theta=4 )
temp2 = stationary.cov( x1, x2, C= Ctest, theta=4, distMat=distMat )
exp1 <- Exp.cov( x1, x2, C= Ctest, theta=4 )
exp2 = Exp.cov( x1, x2, C= Ctest, theta=4, distMat=distMat )

test.for.zero(temp1, temp2, tag="stationary.cov vs stationary.cov with C set and asymmetric distMat given")
test.for.zero(exp1, exp2, tag="Exp.cov vs Exp.cov with C set and asymmetric distMat given")


##### test covariance functions for onlyUpper=TRUE
#

distMat = rdist(x1, x1)
compactDistMat = rdist(x1, compact=TRUE)
out1 = stationary.cov(x1, onlyUpper=TRUE)
exp1 = Exp.cov(x1, onlyUpper=TRUE)
out2 = stationary.cov(x1, onlyUpper=TRUE, distMat=compactDistMat)
exp2 = Exp.cov(x1, onlyUpper=TRUE, distMat=compactDistMat)
out3 = stationary.cov(x1, onlyUpper=TRUE, distMat=distMat)
exp3 = Exp.cov(x1, onlyUpper=TRUE, distMat=distMat)

test.for.zero( out2[upper.tri(out1)], out3[upper.tri(exp1)], tag="onlyUpper=TRUE: stationary.cov versus Exp.cov")
test.for.zero( out2[upper.tri(out1)], out3[upper.tri(out2)], tag="onlyUpper=TRUE: stationary.cov versus stationary.cov with compactDistMat")
test.for.zero( out2[upper.tri(out1)], out3[upper.tri(exp2)], tag="onlyUpper=TRUE: stationary.cov versus Exp.cov with compactDistMat")
test.for.zero( out2[upper.tri(out1)], out3[upper.tri(out3)], tag="onlyUpper=TRUE: stationary.cov versus stationary.cov with matrix distMat")
test.for.zero( out2[upper.tri(out1)], out3[upper.tri(exp3)], tag="onlyUpper=TRUE: stationary.cov versus Exp.cov with matrix distMat")

##### test Exp.cov functions for correct use of p
#

p1 = 1
p2 = 2
p3 = 3
distMat = rdist(x1, x1)

exp1 = Exp.cov(x1, p=p1)
exp2 = Exp.cov(x1, p=p2)
exp2Dist = Exp.cov(x1, p=p2, distMat = distMat)
exp3 = Exp.cov(x1, p=p3)
test.for.zero(exp1^(rdist(x1, x1)^(p2 - p1)), exp2, tag="Testing p=1 v 2")
test.for.zero(exp2^(rdist(x1, x1)^(p3 - p2)), exp3, tag="Testing p=2 v 3")
test.for.zero(exp2, exp2Dist, tag="Testing p=2 v 2 with distMat")

##### test Exp.cov functions for correct use of theta
#

theta1 = 1
theta2 = 2
theta3 = 3
distMat = rdist(x1, x1)

exp1 = Exp.cov(x1, theta=theta1)
exp2 = Exp.cov(x1, thet=theta2)
exp2Dist = Exp.cov(x1, theta=theta2, distMat = distMat)
exp3 = Exp.cov(x1, theta=theta3)
test.for.zero(exp1^(theta1/theta2), exp2, tag="Testing theta=1 v 2")
test.for.zero(exp2^(theta2/theta3), exp3, tag="Testing theta=2 v 3")
test.for.zero(exp2, exp2Dist, tag="Testing theta=2 v 2 with distMat")




cat("end tests of V argument in covariances", fill=TRUE)













