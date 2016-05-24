# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( fields)
options( echo=FALSE)
test.for.zero.flag<- 1

data(COmonthlyMet)
y<- CO.tmin.MAM.climate
good<- !is.na( y)
y<-y[good]
x<- CO.loc[good,]
Z<- CO.elev[good]
out<- mKrig( x,y, Z=Z,  cov.function="stationary.cov", Covariance="Matern",
                    theta=4.0, smoothness=1.0, lambda=.1)

out2<- Krig( x,y, Z=Z,  cov.function="stationary.cov", Covariance="Matern",
                    theta=4.0, smoothness=1.0, lambda=.1, GCV=TRUE)

test.for.zero( predict( out), predict(out2), tag="Full prediction")
test.for.zero( predict( out, drop.Z=TRUE), predict(out2, drop.Z=TRUE), tag=" prediction dropping Z")

xnew<- CO.loc[!good,]
Znew<-  CO.elev[!good]
temp1<- predict( out, xnew=xnew, drop.Z=TRUE)
temp2<- predict( out2, x=xnew, drop.Z=TRUE)
test.for.zero( temp1,temp2, tag="new x's dropping Z")

temp1<- predict( out, xnew=xnew, Z=Znew)
temp2<- predict( out2, x=xnew, Z=Znew)
test.for.zero( temp1,temp2, tag="new x's new Z's")

temp1<- predictSurface( out, nx=20, ny=20, drop.Z=TRUE, extrap=TRUE)
temp2<- predictSurface( out2, nx=20, ny=20, drop.Z=TRUE, extrap=TRUE)
test.for.zero( temp1$z,temp2$z, tag="predicting on surface with drop.Z")


cat("all done with mKrig Z tests", fill=TRUE)
options( echo=TRUE)

