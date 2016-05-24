# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# test of Wendland covariance and  stationary.taper.cov

library( fields)
options( echo=FALSE)
 test.for.zero.flag<- 1

set.seed(123)
x1<- matrix( runif(2*20), ncol=2)
x2<- matrix( runif(2*10), ncol=2)

fields.rdist.near( x1,x2, delta=.75)-> look

temp<- matrix( NA, nrow(x1),nrow(x2))
temp[ look$ind] <- look$ra

temp2<- rdist( x1, x2)
temp2[ temp2> .75] <- NA
#set.panel( 2,1) ; image.plot( temp); image.plot( temp2)

temp[ is.na( temp)]<- 0
temp2[ is.na( temp2)]<- 0
test.for.zero( temp, temp2)


# test of constructing covariance matrix
# and also versions of Wendland function
# default taper is wendland k=2.
DD<- rdist( x1,x2)
temp<- Wendland2.2(DD, theta=.8)
temp2<- Wendland( DD, theta=.8, k=2, dimension=2)

test.for.zero( temp, temp2)




stationary.taper.cov( x1,x2, Taper="Wendland2.2", 
           Taper.args= list( theta=.8), spam.format=FALSE )-> look
temp0<- look

stationary.taper.cov( x1,x2, Taper="Wendland2.2",
           Taper.args= list( theta=.8), spam.format=TRUE )-> look
temp1<-  spam2full( look)

test.for.zero( temp1, temp0)

stationary.taper.cov( x1,x2, Taper="Wendland",
           Taper.args= list( theta=.8, k=2, dimension=2),
                     spam.format=TRUE )-> look
temp1b<-  spam2full( look)

temp2<-  Wendland2.2(DD, theta=.8) * Exponential(DD)
temp3<-  wendland.cov(x1,x2, k=2, theta=.8) * Exponential(DD)
temp4<-  Wendland(DD, k=2, dimension=2, theta=.8)* Exponential(DD)


test.for.zero( temp1, temp0, rel=FALSE)
test.for.zero( temp1b, temp0, rel=FALSE)
test.for.zero( temp2, temp0, rel=FALSE)

test.for.zero( temp2, temp3,rel=FALSE)
test.for.zero( temp2, temp4,rel=FALSE)

set.seed( 256)
rv<- runif( nrow(x2))

# test of multiply 
stationary.taper.cov( x1, x2,  C= rv)-> look
temp2<-stationary.taper.cov( x1,x2)

(as.matrix(temp2))%*%(rv)-> look2
test.for.zero( look, look2)

temp2%*%(rv)-> look2
test.for.zero( look, look2)


cat( "Done with testing Wendland family", fill=TRUE)
options( echo=TRUE)
