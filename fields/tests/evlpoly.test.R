# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( fields)
options( echo=FALSE)
test.for.zero.flag<-1

set.seed( 245)

x<- runif(3)

coef<- runif( 5)
temp<- fields.evlpoly( x, coef)

temp2<- coef[1]

for(  k in (2:5) ){
temp2<- temp2 + coef[k]*x**(k-1)
}

test.for.zero( temp, temp2)


set.seed( 124)
x<-  matrix( runif(12), ncol=3)

fields.mkpoly(x, m=3)-> out

attr( out, "ptab")-> ptab

J<- nrow( ptab)

coef<- runif( J)
temp<- fields.evlpoly2( x, coef, ptab)

temp2<-out%*% coef

test.for.zero( temp,temp2)

fields.derivative.poly( x, m=3, coef)-> temp

fields.mkpoly( cbind( x[,1:2], x[,3]+1e-6), m=3)%*% coef-> temp2
fields.mkpoly( cbind( x[,1:2], x[,3]-1e-6), m=3)%*% coef-> temp3
temp2<- (temp2- temp3)/ 2e-6

test.for.zero( temp[,3], temp2)

cat("Done testing polynomial evaluation",fill=TRUE)

options( echo=FALSE)




