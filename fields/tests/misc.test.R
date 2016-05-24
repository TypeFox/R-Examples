# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( fields)
options( echo=FALSE)
set.seed( 234)
test.for.zero.flag<-1

y<- runif( 30)
lev<- sort(sample( 1:5,30, replace=TRUE))
w<- runif( 30)*.1+1
y<- as.matrix(y)

# compute by loop
hold<- rep( NA, 5)
for( k in 1:5){
  ind<- lev==k
  hold[k]<- sum( y[ind,]*w[ind])/ sum( w[ind])}

look<- fast.1way( lev, y, w)
test.for.zero( look$means, hold, tag="fast.1way means")

# now vectorized case

ytemp<- cbind( y, y-10, y+10)
look2<- fast.1way( lev, ytemp, w)
test.for.zero( look2$means[,2], hold-10, tag="fast.1way vectorized means")


cat("All done with testing misc functions", fill=TRUE)
options(echo=TRUE)
