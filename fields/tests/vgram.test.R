# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# test of vgram

library( fields)
options(echo=FALSE)

data( ozone2)

y<- ozone2$y[16,]
x<- ozone2$lon.lat

vgram( x, y, lon.lat=TRUE)-> out

# compute "by hand"

outer( y, y ,"-")-> hold
hold<- .5*hold^2 
rdist.earth( x,x)-> hold2
col( hold2)> row( hold2)-> upper

hold<- hold[upper]
hold2<- hold2[upper]
 order( hold2)-> od
hold2<- hold2[od]
hold<- hold[od]
ind<- is.na(hold)
hold<- hold[!ind]
hold2<- hold2[!ind]

test.for.zero( hold, out$vgram, tag="vgram single time")


# multiple times including NAs at some times

y<- t(ozone2$y[16:18,])
x<- ozone2$lon.lat[,]

out<- vgram( x, y, lon.lat=TRUE)


N<- nrow( y)

hold<-  cbind(c(outer( y[,1], y[,1],"-")),
         c(outer( y[,2], y[,2],"-") ),
         c(outer(y[,3], y[,3],"-"))  )
hold<- .5*hold^2
hold<- rowMeans( hold, na.rm=TRUE)
hold<- matrix( hold, N,N)

rdist.earth( x,x)-> hold2

col( hold2)> row( hold2)-> upper
hold<- hold[upper]
hold2<- hold2[upper]

order( hold2)-> od
hold2<- hold2[od]
hold<- hold[od]

ind<- is.na(hold)
hold<- hold[!ind]
hold2<- hold2[!ind]

test.for.zero( hold, out$vgram, tag="vgram more than one time point")

# test covariogram versus correlogram
y<- ozone2$y[16,]
x<- ozone2$lon.lat

sigma2 = var(y, na.rm=TRUE)
lookCov = vgram(x, y, lon.lat=TRUE, type="covariogram")
lookCor = vgram(x, y, lon.lat=TRUE, type="correlogram")

test.for.zero(lookCov$vgram*(1/sigma2), lookCor$vgram, tag="correlogram versus covariogram")

# test cross-covariogram versus cross-correlogram

sigma2 = var(y, na.rm=TRUE)
lookCov = crossCoVGram(x, x, y, y, lon.lat=TRUE, type="cross-covariogram")
lookCor = crossCoVGram(x, x, y, y, lon.lat=TRUE, type="cross-correlogram")

test.for.zero(lookCov$vgram*(1/sigma2), lookCor$vgram, tag="correlogram versus covariogram")