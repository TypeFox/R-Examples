# some tests to ensure that modifications don't change
# the global behaviour of Time*DataFrame

#=====================
# TimeInstantDataFrame
library( timetools )

tidf <- TimeInstantDataFrame(as.POSIXct(sprintf('2012-01-%02i', 4:8), 'UTC'),
	'CET', data.frame(un=1:5, two=6:10))
tidf
timezone(tidf) <- 'UTC'
tidf
tidf[3:4,'two']
when(tidf)

tidf <- TimeInstantDataFrame(as.POSIXct(sprintf('2012-01-%02i', 4:8), 'CET'),
	'UTC', data.frame(un=1:5, two=6:10))
tidf
timezone(tidf) <- 'CET'
tidf
tidf[3:4,'two']
when(tidf)

split( tidf, c(rep('C', 3), rep('A', 2)) )

#======================
# TimeIntervalDataFrame

tidf <- TimeIntervalDataFrame(as.POSIXct(sprintf('2012-01-%02i', 4:9), 'UTC'),
	tiemzone='CET', data=data.frame(un=1:5, two=6:10))
tidf
timezone(tidf) <- 'UTC'
tidf
tidf[3:4,'two']
when(tidf)
start(tidf)
end(tidf)

tidf <- TimeIntervalDataFrame(as.POSIXct(sprintf('2012-01-%02i', 4:9), 'CET'),
	timezone='UTC', data=data.frame(un=1:5, two=6:10))
tidf
timezone(tidf) <- 'CET'
tidf
tidf[3:4,'two']
when(tidf)
start(tidf)
end(tidf)

