library('mapmisc')


pdf("colourScaleFile1.pdf",height=12,width=12)
par(mfrow=c(3,3))

someData = SpatialPointsDataFrame(cbind(1:4, 1:4),data=data.frame(y=1:4))

cs=colourScale(someData$y,breaks=1:4,style='fixed',labels=c('a','b','c'))
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)

cs=colourScale(NA,breaks=1:4,style='fixed')
legendBreaks('topleft', cs)


cs=colourScale(NULL,breaks=1:4,style='fixed')
legendBreaks('topleft', cs)

cs=colourScale(breaks=1:4,style='fixed')
legendBreaks('topleft', cs)

cs=colourScale(someData$y,breaks=1:4,style='fixed',labels=c('a','b','c'),col=heat.colors(3), opacity=0.5)
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)

cs=colourScale(x=someData$y,breaks=1:4,
		style='unique',labels=c('a','b','c','d'),col=heat.colors(4), opacity=0.5)
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)
legend('bottomright', fill=cs$col, legend=cs$legend)

cs=colourScale(x=NA,breaks=1:4,style='unique',labels=c('a','b','c','d'),
		col=t(col2rgb(heat.colors(4))), opacity=0.5)
legendBreaks('topleft', cs)
legend('bottomright', fill=cs$col, legend=cs$legend)


cs=colourScale(x=NA,breaks=1:4,style='unique',labels=c('a','b','c','d'),col=t(col2rgb(heat.colors(4))), opacity=0.5,exclude=2)
legendBreaks('topleft', cs)
legend('bottomright', fill=cs$col, legend=cs$legend)

cs=colourScale(x=NA,breaks=1:4,style='unique',labels=c('a','b','c','d'),col=t(col2rgb(heat.colors(4))), opacity=0.5,exclude='a')
legendBreaks('topleft', cs)
legend('bottomright', fill=cs$col, legend=cs$legend)

someData = SpatialPointsDataFrame(cbind(
				sample(1:4, 12,replace=TRUE), 
				sample(1:4, 12,replace=TRUE)),
		data=data.frame(y=sample(0:5, 12,replace=TRUE)))
cs=colourScale(x=someData$y,
		breaks=1:4,style='unique',labels=c('a','b','c','d'),
		col=t(col2rgb(heat.colors(4))), opacity=0.5,exclude='a')
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)
forLegend = na.omit(cs$levels)
legend('bottomright', fill=forLegend$col, legend=forLegend$label)

cs=colourScale(x=someData$y,
		breaks=1:4,style='unique',labels=c('missing','a','b','c','d','e'),
		col=t(col2rgb(terrain.colors(4))), opacity=0.5,exclude='a')
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)
forLegend = na.omit(cs$levels)
legend('bottomright', fill=forLegend$col, legend=forLegend$label)


someData = SpatialPointsDataFrame(cbind(
				sample(1:4, 11,replace=TRUE), 
				sample(1:4, 11,replace=TRUE)),
		data=data.frame(y=0:10))
cs=colourScale(someData$y,breaks=4,style='equal',dec=2)
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)

cs=colourScale(someData$y,breaks=4,style='quantile',dec=2)
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)

cs=colourScale(x=someData$y,breaks=4,style='equal',exclude=0)
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)

cs=colourScale(x=someData$y,breaks=4,style='equal',exclude=c(0,10))
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)

dev.off()

pdf("colourScaleFile2.pdf",height=12,width=12)
par(mfrow=c(3,3))

cs=colourScale(x=someData$y,breaks=4,style='equal',exclude='nothing',dec=2)
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)

myraster = raster(matrix(c(0:7,2), 3, 3))
cs=colourScale(x=myraster,breaks=4,style='equal',dec=2)
plot(myraster, breaks=cs$breaks, col=cs$col,legend=FALSE)
legendBreaks('topright', cs)

cs=colourScale(x=myraster,breaks=4,style='quantile',dec=2)
plot(myraster, breaks=cs$breaks, col=cs$col,legend=FALSE)
legendBreaks('topright', cs)


cs=colourScale(x=myraster,breaks=4,style='unique')
plot(myraster, breaks=cs$breaks, col=cs$col,legend=FALSE)
legendBreaks('topright', cs)
forLegend = na.omit(cs$levels)
legend('bottomright', fill=forLegend$col, legend=forLegend$label)

cs=colourScale(x=myraster,breaks=1:4,style='unique',labels=c('a','b','c','d'))
plot(myraster, breaks=cs$breaks, col=cs$col,legend=FALSE)
legendBreaks('topright', cs)
forLegend = na.omit(cs$levels)
legend('bottomright', fill=forLegend$col, legend=forLegend$label)

cs=colourScale(x=myraster,breaks=1:4,style='fixed')
plot(myraster, breaks=cs$breaks, col=cs$col,legend=FALSE)
legendBreaks('topright', cs)



cs=colourScale(x=myraster,breaks=4,style='equal',exclude=0)
plot(myraster, breaks=cs$breaks, col=cs$col,legend=FALSE)
legendBreaks('topright', cs)

cs=colourScale(x=myraster,breaks=4,style='unique',exclude=0)
legendBreaks('topright', cs)
forLegend = na.omit(cs$levels)
legend('bottomright', fill=forLegend$col, legend=forLegend$label)

someData = SpatialPointsDataFrame(cbind(
				sample(1:4, 10,replace=TRUE), 
				sample(1:4, 10,replace=TRUE)),
		data=data.frame(y=factor(sample(0:4,10,replace=TRUE))))

cs=colourScale(x=someData$y,breaks=4,
		style='thisShouldBeIgnored',exclude='0')
plot(someData, pch=16, col=cs$plot)
legendBreaks('topleft', cs)
forLegend = na.omit(cs$levels)
legend('bottomright', fill=forLegend$col, legend=forLegend$label)

dev.off()

