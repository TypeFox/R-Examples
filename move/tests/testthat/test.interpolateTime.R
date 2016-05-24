context('interpolate time')
test_that('interpolateTime',{
data(leroy)
  data<-leroy
dataP<-spTransform(data, center=T)
	expect_equal(data, interpolateTime(data, timestamps(data), spaceMethod='gr'))
	expect_equal(data, interpolateTime(data, timestamps(data), spaceMethod='rh'))
	expect_equal(dataP, interpolateTime(dataP, timestamps(dataP), spaceMethod='eu'))
	sl<-c(levels(data@sensor),'interpolateTime')
	dataS<-new('Move', data,sensor=factor(data@sensor, levels=sl), sensorUnUsedRecords=factor(data@sensorUnUsedRecords, levels=sl))
	rownames(dataS@coords)<-NULL
  expect_equal(dataS[c(1,n.locs(data)),], interpolateTime(data, 30, spaceMethod='gr')[c(1,30),])
	expect_equal(dataS[c(1,n.locs(data)),], interpolateTime(data, 30, spaceMethod='rh')[c(1,30),])
	expect_equal(spTransform(dataS[c(1,n.locs(data)),], proj4string(dataP)), interpolateTime(dataP, 30, spaceMethod='eu')[c(1,30),])
	expect_equal(40, n.locs(interpolateTime(data, 40, spaceMethod='g')))
	t<-move(0:1,0:1, as.POSIXct(0:1, origin='1970-1-1'))
	tt<-move(0:1,0:1, as.POSIXct(c(0,10), origin='1970-1-1'))
	expect_equivalent(coordinates(interpolateTime(t, as.POSIXct(.1,origin='1970-1-1'))), cbind(.1,.1))
	expect_equal(coordinates(interpolateTime(t, as.POSIXct(.16,origin='1970-1-1'))), coordinates(interpolateTime(tt, as.POSIXct(1.6,origin='1970-1-1'))))
	# prevent lineMidpoint from failing
	d<-data[1:40,]
	crd<-do.call(rbind,lapply(lapply(split(burst(d, 1:39)),move:::lineMidpoint), coordinates))
	crd2<-coordinates(interpolateTime(d, timestamps(d)[-40]+ (timestamps(d)[-1]-timestamps(d)[-40])/2, spaceMethod='g'))
	
	expect_equal(c(crd), c(crd2))
expect_warning(a<-interpolateTime(data, as.difftime(1,units='days')),'Euclidean interpolation seems unsuitable for the longitude latitude projection')
expect_equal(unique(timeLag(a, units='days')),1)
expect_equal(floor(sum(timeLag(data,'days'))),
             sum(timeLag(a,'days')))
}
)