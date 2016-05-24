context('turn angle')
test_that('turnAngle',
{
	set.seed(3425)
	l<-cbind(c(59,60),c(59,60))
	r<-runif(1000, -180,180)
	d<-rexp(length(r), .0001)
	for( i in 1:length(r))
		l<-rbind(l,destPoint(tail(l,1), finalBearing(head(tail(l,2),1),tail(l,1))+r[i], d[i]))
	m<-move(l[,1],l[,2], as.POSIXct(1:nrow(l), origin='1970-1-1'), proj='+proj=longlat +ellps=WGS84')
	expect_equal(turnAngleGc(m),r)
})
