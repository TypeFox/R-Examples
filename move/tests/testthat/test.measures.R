context('measures')
test_that('measures',
{
	x<-move(x=(dx<-c(1,1:7)), y=c(1,1:7), as.POSIXct(c(1:7,10), origin='1970-1-1'))
	d<-c(0,rep(sqrt(2),6))
	t<-c(rep(1, 6), 3)
	expect_equal(distance(x),d)
	expect_equal(distanceSummary(x)$AverDist,mean(d))
	expect_equal(speed(x),d/t)
	expect_equal(speedSummary(x)$AverSpeed, mean(d/t))
	expect_equal(timeLag(x), t)
	expect_equal(timeLag(x, 'secs'), t)
	expect_equal(timeLag(x, 'mins'), t/60)
	expect_warning(xx<-moveStack(list(x,x)))
	expect_equal(distance(xx),list(unnamed=d,unnamed1=d))
	expect_warning(expect_equal(timeLag(xx),list(unnamed=t,unnamed1=t)),'Units not specified this could lead to different units for the time differences between individuals')
	expect_equal(speed(xx),list(unnamed=d/t,unnamed1=d/t))
	expect_warning(turnAngleGc(x),"turnAngleGc is probably not a valid calculation on this projection")
	xx<-x;proj4string(xx)<-"+proj=longlat +ellps=WGS84"
	
	expect_equal(distance(xx),(dd<-distGeo( cbind(dx,dx)[-n.locs(xx),], cbind(dx,dx)[-1,])))
	expect_equal(speed(xx),dd/t)
	expect_equal(length(turnAngleGc(xx)), n.locs(xx)-2)
}
)